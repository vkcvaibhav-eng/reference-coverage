"""Attach WorldClim 2.1 bioclimatic values to every georeferenced record.

Walks data/raw/ recursively for GenBank (.gb), BOLD (.tsv/.jsonl) and GBIF (.csv)
files, pulls out coordinates and collection year, and samples BIO1-BIO19 at each
point. Rasters are downloaded to a temporary folder and discarded; only the
extracted table is written.

Usage:
    python climate_extract.py [10m|5m|2.5m] [--koppen PATH]

--koppen takes a local Köppen-Geiger classification raster (for example the
Beck et al. 1 km product) and adds its class code as a column. Download it
yourself and pass the path; it is not fetched automatically.

Writes data/derived/climate_worldclim_<resolution>.csv
"""
import csv
import glob
import json
import os
import re
import sys
import tempfile
import urllib.request
import zipfile

from Bio import SeqIO

RAW = os.path.join("data", "raw")
OUT = os.path.join("data", "derived")
URL = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_{res}_bio.zip"
NUM = re.compile(r"-?\d+(?:\.\d+)?")
YEAR = re.compile(r"(1[89]\d{2}|20\d{2})")
FIELDS = ["source", "record_id", "taxon", "country", "host", "year", "marker",
          "lat", "lon"]


def year_from(text):
    m = YEAR.search(str(text or ""))
    return m.group(1) if m else ""


def marker_label(text):
    """Normalise a gene or marker string to a short code.

    Only the barcode region of COI can serve as reference coverage for a COI
    library, so COI-5P is kept distinct from COI-3P and from other genes.
    """
    t = str(text or "").upper().replace("_", "-")
    if not t:
        return ""
    if any(k in t for k in ("COI-3P", "COX1-3P", "COXI-3P")):
        return "COI-3P"
    if any(k in t for k in ("CYTB", "CYTOCHROME B", "CYTOCHROME-B")):
        return "CYTB"
    if any(k in t for k in ("COI", "COXI", "COX1", "CO1", "CYTOCHROME C OXIDASE",
                            "CYTOCHROME OXIDASE")):
        return "COI-5P"
    for key in ("16S", "18S", "28S", "ITS2", "ITS1", "ITS"):
        if key in t:
            return key
    return t.split()[0][:12]


def genbank_points():
    for path in glob.glob(os.path.join(RAW, "**", "*.gb"), recursive=True):
        for rec in SeqIO.parse(path, "genbank"):
            src = next((f for f in rec.features if f.type == "source"), None)
            if not src or "lat_lon" not in src.qualifiers:
                continue
            m = re.match(r"([\d.]+)\s*([NS])\s+([\d.]+)\s*([EW])",
                         src.qualifiers["lat_lon"][0])
            if not m:
                continue
            q = src.qualifiers
            genes = {g for f in rec.features if f.type in ("gene", "CDS")
                     for g in f.qualifiers.get("gene", [])}
            products = {p for f in rec.features if f.type == "CDS"
                        for p in f.qualifiers.get("product", [])}
            yield {"source": "GenBank", "record_id": rec.id,
                   "taxon": q.get("organism", [""])[0],
                   "country": q.get("country", q.get("geo_loc_name", [""]))[0],
                   "host": q.get("host", [""])[0],
                   "year": year_from(q.get("collection_date", [""])[0]),
                   "marker": marker_label(" ".join(genes | products)),
                   "lat": float(m.group(1)) * (1 if m.group(2) == "N" else -1),
                   "lon": float(m.group(3)) * (1 if m.group(4) == "E" else -1)}


def bold_rows():
    for path in glob.glob(os.path.join(RAW, "**", "bold_*.tsv"), recursive=True):
        with open(path) as f:
            yield from csv.DictReader(f, delimiter="\t")
    for path in glob.glob(os.path.join(RAW, "**", "bold_*.jsonl"), recursive=True):
        with open(path) as f:
            text = f.read().strip()
        if text.startswith("["):
            yield from json.loads(text)
        else:
            yield from (json.loads(l) for l in text.splitlines() if l.strip())


def bold_points():
    for r in bold_rows():
        # Centroid coordinates are a country or province centre, not a site.
        if "centroid" in (r.get("coord_source") or "").lower():
            continue
        coord = r.get("coord")
        if isinstance(coord, list) and len(coord) == 2:
            lat, lon = coord
        else:
            nums = NUM.findall(str(coord or ""))
            if len(nums) != 2:
                continue
            lat, lon = nums
        sequence = (r.get("nuc") or r.get("nucraw") or "").replace("-", "").strip()
        yield {"source": "BOLD", "record_id": r.get("processid", ""),
               "taxon": r.get("species", ""), "country": r.get("country/ocean", ""),
               "host": r.get("associated_taxa", ""),
               "year": year_from(r.get("collection_date")
                                 or r.get("collection_date_start") or ""),
               "marker": marker_label(r.get("marker_code")) if sequence else "",
               "lat": float(lat), "lon": float(lon)}


def gbif_points():
    for path in glob.glob(os.path.join(RAW, "**", "gbif_*.csv"), recursive=True):
        if path.endswith("_query.txt"):
            continue
        with open(path) as f:
            for r in csv.DictReader(f):
                if not (r.get("decimalLatitude") and r.get("decimalLongitude")):
                    continue
                yield {"source": "GBIF", "record_id": r.get("gbifID", ""),
                       "taxon": r.get("scientificName", ""),
                       "country": r.get("country", ""),
                       "host": r.get("associatedTaxa", ""),
                       "year": r.get("year") or year_from(r.get("eventDate")),
                       "marker": "",
                       "lat": float(r["decimalLatitude"]),
                       "lon": float(r["decimalLongitude"])}


def sample_raster(path, coords, points, column, digits=3):
    import rasterio  # imported here so the module loads without the geospatial stack

    with rasterio.open(path) as src:
        nodata = src.nodata
        for p, val in zip(points, src.sample(coords)):
            v = float(val[0])
            blank = nodata is not None and v == nodata
            p[column] = "" if blank else (round(v, digits) if digits else int(v))


def main():
    res = next((a for a in sys.argv[1:] if not a.startswith("--")), "10m")
    koppen = None
    if "--koppen" in sys.argv:
        i = sys.argv.index("--koppen")
        if i + 1 >= len(sys.argv):
            sys.exit("--koppen needs a path to a classification raster")
        koppen = sys.argv[i + 1]
        if not os.path.exists(koppen):
            sys.exit(f"No such raster: {koppen}")

    points = list(genbank_points()) + list(bold_points()) + list(gbif_points())
    print(f"{len(points)} records with coordinates")
    if not points:
        sys.exit(f"Nothing found under {RAW}/ - run the fetch scripts first.")
    with_year = sum(1 for p in points if p["year"])
    print(f"{with_year} of them carry a collection year "
          f"({with_year / len(points):.0%})")
    markers = {}
    for p in points:
        if p["source"] != "GBIF":
            markers[p["marker"] or "(no sequence)"] = \
                markers.get(p["marker"] or "(no sequence)", 0) + 1
    if markers:
        print("sequence-record markers: "
              + ", ".join(f"{k} {v}" for k, v in
                          sorted(markers.items(), key=lambda kv: -kv[1])))

    tmp = tempfile.mkdtemp()
    url = URL.format(res=res)
    zpath = os.path.join(tmp, "wc.zip")
    print("Downloading", url)
    urllib.request.urlretrieve(url, zpath)
    with zipfile.ZipFile(zpath) as z:
        z.extractall(tmp)

    coords = [(p["lon"], p["lat"]) for p in points]
    for i in range(1, 20):
        hits = glob.glob(os.path.join(tmp, "**", f"wc2.1_{res}_bio_{i}.tif"),
                         recursive=True)
        if not hits:
            sys.exit(f"BIO{i} missing from the WorldClim archive for {res}")
        sample_raster(hits[0], coords, points, f"bio{i}")

    columns = FIELDS + [f"bio{i}" for i in range(1, 20)]
    if koppen:
        print("Sampling", koppen)
        sample_raster(koppen, coords, points, "koppen", digits=0)
        columns.append("koppen")

    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, f"climate_worldclim_{res}.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(points)
    print("Saved", out)


if __name__ == "__main__":
    main()
