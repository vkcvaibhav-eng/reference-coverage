"""Georeference validation for GBIF occurrence records.

Applies GBIF's own issue flags plus local geometric tests, and writes a cleaned
file alongside a report of how many records each rule removed. Records flagged
as photograph-based observations are separated rather than deleted: a spider
mite cannot be identified to species from a photograph, so they belong in a
sensitivity analysis, not in the main occurrence set.

Usage:
    python clean_coords.py [--keep-observations] [--max-uncertainty METRES]

Reads  data/raw/gbif_*.csv
Writes data/raw/<name>_clean.csv and results/cleaning/<name>_report.txt
"""
import csv
import glob
import os
import sys

# GBIF issue flags that make a coordinate unusable for climate extraction.
FATAL_ISSUES = {
    "ZERO_COORDINATE",
    "COORDINATE_OUT_OF_RANGE",
    "COORDINATE_INVALID",
    "COUNTRY_COORDINATE_MISMATCH",
    "PRESUMED_SWAPPED_COORDINATE",
    "PRESUMED_NEGATED_LATITUDE",
    "PRESUMED_NEGATED_LONGITUDE",
}
# Flags that mark a coordinate as derived from an administrative centre.
CENTROID_ISSUES = {
    "COUNTRY_COORDINATE_MISMATCH",
    "COORDINATE_ROUNDED",
}
PHOTO_BASIS = {"HUMAN_OBSERVATION", "MACHINE_OBSERVATION"}
DEFAULT_MAX_UNCERTAINTY = 50000.0


def flt(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    keep_obs = "--keep-observations" in sys.argv
    max_unc = DEFAULT_MAX_UNCERTAINTY
    if "--max-uncertainty" in sys.argv:
        i = sys.argv.index("--max-uncertainty")
        max_unc = float(sys.argv[i + 1])

    paths = [p for p in glob.glob(os.path.join("data", "raw", "gbif_*.csv"))
             if not p.endswith("_clean.csv")]
    if not paths:
        sys.exit("No data/raw/gbif_*.csv found. Run fetch_gbif.py first.")

    os.makedirs(os.path.join("results", "cleaning"), exist_ok=True)
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        with open(path) as f:
            rows = list(csv.DictReader(f))
            columns = rows[0].keys() if rows else []
        if not rows:
            print(f"{path} is empty, skipping")
            continue

        counts = {k: 0 for k in ("no_coordinate", "fatal_issue", "centroid_like",
                                 "equal_lat_lon", "integer_coordinate",
                                 "uncertain", "observation", "kept")}
        kept, observations = [], []
        for r in rows:
            lat, lon = flt(r.get("decimalLatitude")), flt(r.get("decimalLongitude"))
            issues = set((r.get("issues") or "").split(";")) - {""}

            if lat is None or lon is None:
                counts["no_coordinate"] += 1
                continue
            if issues & FATAL_ISSUES:
                counts["fatal_issue"] += 1
                continue
            if lat == 0.0 and lon == 0.0:
                counts["no_coordinate"] += 1
                continue
            if lat == lon:
                counts["equal_lat_lon"] += 1
                continue
            if issues & CENTROID_ISSUES:
                counts["centroid_like"] += 1
                continue
            if lat == int(lat) and lon == int(lon):
                counts["integer_coordinate"] += 1
                continue
            unc = flt(r.get("coordinateUncertaintyInMeters"))
            if unc is not None and unc > max_unc:
                counts["uncertain"] += 1
                continue
            if r.get("basisOfRecord") in PHOTO_BASIS:
                counts["observation"] += 1
                observations.append(r)
                if not keep_obs:
                    continue
            kept.append(r)
        counts["kept"] = len(kept)

        out = os.path.join("data", "raw", f"{name}_clean.csv")
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(columns))
            w.writeheader()
            w.writerows(kept)

        obs_out = os.path.join("data", "raw", f"{name}_observations.csv")
        with open(obs_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(columns))
            w.writeheader()
            w.writerows(observations)

        report = [f"source file                {path}",
                  f"records in                 {len(rows)}",
                  f"maximum uncertainty        {max_unc:g} m",
                  f"photograph records         "
                  f"{'kept in main set' if keep_obs else 'held back'}",
                  ""]
        report += [f"removed, {k:<24} {v}" for k, v in counts.items()
                   if k not in ("kept", "observation")]
        report += [f"photograph-based           {counts['observation']}",
                   f"records out                {counts['kept']}",
                   f"cleaned file               {out}",
                   f"photograph file            {obs_out}"]
        text = "\n".join(report)
        with open(os.path.join("results", "cleaning", f"{name}_report.txt"), "w") as f:
            f.write(text + "\n")
        print(text)
        print()


if __name__ == "__main__":
    main()
