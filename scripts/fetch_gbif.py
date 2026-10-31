"""Download GBIF occurrence records for a species or genus.

Two modes:

  search     the occurrence/search API, no account needed, capped at 100,000
             records. Convenient while developing.
  download   the download API, which needs a free GBIF account and returns a
             citable DOI for a frozen snapshot. Use this for the version of the
             data that goes into the manuscript.

Usage:
    python fetch_gbif.py "Tetranychus urticae" gbif_tetranychus_urticae
    python fetch_gbif.py "Tetranychus" gbif_tetranychus_genus --download

For --download, set GBIF_USER, GBIF_PASSWORD and GBIF_EMAIL in the environment.
The DOI is written to data/raw/<name>_download.txt; cite it in the manuscript.

Writes data/raw/<name>.csv and data/raw/<name>_query.txt
"""
import base64
import csv
import io
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile

API = "https://api.gbif.org/v1"
RAW = os.path.join("data", "raw")
FIELDS = ["gbifID", "scientificName", "acceptedScientificName", "decimalLatitude",
          "decimalLongitude", "coordinateUncertaintyInMeters", "country", "countryCode",
          "stateProvince", "locality", "year", "month", "eventDate", "basisOfRecord",
          "institutionCode", "catalogNumber", "datasetName", "recordedBy",
          "identifiedBy", "associatedTaxa", "habitat", "issues"]
POLL_SECONDS = 20
POLL_LIMIT = 90


def get_json(url, auth=None):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    if auth:
        req.add_header("Authorization", "Basic " + auth)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as err:
            print(f"Retry {attempt + 1}: {err}")
            time.sleep(10)
    sys.exit(f"Failed: {url}")


def taxon_key(taxon):
    match = get_json(f"{API}/species/match?name={urllib.parse.quote(taxon)}")
    key = match.get("usageKey")
    if not key:
        sys.exit(f"GBIF found no taxon for {taxon}: {match}")
    print(f"GBIF taxon: {match.get('scientificName')} "
          f"(key {key}, match {match.get('matchType')}, rank {match.get('rank')})")
    return key


def write_query(name, text):
    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, f"{name}_query.txt"), "w") as f:
        f.write(text + "\n")


def via_search(key, name, taxon):
    path = os.path.join(RAW, f"{name}.csv")
    os.makedirs(RAW, exist_ok=True)
    offset, limit, total = 0, 300, 0
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        while True:
            page = get_json(f"{API}/occurrence/search?taxonKey={key}"
                            f"&occurrenceStatus=PRESENT&limit={limit}&offset={offset}")
            for rec in page.get("results", []):
                row = {k: rec.get(k, "") for k in FIELDS}
                row["issues"] = ";".join(rec.get("issues", []))
                w.writerow(row)
                total += 1
            print(f"{total}/{page.get('count')}")
            if page.get("endOfRecords") or offset + limit >= 100000:
                break
            offset += limit
    write_query(name, f"mode\tsearch\ntaxon\t{taxon}\ntaxonKey\t{key}\n"
                      f"records\t{total}\nretrieved\t{time.strftime('%Y-%m-%d')}\n"
                      f"note\tno citable DOI; re-run with --download for the manuscript")
    print(f"Saved {path} ({total} records)")


def via_download(key, name, taxon):
    user = os.environ.get("GBIF_USER")
    password = os.environ.get("GBIF_PASSWORD")
    email = os.environ.get("GBIF_EMAIL")
    if not (user and password and email):
        sys.exit("--download needs GBIF_USER, GBIF_PASSWORD and GBIF_EMAIL "
                 "in the environment. A free account at gbif.org provides them.")
    auth = base64.b64encode(f"{user}:{password}".encode()).decode()

    body = json.dumps({
        "creator": user,
        "notificationAddresses": [email],
        "sendNotification": False,
        "format": "SIMPLE_CSV",
        "predicate": {"type": "and", "predicates": [
            {"type": "equals", "key": "TAXON_KEY", "value": str(key)},
            {"type": "equals", "key": "OCCURRENCE_STATUS", "value": "present"},
        ]},
    }).encode()
    req = urllib.request.Request(f"{API}/occurrence/download/request", data=body,
                                 method="POST")
    req.add_header("Content-Type", "application/json")
    req.add_header("Authorization", "Basic " + auth)
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            download_key = r.read().decode("utf-8").strip()
    except urllib.error.HTTPError as err:
        sys.exit(f"GBIF refused the download request ({err.code}): "
                 f"{err.read().decode('utf-8', 'replace')[:400]}")
    print(f"Download key {download_key} — GBIF is preparing the archive")

    meta = {}
    for attempt in range(POLL_LIMIT):
        meta = get_json(f"{API}/occurrence/download/{download_key}", auth=auth)
        status = meta.get("status")
        if status in ("SUCCEEDED", "KILLED", "FAILED", "CANCELLED"):
            break
        print(f"  {status} ({(attempt + 1) * POLL_SECONDS}s)")
        time.sleep(POLL_SECONDS)
    if meta.get("status") != "SUCCEEDED":
        sys.exit(f"Download ended as {meta.get('status')}")

    doi = meta.get("doi", "")
    link = meta.get("downloadLink")
    print(f"Ready: DOI {doi}, {meta.get('totalRecords')} records")

    with urllib.request.urlopen(link, timeout=600) as r:
        blob = r.read()
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, f"{name}.csv")
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        member = next((m for m in z.namelist() if m.endswith(".csv")), None)
        if not member:
            sys.exit("No CSV inside the GBIF archive")
        raw = z.read(member).decode("utf-8", "replace")

    reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
    total = 0
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for rec in reader:
            row = {k: rec.get(k, "") for k in FIELDS}
            row["issues"] = rec.get("issue", rec.get("issues", ""))
            # SIMPLE_CSV carries countryCode but not country.
            if not row.get("country"):
                row["country"] = rec.get("countryCode", "")
            w.writerow(row)
            total += 1

    write_query(name, f"mode\tdownload\ntaxon\t{taxon}\ntaxonKey\t{key}\n"
                      f"downloadKey\t{download_key}\ndoi\t{doi}\n"
                      f"records\t{total}\nretrieved\t{time.strftime('%Y-%m-%d')}\n")
    with open(os.path.join(RAW, f"{name}_download.txt"), "w") as f:
        f.write(f"GBIF Occurrence Download https://doi.org/{doi}\n"
                f"accessed {time.strftime('%Y-%m-%d')}\n")
    print(f"Saved {path} ({total} records); DOI recorded in "
          f"{name}_download.txt — cite it in the manuscript")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 2:
        sys.exit('Usage: python fetch_gbif.py "<taxon>" <output_name> [--download]')
    taxon, name = args
    key = taxon_key(taxon)
    if "--download" in sys.argv:
        via_download(key, name, taxon)
    else:
        via_search(key, name, taxon)


if __name__ == "__main__":
    main()
