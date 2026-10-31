"""Fetch the per-record identification verdicts from the preceding study.

The reliability weights used by deficit.py come from the archived code and
results of the reference-library reliability study, deposited at
10.5281/zenodo.23808136. Pulling them from the frozen deposit rather than from
a live repository keeps this analysis reproducible and gives the manuscript a
citable source.

Usage:
    python fetch_zenodo.py [record_id_or_doi]

Writes the wanted tables to data/external/ and the deposit's version metadata
to data/external/zenodo_source.txt
"""
import io
import json
import os
import re
import sys
import urllib.request
import zipfile

DEFAULT = "10.5281/zenodo.23808136"
API = "https://zenodo.org/api/records/{rid}"

# Tables deficit.py and holdout.py look for, by basename.
WANTED = {
    "record_anchor_check.csv",
    "anchor_status.csv",
    "record_checks.csv",
    "cluster_anchors.csv",
}
OUT = os.path.join("data", "external")


def record_id(token):
    m = re.search(r"zenodo\.(\d+)", token) or re.fullmatch(r"(\d+)", token.strip())
    if not m:
        sys.exit(f"Cannot read a Zenodo record id out of {token!r}")
    return m.group(1)


def get(url, binary=False):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    return data if binary else json.loads(data.decode("utf-8"))


def save(name, blob, found):
    path = os.path.join(OUT, name)
    with open(path, "wb") as f:
        f.write(blob)
    found.add(name)
    print(f"  wrote {path} ({len(blob)} bytes)")


def main():
    token = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    rid = record_id(token)
    os.makedirs(OUT, exist_ok=True)

    meta = get(API.format(rid=rid))
    title = meta.get("metadata", {}).get("title", "")
    version = meta.get("metadata", {}).get("version", "")
    doi = meta.get("doi", "")
    print(f"Zenodo record {rid}: {title}")
    print(f"  DOI {doi}  version {version or 'unversioned'}")

    found = set()
    for entry in meta.get("files", []):
        name = entry.get("key", "")
        link = (entry.get("links") or {}).get("self")
        if not link:
            continue
        if name in WANTED:
            save(name, get(link, binary=True), found)
            continue
        if name.lower().endswith(".zip"):
            print(f"  opening {name}")
            blob = get(link, binary=True)
            with zipfile.ZipFile(io.BytesIO(blob)) as z:
                for member in z.namelist():
                    base = os.path.basename(member)
                    if base in WANTED and base not in found:
                        save(base, z.read(member), found)

    missing = WANTED - found
    with open(os.path.join(OUT, "zenodo_source.txt"), "w") as f:
        f.write(f"record_id\t{rid}\ndoi\t{doi}\ntitle\t{title}\n"
                f"version\t{version}\nfetched\t{', '.join(sorted(found))}\n")

    if not found:
        sys.exit("No wanted tables found in the deposit. Check the record id.")
    if missing:
        print(f"Not in the deposit: {', '.join(sorted(missing))}")
    print(f"{len(found)} table(s) in {OUT}/")


if __name__ == "__main__":
    main()
