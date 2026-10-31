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
    "merged_records.csv",      # the accession <-> process-ID crosswalk
}
# Files the analysis cannot run without.
REQUIRED = {"record_anchor_check.csv", "merged_records.csv"}
OUT = os.path.join("data", "external")

# The sequence-side raw data, taken from the same frozen deposit so that the
# reference library this analysis measures is the one the verdicts describe.
# Matched by extension and filename prefix rather than exact name.
RAW_OUT = os.path.join("data", "raw")
RAW_SUFFIXES = (".gb", ".tsv", ".jsonl")
RAW_PREFIXES = ("ncbi_", "bold_")


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


def is_raw_sequence(name):
    base = os.path.basename(name)
    return (base.endswith(RAW_SUFFIXES)
            and base.startswith(RAW_PREFIXES)
            and "query" not in base)


def flatten(member):
    """Full relative path as one filename, so no two members can collide.

    Prefixing only the immediate parent is not enough: every genus files its
    anchor table under a folder called "anchors", so results/anchors/x.csv and
    results/oligonychus/anchors/x.csv both reduce to anchors__x.csv and one
    silently overwrites the other. The whole path is kept instead.
    """
    parts = [p for p in member.replace("\\", "/").split("/")
             if p not in ("", ".", "..")]
    while parts and parts[0] in ("results", "data"):
        parts = parts[1:]
    return "__".join(parts) if parts else os.path.basename(member)


def save(name, blob, found, directory=None):
    directory = directory or OUT
    stem = flatten(name)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, stem)
    if stem in found:
        sys.exit(f"Name collision writing {name} as {stem}; refusing to "
                 f"overwrite. Report this as a bug.")
    with open(path, "wb") as f:
        f.write(blob)
    found.add(stem)
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

    found, raw_found = set(), set()
    for entry in meta.get("files", []):
        name = entry.get("key", "")
        link = (entry.get("links") or {}).get("self")
        if not link:
            continue
        if os.path.basename(name) in WANTED:
            save(name, get(link, binary=True), found)
            continue
        if is_raw_sequence(name):
            save(name, get(link, binary=True), raw_found, RAW_OUT)
            continue
        if name.lower().endswith(".zip"):
            print(f"  opening {name}")
            blob = get(link, binary=True)
            with zipfile.ZipFile(io.BytesIO(blob)) as z:
                for member in z.namelist():
                    if member.endswith("/"):
                        continue
                    base = os.path.basename(member)
                    if base in WANTED:
                        save(member, z.read(member), found)
                    elif is_raw_sequence(member) and base not in raw_found:
                        save(base, z.read(member), raw_found, RAW_OUT)

    if raw_found:
        print(f"{len(raw_found)} sequence file(s) in {RAW_OUT}/: "
              f"{', '.join(sorted(raw_found))}")
    else:
        print(f"No GenBank or BOLD files in the deposit. climate_extract.py will "
              f"find no sequence records unless you put them in {RAW_OUT}/ "
              f"yourself.")

    missing = {w for w in WANTED
               if not any(f == w or f.endswith("__" + w) for f in found)}
    absent_required = {w for w in REQUIRED if w in missing}
    with open(os.path.join(OUT, "zenodo_source.txt"), "w") as f:
        f.write(f"record_id\t{rid}\ndoi\t{doi}\ntitle\t{title}\n"
                f"version\t{version}\nfetched\t{', '.join(sorted(found))}\n")

    if not found:
        sys.exit("No wanted tables found in the deposit. Check the record id.")
    if missing:
        print(f"Not in the deposit: {', '.join(sorted(missing))}")
    if absent_required:
        sys.exit(f"Required input missing from the deposit: "
                 f"{', '.join(sorted(absent_required))}. The analysis resolves "
                 f"records to specimens using the crosswalk and cannot run "
                 f"without it.")
    print(f"{len(found)} table(s) in {OUT}/")


if __name__ == "__main__":
    main()
