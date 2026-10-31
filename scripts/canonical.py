"""Resolve public sequence records to specimens.

GenBank and BOLD hold the same physical specimen under different identifiers,
and a record table that keys on (source, record_id) counts such a specimen
twice: once as an accession and once as a process ID. It also fails to find a
verdict filed under the other identifier.

This module reads the crosswalk published by the preceding study, which carries
`genbank_accession` and `bold_processid` side by side, and builds

    alias -> canonical specimen key

so that duplicates collapse and a verdict filed under any identifier is found.
Where two representations of one specimen disagree on coordinates, the conflict
is reported rather than silently resolved.
"""
import csv
import glob
import os
import re

CROSSWALK_GLOBS = (
    os.path.join("data", "external", "*merged_records.csv"),
    os.path.join("data", "derived", "*merged_records.csv"),
)
ID_FIELDS = ("record_id", "genbank_accession", "bold_processid")
VERSION = re.compile(r"\.\d+$")


def normalise(identifier):
    """Strip a GenBank version suffix so KP123456.1 and KP123456 match."""
    t = (identifier or "").strip()
    return VERSION.sub("", t) if t else ""


def crosswalk_rows():
    """Every crosswalk found, across all genera. Missing ones are an error."""
    rows, files = [], []
    for pattern in CROSSWALK_GLOBS:
        for path in sorted(glob.glob(pattern)):
            files.append(path)
            with open(path) as f:
                rows += list(csv.DictReader(f))
    if not rows:
        raise FileNotFoundError(
            "No accession/process-ID crosswalk found under data/external/. "
            "Records cannot be resolved to specimens without it; run "
            "fetch_zenodo.py first.")
    crosswalk_rows.files = files
    return rows


def build_aliases():
    """alias -> canonical key, using union-find over the crosswalk rows."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            # Prefer an accession-looking key as the representative.
            if len(ra) <= len(rb):
                parent[rb] = ra
            else:
                parent[ra] = rb

    for r in crosswalk_rows():
        ids = [normalise(r.get(f)) for f in ID_FIELDS]
        ids = [i for i in ids if i]
        for other in ids[1:]:
            union(ids[0], other)

    return {k: find(k) for k in parent}


class Registry:
    """Canonical keys and verdict lookup across every known identifier."""

    def __init__(self, verdicts=None, anchor_counts=None):
        self.alias = build_aliases()
        self.verdicts = {normalise(k): v for k, v in (verdicts or {}).items()}
        self.anchor_counts = {normalise(k): v
                              for k, v in (anchor_counts or {}).items()}
        self._by_canonical = {}
        for ident, canon in self.alias.items():
            if ident in self.verdicts:
                self._by_canonical.setdefault(canon, self.verdicts[ident])
            if ident in self.anchor_counts:
                self._by_canonical.setdefault(canon + "#k",
                                              self.anchor_counts[ident])

    @property
    def linked(self):
        return len(self.alias)

    def key(self, record_id):
        n = normalise(record_id)
        return self.alias.get(n, n)

    def verdict(self, record_id):
        n = normalise(record_id)
        if n in self.verdicts:
            return self.verdicts[n]
        return self._by_canonical.get(self.alias.get(n, n))

    def anchor_count(self, record_id):
        n = normalise(record_id)
        if n in self.anchor_counts:
            return self.anchor_counts[n]
        return self._by_canonical.get(self.alias.get(n, n) + "#k")


def collapse(rows, registry, source_rank=("GenBank", "BOLD", "GBIF")):
    """One row per specimen, with a report of what was collapsed.

    Among rows sharing a canonical key, the kept row is the one with a verdict
    if any has one, then by source preference, then by record_id, so the result
    does not depend on input order. Coordinate disagreements between collapsed
    rows are counted and returned rather than hidden.
    """
    rank = {s: i for i, s in enumerate(source_rank)}
    groups = {}
    for r in rows:
        groups.setdefault(registry.key(r.get("record_id", "")), []).append(r)

    kept, collapsed, coord_conflicts, conflict_rows = [], 0, 0, []
    for canon, members in groups.items():
        if len(members) > 1:
            collapsed += len(members) - 1
            coords = {(round(float(m["lat"]), 4), round(float(m["lon"]), 4))
                      for m in members if m.get("lat") not in (None, "")}
            if len(coords) > 1:
                coord_conflicts += 1
                conflict_rows.append({
                    "canonical": canon,
                    "records": ";".join(f"{m.get('source')}:{m.get('record_id')}"
                                        for m in members),
                    "coordinates": ";".join(f"{a},{b}" for a, b in sorted(coords)),
                })
        members.sort(key=lambda m: (
            0 if registry.verdict(m.get("record_id", "")) else 1,
            rank.get(m.get("source"), 99),
            m.get("record_id", ""),
        ))
        best = members[0]
        best["_canonical"] = canon
        best["_merged_from"] = len(members)
        kept.append(best)

    return kept, {"collapsed": collapsed, "coordinate_conflicts": coord_conflicts,
                  "conflict_rows": conflict_rows, "specimens": len(kept)}
