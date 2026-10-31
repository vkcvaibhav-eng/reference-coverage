"""Emit every number the manuscript quotes, computed from the analysis.

Table 1 of an earlier draft was typed by hand and did not match the data it
described. Nothing in the manuscript should be transcribed: this script writes
results/manuscript_numbers.json and a markdown fragment for Table 1, both
derived from the same code path as the analysis, so a figure in the text that
disagrees with the output is a build error rather than a reading error.

Usage:
    python manuscript_numbers.py [--genus NAME ...]
"""
import collections
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from canonical import Registry, collapse
from deficit import (BIO, MIN_RELIABLE, WEIGHTS, SHARED_CLUSTER, anchor_table,
                     assign_roles, label_accuracy, load_anchor_counts,
                     load_records, load_verdicts, newest_climate_table,
                     weights_for)

LABELS = [
    ("AGREES", "Agrees with anchor"),
    (SHARED_CLUSTER, "Agrees, shared cluster"),
    ("NO_ANCHOR", "No anchor"),
    ("UNNAMED_SUGGESTED", "Unnamed, one suggested"),
    ("UNNAMED_SEVERAL_POSSIBLE", "Unnamed, several possible"),
    ("DISAGREES", "Contradicted"),
    ("NO_VERDICT", "No verdict"),
]
MEANINGS = {
    "AGREES": "label concordant with an anchored specimen in the same cluster",
    SHARED_CLUSTER: "cluster holds k anchored species, so the label fits any one",
    "NO_ANCHOR": "no anchored specimen in the cluster; label unverifiable",
    "UNNAMED_SUGGESTED": "record unnamed, one species suggested by clustering",
    "UNNAMED_SEVERAL_POSSIBLE": "record unnamed, several species possible",
    "DISAGREES": "label contradicted by an anchored specimen",
    "NO_VERDICT": "not covered by the anchor check",
}


def summary_value(path, key, cast=float):
    """Pull a labelled value out of a summary text file."""
    if not os.path.exists(path):
        return None
    for line in open(path):
        if line.strip().startswith(key):
            tail = line.split(key, 1)[1].strip()
            for token in tail.replace(",", " ").split():
                try:
                    return cast(token)
                except ValueError:
                    continue
    return None


def genus_block(path, genus, registry):
    rows = load_records(path, genus)
    occ, seq, reclass = assign_roles(rows)
    seq, report = collapse(seq, registry)

    counts = collections.Counter(
        registry.verdict(r["record_id"]) or "NO_VERDICT" for r in seq)
    w, unmatched = weights_for(seq, load_verdicts(), registry=registry)

    weight_shown = {}
    for key, _ in LABELS:
        if key == SHARED_CLUSTER:
            ks = {registry.anchor_count(r["record_id"])
                  for r in seq if registry.verdict(r["record_id"]) == key}
            ks = {k for k in ks if k}
            weight_shown[key] = ("1/k" if len(ks) != 1
                                 else f"1/{ks.pop()} = {1 / max(ks | {2}):.2f}"
                                 if ks else "1/k")
        elif key == "NO_VERDICT":
            weight_shown[key] = "0.25"
        else:
            weight_shown[key] = f"{WEIGHTS[key]:.2f}"

    gbif = sum(1 for r in occ if r.get("source") == "GBIF")
    dd = os.path.join("results", "deficit", genus.lower())
    hd = os.path.join("results", "holdout", genus.lower())
    s = os.path.join(dd, "deficit_summary.txt")

    sweep = []
    sweep_path = os.path.join(hd, "cutoff_sweep.csv")
    if os.path.exists(sweep_path):
        sweep = [r for r in csv.DictReader(open(sweep_path)) if r["observed"]]
    adequate = [r for r in sweep if not r["note"]]

    priorities = []
    pri_path = os.path.join(dd, "priority_list.csv")
    if os.path.exists(pri_path):
        priorities = list(csv.DictReader(open(pri_path)))[:6]

    return {
        "genus": genus,
        "records_used": len(occ) + len(seq),
        "occurrence_rows": len(occ),
        "occurrence_from_gbif": gbif,
        "occurrence_reclassified": len(occ) - gbif,
        "reclassified_by_marker": reclass,
        "specimens": len(seq),
        "rows_before_collapse": len(seq) + report["collapsed"],
        "duplicates_collapsed": report["collapsed"],
        "redundant_share": round(report["collapsed"] /
                                 max(1, len(seq) + report["collapsed"]), 4),
        "coordinate_conflicts": report["coordinate_conflicts"],
        "verdict_counts": {k: int(counts.get(k, 0)) for k, _ in LABELS},
        "with_verdict": len(seq) - unmatched,
        "without_verdict": unmatched,
        "verdict_match_rate": round((len(seq) - unmatched) / max(1, len(seq)), 4),
        "weights_shown": weight_shown,
        "mean_verifiability": round(float(w.mean()), 4),
        "effective_specimens": round(float(w.sum()), 2),
        "deficit_counted": summary_value(s, "deficit index, naive"),
        "deficit_weighted": summary_value(s, "deficit index, weighted"),
        "pc1": summary_value(s, "PCA variance explained"),
        "training_deficit_2015": summary_value(
            os.path.join(hd, "holdout_summary.txt"), "training deficit index"),
        "sweep_adequate": len(adequate),
        "sweep_ratios": [float(r["ratio_to_occurrence_null"]) for r in adequate],
        "priorities": priorities,
    }


def table1_markdown(block):
    head = ("| Verdict | Weight | Meaning | *" + block["genus"] +
            "* specimens |\n|---|---|---|---|\n")
    body = "".join(
        f"| {label} | {block['weights_shown'][key]} | {MEANINGS[key]} | "
        f"{block['verdict_counts'][key]} |\n"
        for key, label in LABELS)
    total = sum(block["verdict_counts"].values())
    body += f"| **Total** | | | **{total}** |\n"
    return head + body


def main():
    genera = [a for a in sys.argv[1:] if not a.startswith("--")] or \
        ["Tetranychus", "Oligonychus"]
    path = newest_climate_table()
    registry = Registry(load_verdicts(), load_anchor_counts())

    at = anchor_table()
    checkable = [r for r in at
                 if r["verdict"] in ("AGREES", SHARED_CLUSTER, "DISAGREES")]
    contradicted = sum(1 for r in checkable if r["verdict"] == "DISAGREES")
    shared = sum(1 for r in checkable if r["verdict"] == SHARED_CLUSTER)

    out = {
        "climate_table": path,
        "crosswalk_aliases": registry.linked,
        "anchor_rows": len(at),
        "checkable": len(checkable),
        "contradicted": contradicted,
        "not_contradicted": len(checkable) - contradicted,
        "shared_cluster_among_not_contradicted": shared,
        "non_contradiction_rate": round(
            1 - contradicted / max(1, len(checkable)), 4),
        "coverage_scale": MIN_RELIABLE,
        "genera": {},
    }
    for g in genera:
        try:
            out["genera"][g] = genus_block(path, g, registry)
        except SystemExit:
            continue

    os.makedirs("results", exist_ok=True)
    with open(os.path.join("results", "manuscript_numbers.json"), "w") as f:
        json.dump(out, f, indent=2)
    for g, block in out["genera"].items():
        with open(os.path.join("results", f"table1_{g.lower()}.md"), "w") as f:
            f.write(table1_markdown(block))

    print(json.dumps({k: v for k, v in out.items() if k != "genera"}, indent=2))
    for g, b in out["genera"].items():
        print(f"\n--- {g} ---")
        for k in ("rows_before_collapse", "specimens", "duplicates_collapsed",
                  "redundant_share", "coordinate_conflicts", "with_verdict",
                  "without_verdict", "verdict_match_rate", "mean_verifiability",
                  "effective_specimens", "deficit_counted", "deficit_weighted",
                  "occurrence_rows", "occurrence_from_gbif",
                  "training_deficit_2015", "sweep_adequate"):
            print(f"  {k:<34} {b[k]}")
        print(f"  verdict_counts {b['verdict_counts']}")
        print(f"  sum of counts  {sum(b['verdict_counts'].values())}")


if __name__ == "__main__":
    main()
