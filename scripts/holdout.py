"""Temporal hold-out test of the sampling-priority map.

Builds the coverage deficit surface from records collected up to a cutoff year,
then asks whether the sequence records collected afterwards actually landed in
the climate space the earlier data flagged as uncovered.

The null model draws the same number of locations from the occurrence records,
in proportion to occurrence density, so the test asks whether real sequencing
effort tracked the gaps better than sampling the recorded range at random.

Usage:
    python holdout.py [cutoff_year] [climate_table.csv]

Outputs to results/holdout/:
    holdout_summary.txt
    holdout_records.csv   post-cutoff sequence records with their uncovered fraction
"""
import csv
import os
import sys

import numpy as np

from deficit import (BIO, GRID, OCCURRENCE_SOURCES, SEQUENCE_SOURCES, WEIGHTS,
                     bilinear, coverage_and_deficit, density_on_grid,
                     load_records, load_verdicts, newest_climate_table, pca2,
                     weights_for)

OUTDIR = os.path.join("results", "holdout")
DRAWS = 2000
SEED = 20261031


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cutoff = int(args[0]) if args else 2015
    path = args[1] if len(args) > 1 else newest_climate_table()
    os.makedirs(OUTDIR, exist_ok=True)

    rows = load_records(path)
    dated = [r for r in rows if (r.get("year") or "").isdigit()]
    if not dated:
        sys.exit("The climate table carries no usable year column. Re-run "
                 "climate_extract.py so that collection years are included.")
    for r in dated:
        r["_year"] = int(r["year"])

    early = [r for r in dated if r["_year"] <= cutoff]
    late_seq = [r for r in dated
                if r["_year"] > cutoff and r.get("source") in SEQUENCE_SOURCES]
    early_occ = [r for r in early if r.get("source") in OCCURRENCE_SOURCES]
    early_seq = [r for r in early if r.get("source") in SEQUENCE_SOURCES]

    if len(early_occ) < 20 or len(early_seq) < 20 or len(late_seq) < 10:
        sys.exit(f"Not enough records either side of {cutoff}: "
                 f"early occurrence {len(early_occ)}, early sequence "
                 f"{len(early_seq)}, late sequence {len(late_seq)}")

    # One PCA over every dated record, so early and late sit in the same space.
    scores, var, _ = pca2([r["_bio"] for r in dated])
    index = {id(r): i for i, r in enumerate(dated)}
    pc = lambda subset: np.array([scores[index[id(r)]] for r in subset])

    pad = 0.05 * (scores.max(axis=0) - scores.min(axis=0))
    lo, hi = scores.min(axis=0) - pad, scores.max(axis=0) + pad
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], GRID),
                         np.linspace(lo[1], hi[1], GRID))

    occ_density, _ = density_on_grid(pc(early_occ), gx, gy)
    w, unmatched = weights_for(early_seq, load_verdicts())
    seq_density, area = density_on_grid(pc(early_seq), gx, gy, weights=w)
    coverage, _, index_value = coverage_and_deficit(occ_density, seq_density, w.sum(), area)

    late_pc = pc(late_seq)
    observed = 1.0 - bilinear(coverage, gx, gy, late_pc[:, 0], late_pc[:, 1])

    # Two nulls, both drawing from the early occurrence cloud.
    #   occurrence  weighted by occurrence density: sampling the recorded range.
    #   sequencing  weighted by early sequence density: carrying on sequencing
    #               wherever sequencing already happened. This one guards the
    #               alternative reading, that any new locality looks novel
    #               against a sparse library whatever the climate.
    occ_pc = pc(early_occ)
    occ_at = 1.0 - bilinear(coverage, gx, gy, occ_pc[:, 0], occ_pc[:, 1])
    rng = np.random.default_rng(SEED)
    nulls = {}
    for label, surface in (("occurrence", occ_density), ("sequencing", seq_density)):
        dens_at = bilinear(surface, gx, gy, occ_pc[:, 0], occ_pc[:, 1])
        p = dens_at / dens_at.sum() if dens_at.sum() > 0 else None
        nulls[label] = np.array(
            [rng.choice(occ_at, size=len(late_seq), replace=True, p=p).mean()
             for _ in range(DRAWS)])

    obs_mean = float(observed.mean())
    pvals = {k: max(float((v >= obs_mean).mean()), 1.0 / DRAWS)
             for k, v in nulls.items()}
    null = nulls["occurrence"]
    pval = pvals["occurrence"]

    with open(os.path.join(OUTDIR, "holdout_records.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["source", "record_id", "taxon", "country", "year",
                     "lat", "lon", "uncovered_fraction"])
        for r, u in zip(late_seq, observed):
            wr.writerow([r.get("source"), r.get("record_id"), r.get("taxon"),
                         r.get("country"), r.get("year"),
                         round(r["_lat"], 4), round(r["_lon"], 4), round(float(u), 4)])

    if all(p < 0.05 for p in pvals.values()):
        verdict = ("post-cutoff sequencing targeted the gaps, against both nulls")
    elif pvals["occurrence"] < 0.05:
        verdict = ("post-cutoff sequencing beat sampling the recorded range at "
                   "random, but not sequencing where sequencing already was: the "
                   "signal may be novelty of locality rather than of climate")
    else:
        verdict = ("post-cutoff sequencing did no better than sampling the "
                   "recorded range at random")
    lines = [
        f"climate table              {path}",
        f"cutoff year                {cutoff}",
        f"dated records              {len(dated)}/{len(rows)}",
        f"training set               occurrence {len(early_occ)}, sequence {len(early_seq)}",
        f"test set                   {len(late_seq)} sequence records after {cutoff}",
        f"verdicts matched (train)   {len(early_seq) - unmatched}/{len(early_seq)}",
        f"PCA variance explained     PC1 {var[0]:.1%}, PC2 {var[1]:.1%}",
        f"training deficit index     {index_value:.3f}",
        "",
        f"mean uncovered fraction at post-{cutoff} sequence locations  {obs_mean:.3f}",
    ] + [
        f"null '{label}'  {v.mean():.3f}"
        f"  [{np.percentile(v, 2.5):.3f}, {np.percentile(v, 97.5):.3f}]"
        f"  one-sided p {pvals[label]:.4f}"
        for label, v in nulls.items()
    ] + [
        "",
        f"reading: {verdict}",
    ]
    text = "\n".join(lines)
    with open(os.path.join(OUTDIR, "holdout_summary.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
