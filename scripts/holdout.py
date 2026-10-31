"""Temporal hold-out test of the sampling-priority map.

Builds the coverage deficit surface from records collected up to a cutoff year,
then asks whether the sequence records collected afterwards landed in the
climate space the earlier data flagged as uncovered.

Two nulls, both drawing from the early occurrence cloud:

    occurrence   weighted by occurrence density — sampling the recorded range.
    sequencing   weighted by early sequence density — carrying on sequencing
                 wherever sequencing already happened. This guards the
                 alternative reading, that any new locality looks novel against
                 a sparse library whatever its climate.

Usage:
    python holdout.py [cutoff_year] [climate_table.csv]
    python holdout.py --sweep [climate_table.csv]

--sweep repeats the test across a range of cutoffs, so the cutoff is reported
as a robustness result rather than chosen.

Outputs to results/holdout/:
    holdout_summary.txt    the single-cutoff report
    holdout_records.csv    post-cutoff sequence records with uncovered fraction
    cutoff_sweep.csv       one row per cutoff (--sweep only)
"""
import csv
import os
import sys

import numpy as np

from deficit import (GRID, WEIGHTS, assign_roles, bilinear,
                     coverage_and_deficit, density_on_grid, load_records,
                     load_verdicts, newest_climate_table, pca2, weights_for)

OUTDIR = os.path.join("results", "holdout")
DRAWS = 2000
SEED = 20261031
SWEEP = (2006, 2008, 2010, 2012, 2014, 2015, 2016, 2018, 2020)
MIN_TRAIN_SEQ = 150   # below this the training library is too thin to be a fair test


def prepare(path):
    """Dated records with roles assigned, and their shared PCA climate space."""
    rows = load_records(path)
    assign_roles(rows)
    dated = [r for r in rows if (r.get("year") or "").isdigit()]
    if not dated:
        sys.exit("The climate table carries no usable year column. Re-run "
                 "climate_extract.py so that collection years are included.")
    for r in dated:
        r["_year"] = int(r["year"])
    scores, var, _ = pca2([r["_bio"] for r in dated])
    position = {id(r): scores[i] for i, r in enumerate(dated)}
    pad = 0.05 * (scores.max(axis=0) - scores.min(axis=0))
    lo, hi = scores.min(axis=0) - pad, scores.max(axis=0) + pad
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], GRID),
                         np.linspace(lo[1], hi[1], GRID))
    return rows, dated, position, var, gx, gy


def run_one(cutoff, dated, position, gx, gy, verdicts):
    """One hold-out test. Returns a result dict, or a dict with 'skipped'."""
    pc = lambda subset: np.array([position[id(r)] for r in subset])
    early = [r for r in dated if r["_year"] <= cutoff]
    early_occ = [r for r in early if r["_role"] == "occurrence"]
    early_seq = [r for r in early if r["_role"] == "sequence"]
    late_seq = [r for r in dated if r["_year"] > cutoff and r["_role"] == "sequence"]

    if len(early_occ) < 20 or len(early_seq) < 20 or len(late_seq) < 10:
        return {"cutoff": cutoff, "skipped": "too few records either side",
                "train_occurrence": len(early_occ), "train_sequence": len(early_seq),
                "test_sequence": len(late_seq)}

    occ_density, _ = density_on_grid(pc(early_occ), gx, gy)
    w, unmatched = weights_for(early_seq, verdicts)
    seq_density, area = density_on_grid(pc(early_seq), gx, gy, weights=w)
    coverage, _, train_index = coverage_and_deficit(occ_density, seq_density,
                                                    w.sum(), area)

    late_pc = pc(late_seq)
    observed = 1.0 - bilinear(coverage, gx, gy, late_pc[:, 0], late_pc[:, 1])
    obs_mean = float(observed.mean())

    occ_pc = pc(early_occ)
    occ_at = 1.0 - bilinear(coverage, gx, gy, occ_pc[:, 0], occ_pc[:, 1])
    rng = np.random.default_rng(SEED)
    nulls, pvals = {}, {}
    for label, surface in (("occurrence", occ_density), ("sequencing", seq_density)):
        dens_at = bilinear(surface, gx, gy, occ_pc[:, 0], occ_pc[:, 1])
        p = dens_at / dens_at.sum() if dens_at.sum() > 0 else None
        draws = np.array(
            [rng.choice(occ_at, size=len(late_seq), replace=True, p=p).mean()
             for _ in range(DRAWS)])
        nulls[label] = draws
        pvals[label] = max(float((draws >= obs_mean).mean()), 1.0 / DRAWS)

    return {"cutoff": cutoff, "skipped": "",
            "train_occurrence": len(early_occ), "train_sequence": len(early_seq),
            "test_sequence": len(late_seq), "unmatched_train": unmatched,
            "thin_training": len(early_seq) < MIN_TRAIN_SEQ,
            "train_deficit_index": train_index, "observed": obs_mean,
            "nulls": nulls, "pvals": pvals,
            "late_seq": late_seq, "observed_per_record": observed}


def report(res, var, path):
    nulls, pvals = res["nulls"], res["pvals"]
    if all(p < 0.05 for p in pvals.values()):
        verdict = "post-cutoff sequencing targeted the gaps, against both nulls"
    elif pvals["occurrence"] < 0.05:
        verdict = ("post-cutoff sequencing beat sampling the recorded range at "
                   "random, but not sequencing where sequencing already was: the "
                   "signal may be novelty of locality rather than of climate")
    else:
        verdict = ("post-cutoff sequencing did no better than sampling the "
                   "recorded range at random")
    lines = [
        f"climate table              {path}",
        f"cutoff year                {res['cutoff']}",
        f"training set               occurrence {res['train_occurrence']}, "
        f"sequence {res['train_sequence']}",
        f"test set                   {res['test_sequence']} sequence records "
        f"after {res['cutoff']}",
        f"verdicts matched (train)   "
        f"{res['train_sequence'] - res['unmatched_train']}/{res['train_sequence']}",
        f"PCA variance explained     PC1 {var[0]:.1%}, PC2 {var[1]:.1%}",
        f"training deficit index     {res['train_deficit_index']:.3f}",
        "",
        f"mean uncovered fraction at post-{res['cutoff']} sequence locations  "
        f"{res['observed']:.3f}",
    ] + [
        f"null '{label}'  {v.mean():.3f}"
        f"  [{np.percentile(v, 2.5):.3f}, {np.percentile(v, 97.5):.3f}]"
        f"  one-sided p {pvals[label]:.4f}"
        for label, v in nulls.items()
    ] + ["", f"reading: {verdict}"]
    if res["thin_training"]:
        lines += [f"caution: only {res['train_sequence']} training sequences, "
                  f"below the {MIN_TRAIN_SEQ} needed for a fair test"]
    return "\n".join(lines)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sweep = "--sweep" in sys.argv
    if sweep:
        path = args[0] if args else newest_climate_table()
    else:
        path = args[1] if len(args) > 1 else newest_climate_table()
    os.makedirs(OUTDIR, exist_ok=True)

    _, dated, position, var, gx, gy = prepare(path)
    verdicts = load_verdicts()

    if sweep:
        out = []
        for cutoff in SWEEP:
            res = run_one(cutoff, dated, position, gx, gy, verdicts)
            if res.get("skipped"):
                print(f"{cutoff}  skipped: {res['skipped']}")
                out.append({"cutoff": cutoff, "train_sequence": res["train_sequence"],
                            "test_sequence": res["test_sequence"],
                            "observed": "", "null_occurrence": "",
                            "p_occurrence": "", "null_sequencing": "",
                            "p_sequencing": "", "ratio_to_occurrence_null": "",
                            "note": res["skipped"]})
                continue
            n_occ = float(res["nulls"]["occurrence"].mean())
            n_seq = float(res["nulls"]["sequencing"].mean())
            out.append({
                "cutoff": cutoff, "train_sequence": res["train_sequence"],
                "test_sequence": res["test_sequence"],
                "observed": round(res["observed"], 4),
                "null_occurrence": round(n_occ, 4),
                "p_occurrence": round(res["pvals"]["occurrence"], 4),
                "null_sequencing": round(n_seq, 4),
                "p_sequencing": round(res["pvals"]["sequencing"], 4),
                "ratio_to_occurrence_null": round(res["observed"] / n_occ, 3)
                if n_occ else "",
                "note": "thin training library" if res["thin_training"] else ""})
            print(f"{cutoff}  train {res['train_sequence']:>5}  "
                  f"test {res['test_sequence']:>5}  "
                  f"observed {res['observed']:.3f}  "
                  f"null {n_occ:.3f}/{n_seq:.3f}  "
                  f"p {res['pvals']['occurrence']:.4f}/"
                  f"{res['pvals']['sequencing']:.4f}"
                  f"{'  (thin)' if res['thin_training'] else ''}")
        with open(os.path.join(OUTDIR, "cutoff_sweep.csv"), "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            wr.writeheader()
            wr.writerows(out)
        usable = [r for r in out if r["observed"] != "" and not r["note"]]
        if usable:
            ratios = [r["ratio_to_occurrence_null"] for r in usable]
            print(f"\n{len(usable)} cutoffs with an adequate training library; "
                  f"observed exceeds the occurrence null in "
                  f"{sum(1 for r in usable if r['observed'] > r['null_occurrence'])}"
                  f" of them")
            print(f"ratio to the occurrence null: {min(ratios):.2f} to "
                  f"{max(ratios):.2f}")
        print(f"\nSaved {os.path.join(OUTDIR, 'cutoff_sweep.csv')}")
        return

    cutoff = int(args[0]) if args else 2015
    res = run_one(cutoff, dated, position, gx, gy, verdicts)
    if res.get("skipped"):
        sys.exit(f"Cutoff {cutoff}: {res['skipped']} (training occurrence "
                 f"{res['train_occurrence']}, training sequence "
                 f"{res['train_sequence']}, test sequence {res['test_sequence']})")

    with open(os.path.join(OUTDIR, "holdout_records.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["source", "record_id", "taxon", "country", "year",
                     "lat", "lon", "uncovered_fraction"])
        for r, u in zip(res["late_seq"], res["observed_per_record"]):
            wr.writerow([r.get("source"), r.get("record_id"), r.get("taxon"),
                         r.get("country"), r.get("year"),
                         round(r["_lat"], 4), round(r["_lon"], 4),
                         round(float(u), 4)])

    text = report(res, var, path)
    with open(os.path.join(OUTDIR, "holdout_summary.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
