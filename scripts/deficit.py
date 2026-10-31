"""Reliability-weighted coverage deficit of a DNA barcode reference library.

Compares where a taxon is recorded (occurrence records) against where the
reference library actually holds trustworthy sequences, in climate space, and
ranks the geographic cells where one new reference specimen would add most.

Usage:
    python deficit.py [climate_table.csv] [--sensitivity]

Defaults to the newest data/derived/climate_*.csv. Identification verdicts are
read from data/external/record_anchor_check.csv when present; records with no
verdict take DEFAULT_WEIGHT and the count is reported.

Outputs to results/deficit/:
    deficit_summary.txt      headline numbers
    climate_space.csv        PCA grid with occurrence density and both deficits
    geographic_deficit.csv   2-degree cells with mean uncovered fraction
    priority_list.csv        ranked cells where a new specimen adds most
    weight_sensitivity.csv   deficit index under alternative weight schemes
"""
import csv
import glob
import os
import sys

import numpy as np
from scipy.stats import gaussian_kde

BIO = [f"bio{i}" for i in range(1, 20)]
SEQUENCE_SOURCES = {"BOLD", "GenBank"}
OCCURRENCE_SOURCES = {"GBIF"}

# Only the barcode region of COI can supply reference coverage for a COI
# library. A georeferenced specimen record whose sequence is COI-3P, another
# gene, or absent is an occurrence record, not reference coverage, and is
# counted on the occurrence side. Records from a table with no marker column
# are left on the sequence side, since the preceding study's query restricted
# them to COI.
BARCODE_MARKERS = {"COI-5P", "COI", ""}

WEIGHTS = {
    "AGREES": 1.00,
    "AGREES_SHARED_CLUSTER": 0.50,
    "NO_ANCHOR": 0.25,
    "UNNAMED_SUGGESTED": 0.25,
    "UNNAMED_SEVERAL_POSSIBLE": 0.10,
    "DISAGREES": 0.00,
}
DEFAULT_WEIGHT = 0.25

# Alternative schemes for --sensitivity, as (weight map, weight for records with
# no verdict). "naive" is the published convention: every record counts fully,
# whatever its label is worth. "unverified_trusted" is the most conservative
# test of the method's own claim, giving unmatched records full credit.
SCHEMES = {
    "naive": ({k: 1.0 for k in WEIGHTS}, 1.00),
    "default": (WEIGHTS, DEFAULT_WEIGHT),
    "unverified_trusted": (WEIGHTS, 1.00),
    "strict": ({"AGREES": 1.0, "AGREES_SHARED_CLUSTER": 0.25, "NO_ANCHOR": 0.0,
                "UNNAMED_SUGGESTED": 0.0, "UNNAMED_SEVERAL_POSSIBLE": 0.0,
                "DISAGREES": 0.0}, 0.0),
    "lenient": ({"AGREES": 1.0, "AGREES_SHARED_CLUSTER": 0.75, "NO_ANCHOR": 0.5,
                 "UNNAMED_SUGGESTED": 0.5, "UNNAMED_SEVERAL_POSSIBLE": 0.25,
                 "DISAGREES": 0.0}, 0.5),
    "anchored_only": ({"AGREES": 1.0, "AGREES_SHARED_CLUSTER": 0.0, "NO_ANCHOR": 0.0,
                       "UNNAMED_SUGGESTED": 0.0, "UNNAMED_SEVERAL_POSSIBLE": 0.0,
                       "DISAGREES": 0.0}, 0.0),
}

GRID = 120          # grid points per PCA axis
CELL = 2.0          # geographic cell size in degrees
N_PRIORITY = 20     # length of the priority list

# How many trustworthy reference sequences from one neighbourhood of climate
# space before that climate counts as represented in the library. There is no
# standard value; it is reported as a parameter and swept by --sensitivity.
MIN_RELIABLE = 5.0
THRESHOLD_SWEEP = (1.0, 3.0, 5.0, 10.0, 20.0)

OUTDIR = os.path.join("results", "deficit")


def newest_climate_table():
    paths = sorted(glob.glob(os.path.join("data", "derived", "climate_*.csv")))
    if not paths:
        paths = sorted(glob.glob(os.path.join("data", "climate_*.csv")))
    if not paths:
        sys.exit("No climate table found. Run climate_extract.py first.")
    return paths[-1]


def genus_of(taxon):
    """Leading word of a taxon string; empty for BOLD BIN labels and blanks."""
    t = (taxon or "").strip()
    if not t or t.startswith("BOLD:"):
        return ""
    return t.split()[0]


def load_records(path, genus=None):
    """Climate table rows that carry a full BIO1-BIO19 vector, deduplicated.

    genus restricts both sides to one genus. A sequence of one genus supplies no
    reference coverage for another, so pooling genera understates the gap; the
    filter is the correct default for a per-genus result.
    """
    seen, rows = set(), []
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (r.get("source"), r.get("record_id"))
            if key in seen:
                continue
            if genus and genus_of(r.get("taxon")).lower() != genus.lower():
                continue
            if any(r.get(b) in (None, "") for b in BIO):
                continue
            try:
                r["_bio"] = [float(r[b]) for b in BIO]
                r["_lat"] = float(r["lat"])
                r["_lon"] = float(r["lon"])
            except (TypeError, ValueError):
                continue
            seen.add(key)
            rows.append(r)
    return rows


def assign_roles(rows):
    """Split records into occurrence and reference-sequence roles.

    Sets r["_role"] on every row. A sequence-source record whose marker is not
    the barcode region is moved to the occurrence side: it locates the taxon but
    adds nothing to a COI reference library. Returns the two lists and a count
    of what was moved, by marker.
    """
    reclassified = {}
    for r in rows:
        if (r.get("source") in SEQUENCE_SOURCES
                and "marker" in r
                and r.get("marker", "") not in BARCODE_MARKERS):
            label = r.get("marker") or "(no sequence)"
            reclassified[label] = reclassified.get(label, 0) + 1
            r["_role"] = "occurrence"
        elif r.get("source") in SEQUENCE_SOURCES:
            r["_role"] = "sequence"
        else:
            r["_role"] = "occurrence"
    return ([r for r in rows if r["_role"] == "occurrence"],
            [r for r in rows if r["_role"] == "sequence"],
            reclassified)


def load_verdicts():
    for path in (os.path.join("data", "external", "record_anchor_check.csv"),
                 os.path.join("data", "derived", "record_anchor_check.csv")):
        if os.path.exists(path):
            with open(path) as f:
                return {r["record_id"]: r["verdict"] for r in csv.DictReader(f)}
    return {}


def pca2(matrix):
    """Two leading principal components of standardised columns."""
    x = np.asarray(matrix, dtype=float)
    mean = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd == 0] = 1.0
    z = (x - mean) / sd
    _, s, vt = np.linalg.svd(z, full_matrices=False)
    scores = z @ vt[:2].T
    var = (s ** 2) / (s ** 2).sum()
    return scores, var[:2], (mean, sd, vt[:2])


def density_on_grid(points, gx, gy, weights=None):
    """Gaussian KDE on a grid, with the effective area of its kernel.

    Returns (density, kernel_area). The density integrates to 1; kernel_area is
    the climate-space area of one bandwidth, which converts a density into an
    expected record count per neighbourhood and is independent of grid spacing.
    """
    shape = gx.shape
    if len(points) < 3:
        return np.zeros(shape), 1.0
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if weights.sum() <= 0:
            return np.zeros(shape), 1.0
        weights = weights / weights.sum()
    try:
        kde = gaussian_kde(np.asarray(points).T, weights=weights)
    except np.linalg.LinAlgError:
        return np.zeros(shape), 1.0
    flat = np.vstack([gx.ravel(), gy.ravel()])
    area = 2.0 * np.pi * float(np.sqrt(np.linalg.det(kde.covariance)))
    return kde(flat).reshape(shape), area


def coverage_and_deficit(occ_density, seq_density, effective_n, kernel_area,
                         min_reliable=MIN_RELIABLE):
    """Absolute coverage and the deficit surface it implies.

    lambda(z) is the effective number of trustworthy reference sequences within
    one bandwidth of a point in climate space, and coverage saturates as that
    number approaches min_reliable:

        lambda(z)   = effective_n * density(z) * kernel_area
        coverage(z) = 1 - exp(-lambda(z) / min_reliable)

    effective_n is the sum of reliability weights, so reducing confidence in the
    labels lowers coverage in absolute terms rather than merely redistributing
    it. The index is the share of recorded occurrence sitting in climate space
    without trustworthy reference sequences: 0 is full coverage, 1 is none.
    """
    lam = float(effective_n) * seq_density * float(kernel_area)
    coverage = 1.0 - np.exp(-lam / float(min_reliable))
    total = occ_density.sum()
    occ_mass = occ_density / total if total > 0 else occ_density
    deficit = occ_mass * (1.0 - coverage)
    index = float(deficit.sum()) if total > 0 else float("nan")
    return coverage, deficit, index


def bilinear(surface, gx, gy, x, y):
    """Sample a grid surface at arbitrary points by nearest-node lookup."""
    xs, ys = gx[0, :], gy[:, 0]
    ix = np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 1)
    iy = np.clip(np.searchsorted(ys, y) - 1, 0, len(ys) - 1)
    return surface[iy, ix]


def weights_for(seq_rows, verdicts, scheme=None, unmatched_weight=DEFAULT_WEIGHT):
    """Reliability weight per sequence record, and how many had no verdict."""
    scheme = WEIGHTS if scheme is None else scheme
    out, unmatched = [], 0
    for r in seq_rows:
        verdict = verdicts.get(r.get("record_id", ""))
        if verdict is None:
            unmatched += 1
            out.append(unmatched_weight)
        else:
            out.append(scheme.get(verdict, unmatched_weight))
    return np.array(out, dtype=float), unmatched


def priority_list(occ_rows, occ_scores, uncovered, pc_bandwidth):
    """Greedy cells, suppressing cells climatically close to one already chosen."""
    cells = {}
    for row, pc, u in zip(occ_rows, occ_scores, uncovered):
        key = (int(np.floor(row["_lat"] / CELL) * CELL),
               int(np.floor(row["_lon"] / CELL) * CELL))
        cells.setdefault(key, {"n": 0, "u": [], "pc": [], "countries": {}})
        cells[key]["n"] += 1
        cells[key]["u"].append(float(u))
        cells[key]["pc"].append(pc)
        country = (row.get("country") or "").strip()
        if country:
            cells[key]["countries"][country] = cells[key]["countries"].get(country, 0) + 1

    live = {}
    for key, c in cells.items():
        live[key] = {"n": c["n"], "u": np.array(c["u"]),
                     "pc": np.array(c["pc"]),
                     "centroid": np.array(c["pc"]).mean(axis=0),
                     "country": max(c["countries"], key=c["countries"].get)
                                if c["countries"] else ""}

    chosen = []
    for _ in range(min(N_PRIORITY, len(live))):
        best, best_score = None, -1.0
        for key, c in live.items():
            score = float(c["u"].sum())
            if score > best_score:
                best, best_score = key, score
        if best is None or best_score <= 0:
            break
        c = live[best]
        chosen.append({"lat_bin": best[0], "lon_bin": best[1],
                       "country": c["country"], "occurrence_records": c["n"],
                       "uncovered_mass": round(best_score, 4),
                       "mean_uncovered_fraction": round(float(c["u"].mean()), 4)})
        # A specimen from the chosen cell also covers climatically similar cells.
        centroid = c["centroid"]
        for key, other in live.items():
            d = np.linalg.norm(other["pc"] - centroid, axis=1)
            other["u"] = other["u"] * (1.0 - np.exp(-0.5 * (d / pc_bandwidth) ** 2))
        del live[best]
    return chosen


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sensitivity = "--sensitivity" in sys.argv
    genus = None
    if "--genus" in sys.argv:
        genus = sys.argv[sys.argv.index("--genus") + 1]
        args = [a for a in args if a.lower() != genus.lower()]
    path = args[0] if args else newest_climate_table()
    outdir = os.path.join(OUTDIR, genus.lower()) if genus else OUTDIR
    os.makedirs(outdir, exist_ok=True)

    rows = load_records(path, genus)
    occ_rows, seq_rows, reclassified = assign_roles(rows)
    if not occ_rows or not seq_rows:
        sys.exit(f"Need both occurrence and sequence records; got "
                 f"{len(occ_rows)} and {len(seq_rows)} from {path}")

    verdicts = load_verdicts()
    scores, var, _ = pca2([r["_bio"] for r in rows])
    occ_mask = np.array([r["_role"] == "occurrence" for r in rows])
    seq_mask = np.array([r["_role"] == "sequence" for r in rows])
    occ_scores, seq_scores = scores[occ_mask], scores[seq_mask]

    pad = 0.05 * (scores.max(axis=0) - scores.min(axis=0))
    lo, hi = scores.min(axis=0) - pad, scores.max(axis=0) + pad
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], GRID),
                         np.linspace(lo[1], hi[1], GRID))

    occ_density, _ = density_on_grid(occ_scores, gx, gy)
    naive_density, naive_area = density_on_grid(seq_scores, gx, gy)
    w, unmatched = weights_for(seq_rows, verdicts)
    weighted_density, wtd_area = density_on_grid(seq_scores, gx, gy, weights=w)

    naive_cov, naive_def, naive_index = coverage_and_deficit(
        occ_density, naive_density, len(seq_scores), naive_area)
    wtd_cov, wtd_def, wtd_index = coverage_and_deficit(
        occ_density, weighted_density, w.sum(), wtd_area)

    # Bounds on the weighted index from the records that carry no verdict:
    # giving them full credit, and giving them none.
    bounds = {}
    for label, unmatched_w in (("unverified trusted", 1.0),
                               ("unverified discarded", 0.0)):
        bw, _ = weights_for(seq_rows, verdicts, WEIGHTS, unmatched_w)
        bd, ba = density_on_grid(seq_scores, gx, gy, weights=bw)
        bounds[label] = coverage_and_deficit(occ_density, bd, bw.sum(), ba)[2]

    # Per-occurrence-record uncovered fraction under the weighted surface.
    uncovered = 1.0 - bilinear(wtd_cov, gx, gy, occ_scores[:, 0], occ_scores[:, 1])
    bandwidth = float(np.mean([gx[0, 1] - gx[0, 0], gy[1, 0] - gy[0, 0]]) * GRID / 12)

    with open(os.path.join(outdir, "climate_space.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["pc1", "pc2", "occurrence_density", "coverage_naive",
                     "coverage_weighted", "deficit_naive", "deficit_weighted"])
        for i in range(GRID):
            for j in range(GRID):
                wr.writerow([round(gx[i, j], 4), round(gy[i, j], 4),
                             f"{occ_density[i, j]:.6g}", f"{naive_cov[i, j]:.4f}",
                             f"{wtd_cov[i, j]:.4f}", f"{naive_def[i, j]:.6g}",
                             f"{wtd_def[i, j]:.6g}"])

    geo = {}
    for row, u in zip(occ_rows, uncovered):
        key = (int(np.floor(row["_lat"] / CELL) * CELL),
               int(np.floor(row["_lon"] / CELL) * CELL))
        geo.setdefault(key, []).append(float(u))
    with open(os.path.join(outdir, "geographic_deficit.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["lat_bin", "lon_bin", "occurrence_records",
                     "mean_uncovered_fraction"])
        for (la, lo_), vals in sorted(geo.items()):
            wr.writerow([la, lo_, len(vals), round(float(np.mean(vals)), 4)])

    picks = priority_list(occ_rows, occ_scores, uncovered, bandwidth)
    with open(os.path.join(outdir, "priority_list.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["rank", "lat_bin", "lon_bin", "country",
                                           "occurrence_records", "uncovered_mass",
                                           "mean_uncovered_fraction"])
        wr.writeheader()
        for i, p in enumerate(picks, 1):
            wr.writerow({"rank": i, **p})

    sens_rows, thresh_rows = [], []
    if sensitivity:
        for name, (scheme, unmatched_w) in SCHEMES.items():
            sw, _ = weights_for(seq_rows, verdicts, scheme, unmatched_w)
            d, area = density_on_grid(seq_scores, gx, gy, weights=sw)
            _, _, idx = coverage_and_deficit(occ_density, d, sw.sum(), area)
            sens_rows.append({"scheme": name,
                              "mean_weight": round(float(sw.mean()), 4),
                              "effective_records": round(float(sw.sum()), 1),
                              "deficit_index": round(idx, 4)})
        with open(os.path.join(outdir, "weight_sensitivity.csv"), "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=["scheme", "mean_weight",
                                               "effective_records", "deficit_index"])
            wr.writeheader()
            wr.writerows(sens_rows)

        for k in THRESHOLD_SWEEP:
            _, _, n_idx = coverage_and_deficit(occ_density, naive_density,
                                               len(seq_scores), naive_area, k)
            _, _, w_idx = coverage_and_deficit(occ_density, weighted_density,
                                               w.sum(), wtd_area, k)
            thresh_rows.append({"min_reliable": k,
                                "deficit_naive": round(n_idx, 4),
                                "deficit_weighted": round(w_idx, 4),
                                "relative_inflation":
                                    round(w_idx / n_idx - 1, 4) if n_idx else ""})
        with open(os.path.join(outdir, "threshold_sensitivity.csv"), "w",
                  newline="") as f:
            wr = csv.DictWriter(f, fieldnames=["min_reliable", "deficit_naive",
                                               "deficit_weighted",
                                               "relative_inflation"])
            wr.writeheader()
            wr.writerows(thresh_rows)

    reclass_note = (", ".join(f"{v} {k}" for k, v in
                              sorted(reclassified.items(), key=lambda kv: -kv[1]))
                    if reclassified else "none")
    lines = [
        f"climate table            {path}",
        f"genus filter             {genus or 'none (genera pooled)'}",
        f"records used             {len(rows)} (occurrence {len(occ_rows)}, sequence {len(seq_rows)})",
        f"moved to occurrence      {reclass_note}",
        f"verdicts matched         {len(seq_rows) - unmatched}/{len(seq_rows)}"
        f" ({unmatched} took the default weight of {DEFAULT_WEIGHT})",
        f"mean reliability weight  {w.mean():.3f}",
        f"effective sequences      {w.sum():.1f} of {len(seq_rows)} records",
        f"coverage threshold       {MIN_RELIABLE:g} reliable sequences per neighbourhood",
        f"PCA variance explained   PC1 {var[0]:.1%}, PC2 {var[1]:.1%}",
        "",
        f"deficit index, naive     {naive_index:.3f}",
        f"deficit index, weighted  {wtd_index:.3f}",
        f"inflation of the gap     {wtd_index - naive_index:+.3f}"
        + (f" ({(wtd_index / naive_index - 1) * 100:+.1f}% relative)"
           if naive_index else ""),
        f"bounds on the weighted index from unverified records:",
        f"  {bounds['unverified trusted']:.3f} if they are trusted in full",
        f"  {bounds['unverified discarded']:.3f} if they are discarded",
        "",
        f"occurrence records in climate space with <10% effective coverage: "
        f"{int((uncovered > 0.9).sum())}/{len(uncovered)}"
        f" ({(uncovered > 0.9).mean():.1%})",
        f"geographic cells ({CELL:g} degrees) holding occurrences: {len(geo)}",
    ]
    if sens_rows:
        lines += ["", "weight sensitivity:"]
        lines += [f"  {r['scheme']:<20} effective sequences {r['effective_records']:>7.1f}"
                  f"  deficit index {r['deficit_index']:.3f}" for r in sens_rows]
    if thresh_rows:
        lines += ["", "coverage-threshold sensitivity (naive -> weighted):"]
        lines += [f"  {r['min_reliable']:>5g}  {r['deficit_naive']:.3f} -> "
                  f"{r['deficit_weighted']:.3f}"
                  f"  ({r['relative_inflation']:+.1%} inflation)"
                  if r["relative_inflation"] != "" else
                  f"  {r['min_reliable']:>5g}  {r['deficit_naive']:.3f} -> "
                  f"{r['deficit_weighted']:.3f}" for r in thresh_rows]
    if picks:
        lines += ["", f"top {min(5, len(picks))} sampling priorities:"]
        lines += [f"  {i}. {p['country'] or 'unknown'} "
                  f"({p['lat_bin']}, {p['lon_bin']})  "
                  f"uncovered mass {p['uncovered_mass']}"
                  for i, p in enumerate(picks[:5], 1)]

    text = "\n".join(l for l in lines if l is not None)
    with open(os.path.join(outdir, "deficit_summary.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
