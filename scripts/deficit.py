"""Verifiable coverage of a DNA barcode reference library, in climate space.

Compares where a taxon is recorded (occurrence records) against where the
library holds reference sequences whose species label can be *verified* against
a morphologically identified specimen, and ranks the geographic cells where one
new reference specimen would add most.

The quantity measured is verifiable coverage, not accuracy. A record whose label
cannot be checked may well be correct — across this anchor table unchecked
labels are right about nine times in ten — but a user of the library has no way
to tell which, so unverifiable records cannot carry the weight of a reference.
The claim is that usable coverage is far smaller than the record count implies,
and it needs no assumption that the labels are wrong.

Usage:
    python deficit.py [climate_table.csv] [--sensitivity] [--genus NAME]

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
from canonical import Registry, collapse

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

# A weight of None is derived from the data rather than chosen: a record whose
# cluster holds k anchored species can be assigned to any one of them, so its
# label is worth 1/k. Every other value is a stated position on what an
# unverifiable record is worth as a reference, not an estimate of how often such
# labels are right.
SHARED_CLUSTER = "AGREES_SHARED_CLUSTER"
WEIGHTS = {
    "AGREES": 1.00,
    SHARED_CLUSTER: None,      # derived as 1/k
    "NO_ANCHOR": 0.25,
    "UNNAMED_SUGGESTED": 0.25,
    "UNNAMED_SEVERAL_POSSIBLE": 0.10,
    "DISAGREES": 0.00,
}
DEFAULT_WEIGHT = 0.25
FALLBACK_SHARED = 0.50     # used when the anchored-species count is unavailable

# Alternative schemes for --sensitivity, as (weight map, weight for records with
# no verdict). "naive" is the published convention: every record counts fully,
# whatever its label is worth. "unverified_trusted" is the most conservative
# test of the method's own claim, giving unmatched records full credit.
SCHEMES = {
    "naive": ({k: 1.0 for k in WEIGHTS}, 1.00),
    "default": (WEIGHTS, DEFAULT_WEIGHT),
    "unverified_trusted": (WEIGHTS, 1.00),
    # The floor: every weight set to the chance the label is simply correct,
    # unverifiable records included at the observed base rate. If the claim
    # survives here it survives the strongest objection to it.
    "non_contradiction": (None, None),
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


def anchor_table():
    """Every per-record anchor check in the deposit, across all genera.

    The preceding study files one table per genus. Reading a single file by
    basename silently drops the others, which is how an entire genus can appear
    to carry no identification evidence when its table exists.
    """
    rows, seen = [], set()
    for pattern in (os.path.join("data", "external", "*record_anchor_check.csv"),
                    os.path.join("data", "derived", "*record_anchor_check.csv")):
        for path in sorted(glob.glob(pattern)):
            if os.path.basename(path) in seen:
                continue
            seen.add(os.path.basename(path))
            with open(path) as f:
                rows += list(csv.DictReader(f))
    return rows


def load_anchor_counts():
    """record_id -> number of anchored species sharing that record's cluster."""
    counts = {}
    for r in anchor_table():
        names = [x.strip() for x in (r.get("anchored_species") or "").split(";")]
        names = [x for x in names if x]
        if names:
            counts[r["record_id"]] = len(names)
    return counts


def label_accuracy(verdicts):
    """Observed share of checkable labels that an anchor did not contradict.

    Estimated only from records whose cluster contains an anchored specimen,
    which are the better-studied part of the library, so applying it to
    unanchored records is an optimistic extrapolation rather than a neutral one.
    Reported as the floor of the claim, never as the working weight.
    """
    checkable = [v for v in verdicts.values()
                 if v in ("AGREES", SHARED_CLUSTER, "DISAGREES")]
    if not checkable:
        return None, 0, 0
    wrong = sum(1 for v in checkable if v == "DISAGREES")
    return 1.0 - wrong / len(checkable), len(checkable), wrong


def load_verdicts():
    return {r["record_id"]: r["verdict"] for r in anchor_table() if r.get("verdict")}


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


def reference_covariance(points):
    """Smoothing covariance, fixed once and reused by every weight scheme.

    SciPy derives its bandwidth from the weighted effective sample size, so a
    weighted and an unweighted density differ in how far credit is spread as
    well as in how much credit there is. Comparing them would then confound the
    two. The bandwidth is therefore computed once from the unweighted points by
    Scott's rule and held fixed across schemes.
    """
    x = np.asarray(points, dtype=float)
    if len(x) < 3:
        return np.eye(2)
    cov = np.cov(x.T)
    factor = len(x) ** (-1.0 / (x.shape[1] + 4))     # Scott
    return np.atleast_2d(cov) * factor ** 2


def density_on_grid(points, gx, gy, weights=None, covariance=None):
    """Gaussian density on a grid under an explicit covariance.

    Returns (density, kernel_area). The density integrates to 1; kernel_area is
    the climate-space area of one bandwidth, which converts a density into an
    expected specimen count per neighbourhood, independent of grid spacing.
    """
    shape = gx.shape
    x = np.asarray(points, dtype=float)
    if len(x) < 3:
        return np.zeros(shape), 1.0
    H = reference_covariance(x) if covariance is None else np.asarray(covariance)
    try:
        inv = np.linalg.inv(H)
        det = np.linalg.det(H)
    except np.linalg.LinAlgError:
        return np.zeros(shape), 1.0
    if det <= 0:
        return np.zeros(shape), 1.0

    if weights is None:
        w = np.ones(len(x))
    else:
        w = np.asarray(weights, dtype=float)
        if w.sum() <= 0:
            return np.zeros(shape), 1.0
    w = w / w.sum()

    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    out = np.zeros(len(nodes))
    norm = 1.0 / (2.0 * np.pi * np.sqrt(det))
    for start in range(0, len(x), 512):                 # chunked to bound memory
        chunk, wc = x[start:start + 512], w[start:start + 512]
        d = nodes[:, None, :] - chunk[None, :, :]
        m = np.einsum("nki,ij,nkj->nk", d, inv, d)
        out += (np.exp(-0.5 * m) * wc).sum(axis=1)
    area = 2.0 * np.pi * float(np.sqrt(det))
    return (out * norm).reshape(shape), area


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


def sample_surface(surface, gx, gy, x, y):
    """Sample a grid surface at arbitrary points by nearest-node lookup."""
    xs, ys = gx[0, :], gy[:, 0]
    ix = np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 1)
    iy = np.clip(np.searchsorted(ys, y) - 1, 0, len(ys) - 1)
    return surface[iy, ix]


def weights_for(seq_rows, verdicts, scheme=None, unmatched_weight=DEFAULT_WEIGHT,
                anchor_counts=None, registry=None):
    """Verifiability weight per sequence record, and how many had no verdict.

    A scheme value of None is derived rather than chosen: for a record whose
    cluster holds k anchored species the weight is 1/k.
    """
    scheme = WEIGHTS if scheme is None else scheme
    anchor_counts = anchor_counts or {}
    out, unmatched = [], 0
    for r in seq_rows:
        rid = r.get("record_id", "")
        verdict = registry.verdict(rid) if registry else verdicts.get(rid)
        if verdict is None:
            unmatched += 1
            out.append(unmatched_weight)
            continue
        value = scheme.get(verdict, unmatched_weight)
        if value is None:
            k = (registry.anchor_count(rid) if registry
                 else anchor_counts.get(rid))
            value = (1.0 / k) if k else FALLBACK_SHARED
        out.append(value)
    return np.array(out, dtype=float), unmatched


def marginal_gain(occ_density, lam, gx, gy, centroid, covariance, min_reliable):
    """Exact fall in the deficit index from adding one verified specimen.

    Adding a specimen of unit weight at `centroid` raises the expected count by
    one kernel, so lambda_new = lambda + A * N(z; centroid, H), where A is the
    kernel area. No approximation and no suppression constant is involved.
    """
    inv = np.linalg.inv(covariance)
    det = np.linalg.det(covariance)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    d = nodes - np.asarray(centroid)[None, :]
    kernel = np.exp(-0.5 * np.einsum("ni,ij,nj->n", d, inv, d)) / (2 * np.pi * np.sqrt(det))
    lam_new = lam + (2 * np.pi * np.sqrt(det)) * kernel.reshape(gx.shape)
    mass = occ_density / occ_density.sum() if occ_density.sum() > 0 else occ_density
    before = float((mass * np.exp(-lam / min_reliable)).sum())
    after = float((mass * np.exp(-lam_new / min_reliable)).sum())
    return before - after


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
    registry = Registry(load_verdicts(), load_anchor_counts())
    seq_rows, collapse_report = collapse(seq_rows, registry)
    rows = occ_rows + seq_rows
    if not occ_rows or not seq_rows:
        sys.exit(f"Need both occurrence and sequence records; got "
                 f"{len(occ_rows)} and {len(seq_rows)} from {path}")

    verdicts = load_verdicts()
    anchor_counts = load_anchor_counts()
    accuracy, checkable, contradicted = label_accuracy(verdicts)
    scores, var, _ = pca2([r["_bio"] for r in rows])
    occ_mask = np.array([r["_role"] == "occurrence" for r in rows])
    seq_mask = np.array([r["_role"] == "sequence" for r in rows])
    occ_scores, seq_scores = scores[occ_mask], scores[seq_mask]

    pad = 0.05 * (scores.max(axis=0) - scores.min(axis=0))
    lo, hi = scores.min(axis=0) - pad, scores.max(axis=0) + pad
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], GRID),
                         np.linspace(lo[1], hi[1], GRID))

    seq_cov = reference_covariance(seq_scores)
    occ_density, _ = density_on_grid(occ_scores, gx, gy)
    naive_density, naive_area = density_on_grid(seq_scores, gx, gy,
                                                covariance=seq_cov)
    w, unmatched = weights_for(seq_rows, verdicts, anchor_counts=anchor_counts, registry=registry)
    weighted_density, wtd_area = density_on_grid(seq_scores, gx, gy, weights=w, covariance=seq_cov)

    naive_cov, naive_def, naive_index = coverage_and_deficit(
        occ_density, naive_density, len(seq_scores), naive_area)
    wtd_cov, wtd_def, wtd_index = coverage_and_deficit(
        occ_density, weighted_density, w.sum(), wtd_area)

    # Bounds on the index. The first two vary only what an unverified record is
    # worth. The third is the floor: every weight replaced by the chance the
    # label is simply correct, which is the strongest objection to the method.
    bounds = {}
    for label, unmatched_w in (("unverified trusted", 1.0),
                               ("unverified discarded", 0.0)):
        bw, _ = weights_for(seq_rows, verdicts, WEIGHTS, unmatched_w,
                            anchor_counts, registry)
        bd, ba = density_on_grid(seq_scores, gx, gy, weights=bw, covariance=seq_cov)
        bounds[label] = coverage_and_deficit(occ_density, bd, bw.sum(), ba)[2]
    if accuracy is not None:
        # Every weight set to the non-contradiction rate, as described: the
        # scheme in which verifiability is ignored and only the chance that a
        # label is uncontradicted is credited.
        acc_scheme = {k: accuracy for k in WEIGHTS}
        bw, _ = weights_for(seq_rows, verdicts, acc_scheme, accuracy, anchor_counts, registry)
        bd, ba = density_on_grid(seq_scores, gx, gy, weights=bw, covariance=seq_cov)
        bounds["non-contradiction"] = coverage_and_deficit(
            occ_density, bd, bw.sum(), ba)[2]

    # Per-occurrence-record uncovered fraction under the weighted surface.
    wtd_cov_matrix = seq_cov
    uncovered = 1.0 - sample_surface(wtd_cov, gx, gy, occ_scores[:, 0], occ_scores[:, 1])
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

    # Rank every occurrence-bearing cell by the exact reduction in the deficit
    # index that one verified specimen from that cell's climate would produce.
    lam = w.sum() * weighted_density * wtd_area
    cells = {}
    for row, pc_pos, u in zip(occ_rows, occ_scores, uncovered):
        key = (int(np.floor(row["_lat"] / CELL) * CELL),
               int(np.floor(row["_lon"] / CELL) * CELL))
        c = cells.setdefault(key, {"n": 0, "u": [], "pc": [], "country": {}})
        c["n"] += 1
        c["u"].append(float(u))
        c["pc"].append(pc_pos)
        country = (row.get("country") or "").strip()
        if country:
            c["country"][country] = c["country"].get(country, 0) + 1

    scored = []
    for key, c in cells.items():
        centroid = np.array(c["pc"]).mean(axis=0)
        scored.append({
            "lat_bin": key[0], "lon_bin": key[1],
            "lat_range": f"{key[0]} to {key[0] + CELL:g}",
            "lon_range": f"{key[1]} to {key[1] + CELL:g}",
            "country": max(c["country"], key=c["country"].get) if c["country"] else "",
            "occurrence_records": c["n"],
            "uncovered_mass": round(float(np.sum(c["u"])), 4),
            "mean_uncovered_fraction": round(float(np.mean(c["u"])), 4),
            "deficit_reduction": round(marginal_gain(
                occ_density, lam, gx, gy, centroid, wtd_cov_matrix,
                MIN_RELIABLE), 6),
        })
    scored.sort(key=lambda r: -r["deficit_reduction"])
    picks = scored[:N_PRIORITY]
    with open(os.path.join(outdir, "priority_list.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["rank", "lat_bin", "lon_bin", "lat_range",
                                           "lon_range", "country",
                                           "occurrence_records", "uncovered_mass",
                                           "mean_uncovered_fraction",
                                           "deficit_reduction"])
        wr.writeheader()
        for i, p in enumerate(picks, 1):
            wr.writerow({"rank": i, **p})

    sens_rows, thresh_rows = [], []
    if sensitivity:
        for name, (scheme, unmatched_w) in SCHEMES.items():
            if scheme is None:            # the label_accuracy floor
                if accuracy is None:
                    continue
                scheme, unmatched_w = acc_scheme, accuracy
            sw, _ = weights_for(seq_rows, verdicts, scheme, unmatched_w,
                                anchor_counts, registry)
            d, area = density_on_grid(seq_scores, gx, gy, weights=sw, covariance=seq_cov)
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
        f"records used             {len(rows)} (occurrence {len(occ_rows)}, "
        f"sequence specimens {len(seq_rows)})",
        f"duplicates collapsed     {collapse_report['collapsed']} rows "
        f"({collapse_report['coordinate_conflicts']} with conflicting coordinates)",
        f"crosswalk aliases        {registry.linked}",
        f"moved to occurrence      {reclass_note}",
        f"verdicts matched         {len(seq_rows) - unmatched}/{len(seq_rows)}"
        f" ({unmatched} took the default weight of {DEFAULT_WEIGHT})",
        f"mean verifiability       {w.mean():.3f}",
        f"effective sequences      {w.sum():.1f} of {len(seq_rows)} records",
        f"coverage threshold       {MIN_RELIABLE:g} verifiable sequences per neighbourhood",
        f"shared-cluster weight    derived as 1/k from the anchored species per "
        f"cluster ({len(anchor_counts)} records)",
        f"PCA variance explained   PC1 {var[0]:.1%}, PC2 {var[1]:.1%}",
        "",
        f"non-contradiction rate  {accuracy:.3f} of {checkable} checkable labels "
        f"were not contradicted ({contradicted} were)" if accuracy is not None else "",
        "",
        f"deficit index, naive     {naive_index:.3f}",
        f"deficit index, weighted  {wtd_index:.3f}",
        f"inflation of the gap     {wtd_index - naive_index:+.3f}"
        + (f" ({(wtd_index / naive_index - 1) * 100:+.1f}% relative)"
           if naive_index else ""),
        "bounds:",
        f"  {bounds['unverified trusted']:.3f} unverified records trusted in full",
        f"  {bounds['unverified discarded']:.3f} unverified records discarded",
    ] + ([
        f"  {bounds['non-contradiction']:.3f} FLOOR: every weight set to the "
        f"non-contradiction rate",
    ] if "labels merely correct" in bounds else []) + [
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
        lines += ["", f"top {min(5, len(picks))} priorities by exact deficit reduction:"]
        lines += [f"  {i}. {p['country'] or 'unknown'} "
                  f"({p['lat_range']}, {p['lon_range']})  "
                  f"deficit reduction {p['deficit_reduction']:.5f}"
                  for i, p in enumerate(picks[:5], 1)]

    text = "\n".join(l for l in lines if l is not None)
    with open(os.path.join(outdir, "deficit_summary.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
