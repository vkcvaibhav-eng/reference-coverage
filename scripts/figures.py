"""Manuscript figures.

Usage:
    python figures.py [--genus NAME]

Reads results/deficit[/<genus>]/ and results/holdout[/<genus>]/ and writes
figures to results/figures[/<genus>]/ as PDF and 600 dpi PNG.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "savefig.bbox": "tight",
})
CM = 1 / 2.54
SINGLE, DOUBLE = 8.5 * CM, 17.5 * CM


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=600)
    plt.close(fig)
    print(f"  {name}")


def fig_climate_space(dd, outdir):
    """Occurrence density, the two coverage surfaces, and the deficit."""
    rows = read(os.path.join(dd, "climate_space.csv"))
    n = int(np.sqrt(len(rows)))
    g = lambda k: np.array([float(r[k]) for r in rows]).reshape(n, n)
    pc1, pc2 = g("pc1"), g("pc2")
    ext = [pc1.min(), pc1.max(), pc2.min(), pc2.max()]

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, DOUBLE / 2.9))
    fig.subplots_adjust(wspace=0.55)
    panels = [
        ("occurrence_density", "Occurrence density", "Greys", None),
        ("coverage_weighted", "Verifiable coverage", "viridis", (0, 1)),
        ("deficit_weighted", "Coverage deficit", "magma", None),
    ]
    for ax, (key, title, cmap, lim) in zip(axes, panels):
        z = g(key)
        kw = {"vmin": lim[0], "vmax": lim[1]} if lim else {}
        im = ax.imshow(z, origin="lower", extent=ext, aspect="auto",
                       cmap=cmap, **kw)
        ax.contour(pc1, pc2, g("occurrence_density"), levels=4,
                   colors="white", linewidths=0.4, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel("PC1 of BIO1–BIO19")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.06)
        cb.ax.tick_params(labelsize=6, width=0.5, pad=1.5)
        cb.outline.set_linewidth(0.5)
    axes[0].set_ylabel("PC2 of BIO1–BIO19")
    for ax in axes[1:]:                      # one y axis is enough for three panels
        ax.tick_params(labelleft=False)
    for ax, lab in zip(axes, "abc"):
        ax.text(0.02, 0.97, f"({lab})", transform=ax.transAxes, va="top",
                fontsize=8, fontweight="bold", color="black",
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.2))
    save(fig, outdir, "Fig1_climate_space")


def fig_map(dd, outdir):
    """Geographic deficit, with the sampling priorities marked."""
    cells = read(os.path.join(dd, "geographic_deficit.csv"))
    pri = read(os.path.join(dd, "priority_list.csv"))

    lat = np.array([float(c["lat_bin"]) for c in cells])
    lon = np.array([float(c["lon_bin"]) for c in cells])
    unc = np.array([float(c["mean_uncovered_fraction"]) for c in cells])
    nrec = np.array([float(c["occurrence_records"]) for c in cells])

    fig, ax = plt.subplots(figsize=(DOUBLE, DOUBLE / 2.3))
    sc = ax.scatter(lon + 1, lat + 1, c=unc, s=8 + 22 * np.sqrt(nrec / nrec.max()),
                    cmap="magma_r", vmin=0, vmax=1, linewidths=0.2,
                    edgecolors="0.3")
    top = pri[:10]
    ax.scatter([float(p["lon_bin"]) + 1 for p in top],
               [float(p["lat_bin"]) + 1 for p in top],
               s=70, facecolors="none", edgecolors="#1f78b4", linewidths=1.0)
    for i, p in enumerate(top[:5], 1):
        ax.annotate(str(i), (float(p["lon_bin"]) + 1, float(p["lat_bin"]) + 1),
                    textcoords="offset points", xytext=(6, 4), fontsize=7,
                    color="#1f78b4", fontweight="bold")
    ax.axhline(0, color="0.8", lw=0.4)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 80)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-60, 81, 30))
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("Mean uncovered fraction", fontsize=7)
    cb.ax.tick_params(labelsize=6, width=0.5)
    cb.outline.set_linewidth(0.5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", mfc="none",
                              mec="#1f78b4", ms=7, label="Top 10 priorities")],
              loc="lower left", frameon=False)
    save(fig, outdir, "Fig2_geographic_deficit")


def fig_sensitivity(dd, outdir):
    """Weight schemes and the coverage threshold."""
    ws = read(os.path.join(dd, "weight_sensitivity.csv"))
    ts = read(os.path.join(dd, "threshold_sensitivity.csv"))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE / 2.7))

    order = ["naive", "label_accuracy", "unverified_trusted", "lenient",
             "default", "strict", "anchored_only"]
    ws = sorted(ws, key=lambda r: order.index(r["scheme"])
                if r["scheme"] in order else 99)
    labels = [r["scheme"].replace("_", " ") for r in ws]
    idx = [float(r["deficit_index"]) for r in ws]
    eff = [float(r["effective_records"]) for r in ws]
    naive = next(float(r["deficit_index"]) for r in ws if r["scheme"] == "naive")

    y = np.arange(len(ws))
    cols = ["0.55" if r["scheme"] == "naive" else "#2b5d8a" for r in ws]
    a1.barh(y, idx, color=cols, height=0.62)
    a1.axvline(naive, color="0.25", ls="--", lw=0.7)
    a1.set_yticks(y)
    a1.set_yticklabels(labels)
    a1.invert_yaxis()
    a1.set_xlabel("Coverage deficit index")
    for yi, (v, e) in enumerate(zip(idx, eff)):
        a1.text(v + 0.006, yi, f"{v:.3f}  (n={e:.0f})", va="center", fontsize=6)
    a1.set_xlim(0, max(idx) * 1.33)

    k = [float(r["min_reliable"]) for r in ts]
    dn = [float(r["deficit_naive"]) for r in ts]
    dw = [float(r["deficit_weighted"]) for r in ts]
    a2.plot(k, dn, "o-", color="0.45", ms=3.5, lw=1.0, label="Record count")
    a2.plot(k, dw, "s-", color="#2b5d8a", ms=3.5, lw=1.0,
            label="Verifiability weighted")
    a2.fill_between(k, dn, dw, color="#2b5d8a", alpha=0.12)
    for xi, a, b in zip(k, dn, dw):
        a2.annotate(f"+{(b / a - 1) * 100:.0f}%", (xi, b),
                    textcoords="offset points", xytext=(0, 5),
                    ha="center", fontsize=6, color="#2b5d8a")
    a2.set_xscale("log")
    a2.set_xticks(k)
    a2.set_xticklabels([f"{int(v)}" for v in k])
    a2.set_xlabel("Verifiable sequences per neighbourhood ($k$)")
    a2.set_ylabel("Coverage deficit index")
    a2.legend(frameon=False, loc="upper left")
    a2.set_ylim(0, max(dw) * 1.2)

    for ax, lab in ((a1, "a"), (a2, "b")):
        ax.text(-0.12, 1.08, f"({lab})", transform=ax.transAxes, va="top",
                fontsize=8, fontweight="bold")
    save(fig, outdir, "Fig3_sensitivity")


def fig_holdout(hd, outdir, min_train=150):
    """Observed uncovered fraction against both nulls, by cutoff year."""
    rows = [r for r in read(os.path.join(hd, "cutoff_sweep.csv")) if r["observed"]]
    yr = np.array([int(r["cutoff"]) for r in rows])
    obs = np.array([float(r["observed"]) for r in rows])
    n_occ = np.array([float(r["null_occurrence"]) for r in rows])
    n_seq = np.array([float(r["null_sequencing"]) for r in rows])
    thin = np.array([bool(r["note"]) for r in rows])
    ntr = np.array([int(r["train_sequence"]) for r in rows])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(DOUBLE, DOUBLE / 2.7),
                                 gridspec_kw={"width_ratios": [1.35, 1]})

    a1.plot(yr, n_occ, "^--", color="#b2182b", ms=4, lw=1.0,
            label="Null: recorded range")
    a1.plot(yr, n_seq, "v--", color="#4d9221", ms=4, lw=1.0,
            label="Null: earlier reference records")
    a1.plot(yr[~thin], obs[~thin], "o-", color="#2b5d8a", ms=5, lw=1.3,
            label="Observed")
    a1.plot(yr[thin], obs[thin], "o-", color="#2b5d8a", ms=5, lw=1.3,
            mfc="white", alpha=0.55)
    if thin.any():
        a1.axvspan(yr.min() - 0.6, yr[thin].max() + 0.6, color="0.9", zorder=0)
        a1.text(yr[thin].mean(), 0.02, f"training library\n< {min_train} sequences",
                ha="center", fontsize=6, color="0.35")
    a1.set_xticks(yr)
    a1.set_xticklabels([str(v) for v in yr], rotation=45, ha="right")
    a1.set_xlabel("Cutoff year")
    a1.set_ylabel("Mean uncovered fraction")
    a1.set_ylim(0, 1.0)
    a1.legend(frameon=False, loc="upper right")

    ratio = obs / n_occ
    a2.axhline(1.0, color="0.25", ls="--", lw=0.8)
    a2.plot(yr[~thin], ratio[~thin], "o-", color="#2b5d8a", ms=5, lw=1.3)
    a2.plot(yr[thin], ratio[thin], "o", color="#2b5d8a", ms=5, mfc="white",
            alpha=0.55)
    a2.set_xticks(yr)
    a2.set_xticklabels([str(v) for v in yr], rotation=45, ha="right")
    a2.set_xlabel("Cutoff year")
    a2.set_ylabel("Observed / null (recorded range)")
    a2.set_ylim(0, max(1.35, ratio.max() * 1.12))
    a2.text(0.5, 0.93, "later-collected records in less-covered climate",
            transform=a2.transAxes, ha="center", fontsize=6, color="0.35")
    a2.text(0.5, 0.06, "later-collected records in better-covered climate",
            transform=a2.transAxes, ha="center", fontsize=6, color="0.35")
    for x, y_, n in zip(yr, ratio, ntr):
        a2.annotate(f"{n}", (x, y_), textcoords="offset points",
                    xytext=(0, -11), ha="center", fontsize=5.5, color="0.45")

    for ax, lab in ((a1, "a"), (a2, "b")):
        ax.text(-0.13, 1.08, f"({lab})", transform=ax.transAxes, va="top",
                fontsize=8, fontweight="bold")
    save(fig, outdir, "Fig4_holdout")


def main():
    genus = None
    if "--genus" in sys.argv:
        genus = sys.argv[sys.argv.index("--genus") + 1].lower()
    dd = os.path.join("results", "deficit", genus) if genus \
        else os.path.join("results", "deficit")
    hd = os.path.join("results", "holdout", genus) if genus \
        else os.path.join("results", "holdout")
    outdir = os.path.join("results", "figures", genus) if genus \
        else os.path.join("results", "figures")
    print(f"figures from {dd} and {hd}:")
    fig_climate_space(dd, outdir)
    fig_map(dd, outdir)
    fig_sensitivity(dd, outdir)
    if os.path.exists(os.path.join(hd, "cutoff_sweep.csv")):
        fig_holdout(hd, outdir)
    print(f"Saved to {outdir}/")


if __name__ == "__main__":
    main()
