"""Figures."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# One-hue sequential ramp: light = no signal, dark = strong signal.
CMAP = LinearSegmentedColormap.from_list(
    "seq_blue", ["#f7f9fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                 "#256abf", "#184f95", "#0d366b"])
INK, MUTED, RULE = "#1f1f1e", "#6b6a66", "#8a8984"

plt.rcParams.update({"font.size": 9, "axes.titleweight": "bold", "axes.titlesize": 11})


def heatmap(res: pd.DataFrame, out_path, classes: pd.Series | None = None,
            title="Which brain cell types are enriched for disease risk genes?"):
    """Cell types (rows, grouped by class) x traits. Colour = -log10 P; * FDR<5%, ** FDR<1%."""
    p = res.pivot(index="cell_type", columns="trait", values="p")
    fdr = res.pivot(index="cell_type", columns="trait", values="fdr")
    traits = [t for t in res["trait"].drop_duplicates() if "control" not in t.lower()] + \
             [t for t in res["trait"].drop_duplicates() if "control" in t.lower()]
    p, fdr = p[traits], fdr[traits]

    if classes is not None:
        cls = classes.reindex(p.index).fillna("Other")
        first = {"neuronal cells": 0, "glial cells": 1}
        order = sorted(p.index, key=lambda ct: (first.get(cls[ct].lower(), 2), cls[ct], ct))
        cls = cls[order]
    else:
        order, cls = sorted(p.index), None
    p, fdr = p.loc[order], fdr.loc[order]
    logp = -np.log10(p)
    vmax = max(4.0, float(np.ceil(np.nanmax(logp.to_numpy()))))

    nr, nc = p.shape
    # Fixed physical sizes (inches) so the figure looks the same for any data size.
    left, cell_w, cell_h, right, top, bottom = 2.7, 0.55, 0.24, 2.0, 1.7, 1.0
    W, H = left + cell_w * nc + right, top + cell_h * nr + bottom
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([left / W, bottom / H, cell_w * nc / W, cell_h * nr / H])
    im = ax.imshow(logp.to_numpy(), cmap=CMAP, vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(-.5, nc), minor=True)
    ax.set_yticks(np.arange(-.5, nr), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    for i in range(nr):
        for j in range(nc):
            q = fdr.iat[i, j]
            mark = "**" if q < 0.01 else ("*" if q < 0.05 else "")
            if mark:
                ax.text(j, i, mark, ha="center", va="center", fontsize=10, fontweight="bold",
                        color="white" if logp.iat[i, j] / vmax > 0.55 else INK)

    # Separate and label the cell classes on the right-hand side.
    if cls is not None:
        bounds = np.flatnonzero(cls.to_numpy()[1:] != cls.to_numpy()[:-1]) + 1
        for b in bounds:
            ax.axhline(b - 0.5, color=RULE, linewidth=1)
        starts = np.r_[0, bounds]
        ends = np.r_[bounds, nr]
        for s, e in zip(starts, ends):
            ax.text(nc - 0.35, (s + e - 1) / 2, cls.iloc[s], ha="left", va="center",
                    fontsize=8.5, color=MUTED, fontstyle="italic")

    ax.xaxis.tick_top()
    ax.set_xticks(range(nc), p.columns, rotation=40, ha="left", rotation_mode="anchor")
    ax.set_yticks(range(nr), p.index)
    fig.text(0.15 / W, 1 - 0.3 / H, title, fontsize=12, fontweight="bold", color=INK,
             va="top")

    cax = fig.add_axes([left / W, 0.45 / H, max(2.4, cell_w * nc) / W, 0.12 / H])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("−log10(P), one-sided Fisher's exact test", color=MUTED, fontsize=8)
    cb.outline.set_visible(False)
    cax.tick_params(length=0, labelsize=8, colors=MUTED)
    fig.text((left + max(2.4, cell_w * nc) + 0.3) / W, 0.51 / H,
             "* FDR < 5%\n** FDR < 1%", fontsize=8, color=MUTED, va="center")
    fig.savefig(out_path, dpi=200, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
