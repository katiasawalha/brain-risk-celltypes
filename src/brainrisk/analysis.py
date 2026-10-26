"""Cell-type specificity and enrichment of disease risk genes.

Specificity (Skene et al. 2018): the share of a gene's total expression, summed
over all cell types, that falls in one cell type. 1.0 = expressed only there.
Each cell type's "signature" is its top 10% most specific genes.

Enrichment: for each trait and cell type, a one-sided Fisher's exact test asks
whether the trait's risk genes fall into the cell type's signature more often
than chance, using all expressed genes as the background.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def specificity(expr: pd.DataFrame, min_ntpm: float = 1.0) -> pd.DataFrame:
    """Genes x cell types; keeps genes reaching `min_ntpm` in at least one cell type."""
    expr = expr.loc[(expr >= min_ntpm).any(axis=1)]
    return expr.div(expr.sum(axis=1), axis=0)


def signatures(spec: pd.DataFrame, top_frac: float = 0.10) -> pd.DataFrame:
    """Boolean genes x cell types: top `top_frac` most specific genes of each cell type."""
    n = max(1, int(round(top_frac * len(spec))))
    out = pd.DataFrame(False, index=spec.index, columns=spec.columns)
    for ct in spec.columns:
        out.loc[spec[ct].nlargest(n).index, ct] = True
    return out


def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted P values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def enrichment(sigs: pd.DataFrame, risk: dict, gene_names: pd.Series | None = None,
               min_genes: int = 5) -> pd.DataFrame:
    """Fisher's exact test for every trait x cell type.

    Returns one row per test with the overlap size, odds ratio, P, FDR and the
    overlapping genes (named, if `gene_names` is given) so hits can be inspected.
    """
    background = sigs.index
    rows = []
    for trait, genes in risk.items():
        g = background.intersection(list(genes))
        if len(g) < min_genes:
            print(f"  skipping {trait}: only {len(g)} risk genes expressed in the brain data")
            continue
        is_risk = background.isin(g)
        for ct in sigs.columns:
            in_sig = sigs[ct].to_numpy()
            a = int(np.sum(is_risk & in_sig))          # risk, in signature
            b = int(np.sum(is_risk & ~in_sig))         # risk, not in signature
            c = int(np.sum(~is_risk & in_sig))         # not risk, in signature
            d = int(np.sum(~is_risk & ~in_sig))
            odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            overlap = background[is_risk & in_sig]
            if gene_names is not None:
                overlap = [gene_names.get(x, x) for x in overlap]
            rows.append({"trait": trait, "cell_type": ct, "n_risk_genes": len(g),
                         "overlap": a, "expected": round(len(g) * in_sig.mean(), 2),
                         "odds_ratio": odds, "p": p,
                         "genes": ",".join(sorted(map(str, overlap)))})
    res = pd.DataFrame(rows)
    res["fdr"] = bh_fdr(res["p"])
    return res.sort_values(["trait", "p"]).reset_index(drop=True)
