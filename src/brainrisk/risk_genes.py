"""Turn GWAS Catalog associations into a list of risk genes per trait.

A trait's risk genes are the genes mapped to its genome-wide significant variants
(P < 5e-8). A variant inside a gene maps to that gene; an intergenic variant maps
to the nearest gene on either side. This is the GWAS Catalog's own mapping
(Ensembl IDs), so no gene-name matching is needed.

Variants in the MHC region (chr6:25-34 Mb) are skipped: it holds hundreds of
immune genes in strong linkage, and one schizophrenia signal there would add
dozens of immune genes and bias the result toward microglia.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

GENOME_WIDE = 5e-8
MHC = ("6", 25_000_000, 34_000_000)

# Each trait is matched by ontology ID (EFO / MONDO) or by its exact catalog label.
# Both are listed because studies of the same disease are sometimes annotated
# with different but equivalent terms.
TRAITS = {
    "Schizophrenia": {"ids": ["EFO_0000692", "MONDO_0005090"],
                      "labels": ["schizophrenia"]},
    "Bipolar disorder": {"ids": ["EFO_0000289", "MONDO_0004985"],
                         "labels": ["bipolar disorder"]},
    "Major depression": {"ids": ["MONDO_0002009", "EFO_0003761"],
                         "labels": ["major depressive disorder", "unipolar depression"]},
    "Autism": {"ids": ["EFO_0003756"],
               "labels": ["autism spectrum disorder", "autism"]},
    "Alzheimer's disease": {"ids": ["MONDO_0004975"],
                            "labels": ["alzheimer disease", "alzheimer's disease"]},
    "Parkinson's disease": {"ids": ["MONDO_0005180"],
                            "labels": ["parkinson disease", "parkinson's disease"]},
    # Negative control: a non-brain trait. No brain cell type should stand out.
    "Height (control)": {"ids": ["EFO_0004339"],
                         "labels": ["body height", "height"]},
}

_ENSG = re.compile(r"ENSG\d+")


def _split(value) -> list[str]:
    if not isinstance(value, str):
        return []
    return [v.strip() for v in re.split(r",", value) if v.strip()]


def matches_trait(uris: str, labels: str, spec: dict) -> bool:
    ids = {u.rsplit("/", 1)[-1] for u in _split(uris)}
    names = {l.lower() for l in _split(labels)}
    return bool(ids & set(spec["ids"]) or names & set(spec["labels"]))


def genes_for_row(row) -> list[str]:
    """Ensembl genes for one association: the gene(s) the variant is in, else the nearest."""
    inside = _ENSG.findall(str(row.get("SNP_GENE_IDS", "")))
    if inside:
        return inside
    up, down = _ENSG.findall(str(row.get("UPSTREAM_GENE_ID", ""))), \
        _ENSG.findall(str(row.get("DOWNSTREAM_GENE_ID", "")))
    d_up = pd.to_numeric(row.get("UPSTREAM_GENE_DISTANCE"), errors="coerce")
    d_down = pd.to_numeric(row.get("DOWNSTREAM_GENE_DISTANCE"), errors="coerce")
    if up and down:
        if np.isnan(d_up) or np.isnan(d_down):
            return up[:1] + down[:1]
        return up[:1] if d_up <= d_down else down[:1]
    return (up or down)[:1]


def in_mhc(row) -> bool:
    chrom = str(row.get("CHR_ID", "")).split(";")[0].strip()
    pos = pd.to_numeric(str(row.get("CHR_POS", "")).split(";")[0], errors="coerce")
    return chrom == MHC[0] and not np.isnan(pos) and MHC[1] <= pos <= MHC[2]


def risk_genes(catalog: pd.DataFrame, traits: dict = TRAITS,
               p_threshold: float = GENOME_WIDE, exclude_mhc: bool = True) -> dict[str, set]:
    """{trait: set of Ensembl gene IDs}."""
    p = pd.to_numeric(catalog["P-VALUE"], errors="coerce")
    sig = catalog[p < p_threshold]
    out = {}
    for trait, spec in traits.items():
        hit = sig[[matches_trait(u, l, spec) for u, l in
                   zip(sig["MAPPED_TRAIT_URI"], sig["MAPPED_TRAIT"])]]
        genes = set()
        for _, row in hit.iterrows():
            if exclude_mhc and in_mhc(row):
                continue
            genes.update(genes_for_row(row))
        out[trait] = genes
    return out
