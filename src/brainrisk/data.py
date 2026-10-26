"""Download and load the two public datasets.

1. Human Protein Atlas (v24) single-nucleus brain atlas: average expression (nTPM)
   of every gene in 34 brain cell types, built from ~2.5 million nuclei across 11
   brain regions (Siletti et al. 2023, Science). About 10 MB.
2. NHGRI-EBI GWAS Catalog: every published genome-wide significant association,
   with the genes each variant maps to. One-time download.
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import pandas as pd

HPA_URL = "https://v24.proteinatlas.org/download/tsv/rna_single_nuclei_cluster_type.tsv.zip"
HPA_CLASSES_URL = ("https://v24.proteinatlas.org/download/tsv/"
                   "rna_single_nuclei_cluster_type_cluster_types.tsv.zip")
GWAS_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
            "gwas-catalog-associations_ontology-annotated-full.zip")

GWAS_COLUMNS = ["CHR_ID", "CHR_POS", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "P-VALUE", "SNPS", "SNP_GENE_IDS",
                "UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", "UPSTREAM_GENE_DISTANCE",
                "DOWNSTREAM_GENE_DISTANCE", "MAPPED_GENE", "STUDY ACCESSION", "PUBMEDID"]


def download(url: str, dest: Path) -> Path:
    """Download `url` to `dest` unless it is already there."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"Downloading {url}")

    def progress(blocks, block_size, total):
        done = blocks * block_size
        if total > 0:
            sys.stdout.write(f"\r  {min(done, total) / 1e6:6.1f} / {total / 1e6:.1f} MB")
        else:
            sys.stdout.write(f"\r  {done / 1e6:6.1f} MB")
        sys.stdout.flush()

    opener = urllib.request.build_opener()
    opener.addheaders = [("User-Agent", "brain-risk-celltypes (research script)")]
    urllib.request.install_opener(opener)
    urllib.request.urlretrieve(url, tmp, reporthook=progress)
    print()
    tmp.rename(dest)
    return dest


def load_hpa_expression(path) -> pd.DataFrame:
    """Genes (Ensembl ID) x cell types matrix of nTPM."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    needed = {"Gene", "Cluster type", "nTPM"}
    if not needed <= set(df.columns):
        raise ValueError(f"Unexpected HPA columns: {list(df.columns)}")
    mat = df.pivot_table(index="Gene", columns="Cluster type", values="nTPM", aggfunc="mean")
    return mat.fillna(0.0)


def load_hpa_gene_names(path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", usecols=["Gene", "Gene name"]).drop_duplicates("Gene")
    return df.set_index("Gene")["Gene name"]


def load_hpa_classes(path) -> pd.Series:
    """Cell type -> broad class (neuronal, glial, vascular, ...)."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    type_col = next(c for c in df.columns if c.lower() == "cluster type")
    class_col = next(c for c in df.columns if c != type_col)
    return df.set_index(type_col)[class_col]


def load_gwas_catalog(path) -> pd.DataFrame:
    """Only the columns we need, read as text (the full file is large)."""
    return pd.read_csv(path, sep="\t", usecols=lambda c: c in GWAS_COLUMNS, dtype=str,
                       low_memory=False)
