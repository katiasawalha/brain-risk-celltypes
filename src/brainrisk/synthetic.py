"""Fake HPA and GWAS Catalog files with a known answer, in the real file formats.

Planted truth:
  * "Alzheimer's disease" risk genes are enriched in microglia-specific genes
  * "Schizophrenia" risk genes are enriched in excitatory-neuron-specific genes
  * "Height (control)" has random risk genes (no cell type should be enriched)
The file also contains rows that must be ignored: non-significant hits, hits in
the MHC region, and unrelated traits.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

CELL_TYPES = {
    "Upper-layer intratelencephalic": "Neuronal cells",
    "Deep-layer intratelencephalic": "Neuronal cells",
    "MGE interneuron": "Neuronal cells",
    "CGE interneuron": "Neuronal cells",
    "Medium spiny neuron": "Neuronal cells",
    "Astrocyte": "Glial cells",
    "Oligodendrocyte": "Glial cells",
    "Oligodendrocyte precursor cell": "Glial cells",
    "Central nervous system macrophage": "Glial cells",
    "Endothelial cell": "Endothelial and mural cells",
    "Pericyte": "Endothelial and mural cells",
}
PLANTED = {"Alzheimer's disease": ("MONDO_0004975", "Alzheimer disease",
                                   "Central nervous system macrophage"),
           "Schizophrenia": ("EFO_0000692", "schizophrenia", "Upper-layer intratelencephalic"),
           "Height (control)": ("EFO_0004339", "body height", None)}


def _zip_tsv(df: pd.DataFrame, path: Path, inner: str):
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(inner, df.to_csv(sep="\t", index=False))


def write_fake_data(data_dir, n_genes=6000, n_markers=300, n_risk=80, seed=0):
    rng = np.random.default_rng(seed)
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    genes = np.array([f"ENSG{i:011d}" for i in range(n_genes)])
    names = np.array([f"GENE{i}" for i in range(n_genes)])
    cts = list(CELL_TYPES)

    # --- expression -------------------------------------------------------
    base = rng.lognormal(1.5, 1.5, n_genes)
    perm = rng.permutation(n_genes)
    markers = {}
    rows = []
    for k, ct in enumerate(cts):
        idx = perm[k * n_markers:(k + 1) * n_markers]
        markers[ct] = idx
        expr = base * rng.lognormal(0, 0.4, n_genes)
        expr[idx] *= rng.uniform(5, 30, len(idx))
        rows.append(pd.DataFrame({"Gene": genes, "Gene name": names, "Cluster type": ct,
                                  "nTPM": expr.round(1)}))
    _zip_tsv(pd.concat(rows), data_dir / "rna_single_nuclei_cluster_type.tsv.zip",
             "rna_single_nuclei_cluster_type.tsv")
    _zip_tsv(pd.DataFrame({"Cluster type": cts, "Cell type group": list(CELL_TYPES.values())}),
             data_dir / "rna_single_nuclei_cluster_type_cluster_types.tsv.zip",
             "rna_single_nuclei_cluster_type_cluster_types.tsv")

    # --- GWAS catalog -----------------------------------------------------
    recs = []

    def add(trait_id, label, gene_idx, p, chrom="1", pos=1_000_000, intergenic=False):
        g = genes[gene_idx]
        recs.append({
            "CHR_ID": chrom, "CHR_POS": str(pos), "MAPPED_GENE": names[gene_idx],
            "SNP_GENE_IDS": "" if intergenic else g,
            "UPSTREAM_GENE_ID": g if intergenic else "",
            "DOWNSTREAM_GENE_ID": genes[(gene_idx + 1) % n_genes] if intergenic else "",
            "UPSTREAM_GENE_DISTANCE": "500" if intergenic else "",
            "DOWNSTREAM_GENE_DISTANCE": "90000" if intergenic else "",
            "SNPS": f"rs{len(recs) + 1}", "P-VALUE": f"{p:.0E}",
            "MAPPED_TRAIT": label, "MAPPED_TRAIT_URI": f"http://www.ebi.ac.uk/efo/{trait_id}",
            "STUDY ACCESSION": "GCST000001", "PUBMEDID": "1",
        })

    for trait, (tid, label, ct) in PLANTED.items():
        if ct is None:
            chosen = rng.choice(n_genes, n_risk, replace=False)
        else:
            half = n_risk // 2
            chosen = np.r_[rng.choice(markers[ct], half, replace=False),
                           rng.choice(n_genes, n_risk - half, replace=False)]
        for j, gi in enumerate(chosen):
            add(tid, label, gi, p=1e-10, intergenic=(j % 4 == 0))
        # decoys that must be ignored
        mg = markers[cts[-1]]  # pericyte markers
        for gi in mg[:40]:
            add(tid, label, gi, p=1e-3)                                  # not significant
            add(tid, label, gi, p=1e-12, chrom="6", pos=30_000_000)       # MHC
    for gi in markers[cts[-1]][:60]:
        add("EFO_0000001", "some other trait", gi, p=1e-20)              # unrelated trait

    _zip_tsv(pd.DataFrame(recs), data_dir / "gwas-catalog-associations_ontology-annotated-full.zip",
             "gwas-catalog-associations_ontology-annotated-full.tsv")
    return markers
