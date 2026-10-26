import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from brainrisk import analysis, risk_genes, synthetic  # noqa: E402
import run  # noqa: E402


# ------------------------------------------------------------ specificity

def test_specificity_rows_sum_to_one():
    expr = pd.DataFrame({"A": [10, 0, 5, 0.1], "B": [0, 10, 5, 0.2]}, index=list("wxyz"))
    spec = analysis.specificity(expr)
    assert "z" not in spec.index                     # never reaches 1 nTPM: dropped
    assert np.allclose(spec.sum(axis=1), 1)
    assert spec.loc["w", "A"] == 1 and spec.loc["y", "A"] == 0.5


def test_signatures_take_top_fraction():
    rng = np.random.default_rng(0)
    spec = pd.DataFrame(rng.random((300, 4)), columns=list("ABCD"))
    assert (analysis.signatures(spec, 0.1).sum() == 30).all()


def test_bh_fdr_matches_r():
    # R: p.adjust(c(0.01, 0.04, 0.03, 0.2), "BH")
    assert np.allclose(analysis.bh_fdr([0.01, 0.04, 0.03, 0.2]), [0.04, 0.16 / 3, 0.16 / 3, 0.2])


def test_fisher_null_is_calibrated():
    """Random risk genes should give P < 0.05 about 5% of the time (or less: Fisher is conservative)."""
    rng = np.random.default_rng(1)
    genes = pd.Index([f"g{i}" for i in range(3000)])
    sigs = pd.DataFrame({f"c{k}": rng.random(3000) < 0.1 for k in range(10)}, index=genes)
    risk = {f"t{k}": set(rng.choice(genes, 60, replace=False)) for k in range(40)}
    res = analysis.enrichment(sigs, risk)
    assert (res["p"] < 0.05).mean() < 0.08


# ------------------------------------------------------------ GWAS parsing

def row(**kw):
    base = {"SNP_GENE_IDS": "", "UPSTREAM_GENE_ID": "", "DOWNSTREAM_GENE_ID": "",
            "UPSTREAM_GENE_DISTANCE": "", "DOWNSTREAM_GENE_DISTANCE": "",
            "CHR_ID": "1", "CHR_POS": "100"}
    return {**base, **kw}


def test_gene_mapping_inside_and_nearest():
    assert risk_genes.genes_for_row(row(SNP_GENE_IDS="ENSG01, ENSG02")) == ["ENSG01", "ENSG02"]
    near = row(UPSTREAM_GENE_ID="ENSG10", DOWNSTREAM_GENE_ID="ENSG20",
               UPSTREAM_GENE_DISTANCE="9000", DOWNSTREAM_GENE_DISTANCE="300")
    assert risk_genes.genes_for_row(near) == ["ENSG20"]


def test_mhc_detection():
    assert risk_genes.in_mhc(row(CHR_ID="6", CHR_POS="29000000"))
    assert not risk_genes.in_mhc(row(CHR_ID="6", CHR_POS="50000000"))
    assert not risk_genes.in_mhc(row(CHR_ID="1", CHR_POS="29000000"))


def test_trait_matching_by_id_or_label():
    spec = risk_genes.TRAITS["Alzheimer's disease"]
    assert risk_genes.matches_trait("http://purl.obolibrary.org/obo/MONDO_0004975", "x", spec)
    assert risk_genes.matches_trait("http://www.ebi.ac.uk/efo/EFO_9", "Alzheimer disease", spec)
    assert not risk_genes.matches_trait("http://www.ebi.ac.uk/efo/EFO_9",
                                        "family history of Alzheimer disease", spec)


# ------------------------------------------------------------ end to end

def test_end_to_end_recovers_planted_signal(tmp_path):
    synthetic.write_fake_data(tmp_path / "data", seed=5)
    res = run.run(tmp_path / "data", tmp_path / "out", download=False)
    hits = set(map(tuple, res.loc[res["fdr"] < 0.05, ["trait", "cell_type"]].to_numpy()))
    assert hits == {("Alzheimer's disease", "Central nervous system macrophage"),
                    ("Schizophrenia", "Upper-layer intratelencephalic")}
    # decoy rows (non-significant, MHC, other traits) all pointed at pericytes
    assert not res.loc[res["cell_type"] == "Pericyte", "fdr"].lt(0.05).any()
    for f in ["heatmap.png", "enrichment.tsv", "top_hits.md", "risk_genes.tsv"]:
        assert (tmp_path / "out" / f).stat().st_size > 0
