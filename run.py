"""Which brain cell types are enriched for disease risk genes?

    python run.py          # real data: downloads once, then runs in seconds
    python run.py --demo   # synthetic data with a known answer (no internet needed)

Outputs in results/ (or results/demo/):
    heatmap.png        main figure
    enrichment.tsv     every trait x cell type test, with the overlapping genes
    top_hits.md        significant results in a readable table
    risk_genes.tsv     the risk genes used for each trait
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from brainrisk import analysis, data, plotting, risk_genes, synthetic  # noqa: E402


def run(data_dir: Path, out_dir: Path, top_frac: float = 0.10, download: bool = True):
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "expr": data_dir / "rna_single_nuclei_cluster_type.tsv.zip",
        "classes": data_dir / "rna_single_nuclei_cluster_type_cluster_types.tsv.zip",
        "gwas": data_dir / "gwas-catalog-associations_ontology-annotated-full.zip",
    }
    if download:
        data.download(data.HPA_URL, files["expr"])
        data.download(data.HPA_CLASSES_URL, files["classes"])
        data.download(data.GWAS_URL, files["gwas"])

    print("Loading brain expression atlas ...")
    expr = data.load_hpa_expression(files["expr"])
    names = data.load_hpa_gene_names(files["expr"])
    classes = data.load_hpa_classes(files["classes"])
    print(f"  {expr.shape[0]:,} genes x {expr.shape[1]} cell types")

    spec = analysis.specificity(expr)
    sigs = analysis.signatures(spec, top_frac=top_frac)

    print("Loading GWAS Catalog ...")
    catalog = data.load_gwas_catalog(files["gwas"])
    risk = risk_genes.risk_genes(catalog)
    for trait, genes in risk.items():
        n_expr = len(sigs.index.intersection(list(genes)))
        print(f"  {trait:<22} {len(genes):5d} risk genes ({n_expr} expressed in brain data)")
    pd.DataFrame([(t, g, names.get(g, "")) for t, gs in risk.items() for g in sorted(gs)],
                 columns=["trait", "gene_id", "gene_name"]).to_csv(
        out_dir / "risk_genes.tsv", sep="\t", index=False)

    print("Testing enrichment ...")
    res = analysis.enrichment(sigs, risk, gene_names=names)
    res.to_csv(out_dir / "enrichment.tsv", sep="\t", index=False)
    plotting.heatmap(res, out_dir / "heatmap.png", classes=classes)
    write_top_hits(res, out_dir / "top_hits.md")

    hits = res[res["fdr"] < 0.05]
    print(f"\n{len(hits)} significant trait x cell type pairs (FDR < 5%):")
    if len(hits):
        print(hits[["trait", "cell_type", "overlap", "expected", "odds_ratio", "fdr"]]
              .to_string(index=False))
    print(f"\nResults written to {out_dir}")
    return res


def write_top_hits(res: pd.DataFrame, path: Path, max_genes: int = 12):
    hits = res[res["fdr"] < 0.05].sort_values("p")
    lines = ["| Trait | Cell type | Risk genes in signature (expected) | Odds ratio | FDR | Example genes |",
             "|---|---|---|---|---|---|"]
    for _, r in hits.iterrows():
        genes = r["genes"].split(",")
        shown = ", ".join(genes[:max_genes]) + (" ..." if len(genes) > max_genes else "")
        lines.append(f"| {r['trait']} | {r['cell_type']} | {r['overlap']} ({r['expected']:.1f}) "
                     f"| {r['odds_ratio']:.2f} | {r['fdr']:.1e} | {shown} |")
    if hits.empty:
        lines.append("| – | no significant results | | | | |")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true", help="run on synthetic data")
    ap.add_argument("--data-dir", type=Path, default=ROOT / "data")
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--top-frac", type=float, default=0.10)
    args = ap.parse_args()

    if args.demo:
        demo_data = ROOT / "data" / "demo"
        synthetic.write_fake_data(demo_data)
        run(demo_data, args.out_dir or ROOT / "results" / "demo", args.top_frac, download=False)
    else:
        run(args.data_dir, args.out_dir or ROOT / "results", args.top_frac)


if __name__ == "__main__":
    main()
