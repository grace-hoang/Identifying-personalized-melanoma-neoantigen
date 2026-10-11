#!/usr/bin/env python3
"""Salmon quant.sf -> transcript and gene TPM tables (for vcf-expression-annotator + ranking)."""
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant", required=True)
    ap.add_argument("--tx2gene", required=True)
    ap.add_argument("--out-tx", required=True)
    ap.add_argument("--out-gene", required=True)
    a = ap.parse_args()

    q = pd.read_csv(a.quant, sep="\t")
    q["transcript_id"] = q["Name"].astype(str).str.split("|").str[0]     # tolerate un-split GENCODE headers
    t2g = pd.read_csv(a.tx2gene, sep="\t").drop_duplicates("transcript_id")
    m = q.merge(t2g, on="transcript_id", how="left")
    unmapped = int(m["gene_id"].isna().sum())
    if unmapped:
        print(f"note: {unmapped} transcripts without gene mapping")
    m["gene_id"] = m["gene_id"].fillna("NA")
    m["gene_name"] = m["gene_name"].fillna("")
    tx = m[["transcript_id", "gene_id", "gene_name", "TPM", "NumReads"]].rename(
        columns={"TPM": "tpm", "NumReads": "num_reads"})
    tx.to_csv(a.out_tx, sep="\t", index=False)
    g = (tx.groupby(["gene_id", "gene_name"], as_index=False)
           .agg(tpm=("tpm", "sum"), num_reads=("num_reads", "sum")))
    g.to_csv(a.out_gene, sep="\t", index=False)
    print(f"transcripts={len(tx)} genes={len(g)} total_TPM={tx['tpm'].sum():.1f}")


if __name__ == "__main__":
    main()
