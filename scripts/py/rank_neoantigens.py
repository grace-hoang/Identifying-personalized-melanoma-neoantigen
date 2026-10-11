#!/usr/bin/env python3
"""Rank pVACseq class-I candidates with a transparent composite score and pick units for design.

Input : pVACseq <sample>.all_epitopes.tsv   (column names are matched by regex so that minor
        pVACtools version differences do not break the script; missing columns are tolerated)
Output: <out>/<P>.ranked_epitopes.tsv   every epitope row with score, tier, fail reasons
        <out>/<P>.ranked_variants.tsv   best epitope per variant
        <out>/<P>.design_input.tsv      selected units for stage 12
        <out>/<P>.rank_summary.json

Score components (each 0..1; the composite is the weight-normalised mean of the AVAILABLE ones):
  binding        0.30  IC50 (10 nM -> 1, 5000 nM -> 0) averaged with percentile (0 -> 1, 2 -> 0)
  agretopicity   0.15  WT/MT IC50 fold change (1x -> 0, >=16x -> 1); frameshift/novel sequence -> 1
  clonality      0.15  tumour DNA VAF (0.5, or purity/2 when --purity is given, -> 1)
  expression     0.20  log10(1+TPM)/2   (100 TPM -> 1), transcript TPM preferred over gene TPM
  rna_support    0.10  tumour RNA VAF (>=0.25 -> 1)
  hla_breadth    0.05  number of HLA alleles predicted to present the same mutant peptide (>=2 -> 1)
  variant_class  0.05  frameshift 1.0, in-frame indel 0.8, missense 0.6
Tier A: binds AND passes every quality filter that can be evaluated; B: binds but fails >=1 filter; C: no binding.
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

WEIGHTS = {"binding": 0.30, "agretopicity": 0.15, "clonality": 0.15, "expression": 0.20,
           "rna_support": 0.10, "hla_breadth": 0.05, "variant_class": 0.05}


def find_col(df, patterns):
    for pat in patterns:
        rx = re.compile(pat, re.I)
        for c in df.columns:
            if rx.search(c):
                return c
    return None


def num(df, col):
    if col is None:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(df[col], errors="coerce")


def clip01(x):
    return np.clip(x, 0.0, 1.0)


def weighted_mean(components: pd.DataFrame) -> pd.Series:
    w = pd.Series(WEIGHTS)[components.columns]
    avail = components.notna()
    num_ = (components.fillna(0) * w).sum(axis=1)
    den = (avail * w).sum(axis=1)
    return (num_ / den.replace(0, np.nan)).fillna(0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-epitopes", required=True)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ic50-max", type=float, default=500)
    ap.add_argument("--percentile-max", type=float, default=2.0)
    ap.add_argument("--min-fold-change", type=float, default=1.0)
    ap.add_argument("--min-tumor-vaf", type=float, default=0.05)
    ap.add_argument("--min-tumor-dp", type=float, default=10)
    ap.add_argument("--max-normal-vaf", type=float, default=0.02)
    ap.add_argument("--min-expr-tpm", type=float, default=1.0)
    ap.add_argument("--top-n", type=int, default=30)
    ap.add_argument("--select-n", type=int, default=20)
    ap.add_argument("--purity", type=float, default=None)
    a = ap.parse_args()

    df = pd.read_csv(a.all_epitopes, sep="\t", dtype=str).replace({"NA": np.nan, "": np.nan, "None": np.nan})
    if df.empty:
        sys.exit("all_epitopes.tsv is empty — pVACseq found no candidates")

    col = {
        "chrom": find_col(df, [r"^Chromosome$"]), "start": find_col(df, [r"^Start$"]),
        "ref": find_col(df, [r"^Reference$"]), "alt": find_col(df, [r"^Variant$"]),
        "gene": find_col(df, [r"^Gene Name$", r"^Gene$"]), "tx": find_col(df, [r"^Transcript$"]),
        "vtype": find_col(df, [r"^Variant Type$"]), "mutation": find_col(df, [r"^Mutation$"]),
        "ppos": find_col(df, [r"^Protein Position$"]), "hgvsp": find_col(df, [r"^HGVSp$"]),
        "hla": find_col(df, [r"^HLA Allele$"]), "mt": find_col(df, [r"^MT Epitope Seq$"]),
        "wt": find_col(df, [r"^WT Epitope Seq$"]), "mutpos": find_col(df, [r"^Mutation Position$"]),
        "ic50_mt": find_col(df, [r"^Median MT IC50 Score$", r"^Median MT Score$", r"^Best MT IC50 Score$", r"^Best MT Score$"]),
        "ic50_wt": find_col(df, [r"^Median WT IC50 Score$", r"^Median WT Score$",
                                 r"^Corresponding WT IC50 Score$", r"^Corresponding WT Score$"]),
        "fold": find_col(df, [r"^Median Fold Change$", r"^Corresponding Fold Change$", r"Fold Change"]),
        "dna_vaf": find_col(df, [r"^Tumor DNA VAF$"]), "dna_dp": find_col(df, [r"^Tumor DNA Depth$"]),
        "rna_vaf": find_col(df, [r"^Tumor RNA VAF$"]), "norm_vaf": find_col(df, [r"^Normal VAF$"]),
        "gx": find_col(df, [r"^Gene Expression$"]), "tx_expr": find_col(df, [r"^Transcript Expression$"]),
    }
    for need in ("hla", "mt"):
        if col[need] is None:
            sys.exit(f"cannot find required column '{need}'. Columns present: {list(df.columns)}")
    pct_cols = [c for c in df.columns if re.search(r"MT.*Percentile", c, re.I) and not re.search(r"WT", c)]
    pref = [c for c in pct_cols if c.lower().startswith("median")] or pct_cols
    print("matched columns:", {k: v for k, v in col.items() if v}, "| MT percentile columns:", pref, file=sys.stderr)

    ic50 = num(df, col["ic50_mt"])
    wt50 = num(df, col["ic50_wt"])
    pct = pd.concat([num(df, c) for c in pref], axis=1).min(axis=1) if pref else pd.Series(np.nan, index=df.index)
    fold = num(df, col["fold"])
    if fold.isna().all() and not wt50.isna().all():
        fold = wt50 / ic50
    vaf, dp, rvaf, nvaf = num(df, col["dna_vaf"]), num(df, col["dna_dp"]), num(df, col["rna_vaf"]), num(df, col["norm_vaf"])
    tpm = num(df, col["tx_expr"])
    tpm = tpm.where(tpm.notna(), num(df, col["gx"]))

    vtype = df[col["vtype"]].fillna("").str.lower() if col["vtype"] else pd.Series("", index=df.index)
    is_fs = vtype.str.contains("frameshift|fs")
    is_indel = vtype.str.contains("inframe|ins|del") & ~is_fs

    # ---- components ----
    s_ic = clip01((math.log10(5000) - np.log10(ic50.clip(lower=1))) / (math.log10(5000) - math.log10(10)))
    s_pc = clip01(1 - pct / a.percentile_max)
    binding = pd.concat([s_ic, s_pc], axis=1).mean(axis=1, skipna=True)
    agre = clip01(np.log2(fold.clip(lower=1e-3)) / 4.0)
    agre = agre.where(~is_fs, 1.0)
    vaf_target = 0.5 if not a.purity else max(0.05, a.purity / 2.0)
    clonality = clip01(vaf / vaf_target)
    expression = clip01(np.log10(1 + tpm) / 2.0)
    rna_sup = clip01(rvaf / 0.25)
    key_variant = [c for c in (col["chrom"], col["start"], col["ref"], col["alt"], col["tx"]) if c]
    binder = (ic50 <= a.ic50_max) | ((ic50.isna()) & (pct <= a.percentile_max))
    breadth = (df.assign(_b=binder).groupby(df[col["mt"]])["_b"].transform("sum")
               if col["mt"] else pd.Series(1, index=df.index))
    hla_breadth = clip01(breadth.clip(upper=2) / 2.0)
    vclass = pd.Series(0.6, index=df.index).where(~is_indel, 0.8).where(~is_fs, 1.0)

    comps = pd.DataFrame({"binding": binding, "agretopicity": agre, "clonality": clonality, "expression": expression,
                          "rna_support": rna_sup, "hla_breadth": hla_breadth, "variant_class": vclass})
    score = weighted_mean(comps)

    # ---- tiers ----
    bind_ok = binder & ((pct <= a.percentile_max) | pct.isna())
    reasons = pd.Series("", index=df.index)

    def flag(mask, text):
        nonlocal reasons
        reasons = reasons.where(~mask, reasons + text + ";")

    flag(~bind_ok, "weak_binding")
    flag((fold < a.min_fold_change) & ~is_fs, "fold_change<min")
    flag(vaf < a.min_tumor_vaf, "low_tumor_vaf")
    flag(dp < a.min_tumor_dp, "low_tumor_depth")
    flag(nvaf > a.max_normal_vaf, "normal_vaf_high")
    flag(tpm < a.min_expr_tpm, "low_expression")
    tier = np.where(~bind_ok, "C", np.where(reasons == "", "A", "B"))

    out = pd.DataFrame({
        "patient": a.patient,
        "gene": df[col["gene"]] if col["gene"] else "",
        "chrom": df[col["chrom"]] if col["chrom"] else "", "start": df[col["start"]] if col["start"] else "",
        "ref": df[col["ref"]] if col["ref"] else "", "alt": df[col["alt"]] if col["alt"] else "",
        "transcript": df[col["tx"]] if col["tx"] else "", "variant_type": df[col["vtype"]] if col["vtype"] else "",
        "mutation": df[col["mutation"]] if col["mutation"] else "", "protein_position": df[col["ppos"]] if col["ppos"] else "",
        "hgvsp": df[col["hgvsp"]] if col["hgvsp"] else "",
        "hla_allele": df[col["hla"]], "mt_peptide": df[col["mt"]], "wt_peptide": df[col["wt"]] if col["wt"] else "",
        "mutation_position": df[col["mutpos"]] if col["mutpos"] else "",
        "ic50_mt_nM": ic50, "percentile_mt": pct, "fold_change": fold,
        "tumor_dna_vaf": vaf, "tumor_dna_depth": dp, "tumor_rna_vaf": rvaf, "normal_vaf": nvaf, "expression_tpm": tpm,
        "score": score.round(4), "tier": tier, "fail_reasons": reasons.str.rstrip(";"),
    })
    for k in comps.columns:
        out["s_" + k] = comps[k].round(3)
    out = out.sort_values(["tier", "score"], ascending=[True, False]).reset_index(drop=True)
    out["rank_epitope"] = np.arange(1, len(out) + 1)

    os.makedirs(a.out_dir, exist_ok=True)
    out.to_csv(f"{a.out_dir}/{a.patient}.ranked_epitopes.tsv", sep="\t", index=False)

    # ---- variant level ----
    vkeys = ["chrom", "start", "ref", "alt", "transcript"]
    out["variant_id"] = out[vkeys].astype(str).agg("|".join, axis=1)
    best = out.sort_values(["score"], ascending=False).groupby("variant_id", as_index=False).head(1)
    nb = out[out.tier != "C"].groupby("variant_id").size().rename("n_binding_epitopes")
    best = best.merge(nb, on="variant_id", how="left").fillna({"n_binding_epitopes": 0})
    tier_order = {"A": 0, "B": 1, "C": 2}
    best = best.assign(_t=best.tier.map(tier_order)).sort_values(["_t", "score"], ascending=[True, False]).drop(columns="_t")
    best["rank_variant"] = np.arange(1, len(best) + 1)
    best.to_csv(f"{a.out_dir}/{a.patient}.ranked_variants.tsv", sep="\t", index=False)
    best.head(a.top_n).to_csv(f"{a.out_dir}/{a.patient}.top_candidates.tsv", sep="\t", index=False)

    # ---- selection for design: tier A first, then B, no duplicate / nested peptides ----
    chosen, seqs = [], []
    for _, r in best[best.tier != "C"].iterrows():
        p = str(r["mt_peptide"])
        if not p or p == "nan" or re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", p):
            continue
        if any(p in s or s in p for s in seqs):
            continue
        chosen.append(r)
        seqs.append(p)
        if len(chosen) >= a.select_n:
            break
    sel = pd.DataFrame(chosen)
    if sel.empty:
        sel = pd.DataFrame(columns=best.columns)
    sel = sel.assign(design_rank=np.arange(1, len(sel) + 1))
    sel.to_csv(f"{a.out_dir}/{a.patient}.design_input.tsv", sep="\t", index=False)

    summary = {
        "patient": a.patient, "epitope_rows": int(len(out)), "variants": int(len(best)),
        "tier_counts_epitopes": out.tier.value_counts().to_dict(),
        "tier_counts_variants": best.tier.value_counts().to_dict(),
        "selected_for_design": int(len(sel)),
        "columns_used": {k: v for k, v in col.items() if v},
        "columns_missing": [k for k, v in col.items() if not v],
    }
    with open(f"{a.out_dir}/{a.patient}.rank_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    if sel.empty:
        print("WARNING: no candidate passed binding thresholds — design stage will have nothing to build", file=sys.stderr)


if __name__ == "__main__":
    main()
