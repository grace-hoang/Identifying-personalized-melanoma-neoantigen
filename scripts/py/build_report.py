#!/usr/bin/env python3
"""Self-contained per-patient HTML report. Every section degrades gracefully if its inputs are missing."""
import argparse
import datetime
import glob
import gzip
import html
import json
import os
import re

import pandas as pd

CSS = """body{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:2rem auto;max-width:1150px;color:#1d2330;line-height:1.45}
h1{margin-bottom:.2rem}h2{border-bottom:2px solid #dfe3ea;padding-bottom:.2rem;margin-top:2rem}
table{border-collapse:collapse;margin:.6rem 0;font-size:.86rem}th,td{border:1px solid #d5dae3;padding:.25rem .55rem;text-align:left}
th{background:#f1f4f9}.warn{background:#fff6e0;border:1px solid #f0d18a;padding:.6rem .9rem;border-radius:6px}
.ok{color:#17803d}.bad{color:#b3261e}code,pre{background:#f4f6fa;padding:.1rem .3rem;border-radius:4px;font-size:.8rem}
pre{padding:.6rem;overflow-x:auto;white-space:pre-wrap;word-break:break-all}.small{color:#5b6475;font-size:.8rem}"""


def esc(x):
    return html.escape("" if x is None else str(x))


def table(df, max_rows=None, floatfmt="{:.3g}"):
    if df is None or len(df) == 0:
        return "<p class='small'>no data</p>"
    d = df.head(max_rows) if max_rows else df
    head = "".join(f"<th>{esc(c)}</th>" for c in d.columns)
    body = ""
    for _, r in d.iterrows():
        cells = []
        for v in r:
            cells.append(f"<td>{esc(floatfmt.format(v) if isinstance(v, float) and v == v else ('' if v != v else v))}</td>")
        body += "<tr>" + "".join(cells) + "</tr>"
    return f"<table><tr>{head}</tr>{body}</table>"


def section(title, content):
    return f"<h2>{esc(title)}</h2>\n{content}\n"


def safe(fn, *a, **k):
    try:
        return fn(*a, **k)
    except Exception as e:  # noqa
        return f"<p class='small'>section unavailable: {esc(e)}</p>"


def fastp_summary(res, patient):
    rows = []
    for j in sorted(glob.glob(f"{res}/qc/fastp/{patient}_*.fastp.json")):
        d = json.load(open(j))
        b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
        rows.append({"library": os.path.basename(j).replace(".fastp.json", ""), "reads_before": b["total_reads"],
                     "reads_after": a["total_reads"], "Q30_after_%": round(100 * a["q30_rate"], 1),
                     "GC_after_%": round(100 * a["gc_content"], 1),
                     "dup_rate_%": round(100 * d.get("duplication", {}).get("rate", 0), 1)})
    return table(pd.DataFrame(rows), floatfmt="{:,.1f}")


def flagstat(res, patient):
    rows = []
    for s in ("Tumor", "Normal"):
        f = f"{res}/alignment/{patient}/{patient}_{s}.flagstat.txt"
        if not os.path.exists(f):
            continue
        t = open(f).read()
        tot = int(re.search(r"^(\d+) \+ \d+ in total", t, re.M).group(1))
        mp = re.search(r"^(\d+) \+ \d+ mapped \(([\d.]+)%", t, re.M)
        dup = re.search(r"^(\d+) \+ \d+ duplicates", t, re.M)
        row = {"sample": s, "total_records": tot, "mapped_%": float(mp.group(2)),
               "duplicate_%": round(100 * int(dup.group(1)) / tot, 1) if dup else None}
        ms = f"{res}/qc/mosdepth/{patient}_{s}.mosdepth.summary.txt"
        if os.path.exists(ms):
            m = pd.read_csv(ms, sep="\t")
            r = m[m["chrom"] == ("total_region" if "total_region" in set(m["chrom"]) else "total")].iloc[0]
            row["mean_depth_" + ("targets" if "total_region" in set(m["chrom"]) else "genome-wide")] = float(r["mean"])
        rows.append(row)
    return table(pd.DataFrame(rows), floatfmt="{:,.2f}")


def hla(res, patient):
    f = f"{res}/hla/{patient}/{patient}.hla_alleles.tsv"
    return table(pd.read_csv(f, sep="\t")) if os.path.exists(f) else "<p class='small'>no HLA result</p>"


def variants(res, patient):
    f = f"{res}/variants/{patient}/{patient}.somatic.pass.norm.vep.vcf.gz"
    if not os.path.exists(f):
        return "<p class='small'>no annotated VCF</p>"
    fmt, cons, vt = None, {}, {"SNV": 0, "indel": 0}
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format:")[1].strip().rstrip('">').split("|")
            if line.startswith("#"):
                continue
            c = line.split("\t")
            vt["SNV" if len(c[3]) == len(c[4]) else "indel"] += 1
            m = re.search(r"CSQ=([^;\t]+)", c[7])
            if m and fmt:
                d = dict(zip(fmt, m.group(1).split(",")[0].split("|")))
                for t in d.get("Consequence", "").split("&")[:1]:
                    cons[t] = cons.get(t, 0) + 1
    t1 = pd.DataFrame([{"class": k, "count": v} for k, v in vt.items()])
    t2 = pd.DataFrame(sorted(({"consequence": k, "count": v} for k, v in cons.items()), key=lambda x: -x["count"]))
    return f"<p>PASS somatic variants: <b>{sum(vt.values())}</b></p>" + table(t1) + table(t2, 15)


def neo(res, patient, algorithms):
    d = f"{res}/neoantigen/{patient}"
    out = ""
    sj = f"{d}/{patient}.rank_summary.json"
    if os.path.exists(sj):
        s = json.load(open(sj))
        out += (f"<p>{s['epitope_rows']} epitope–HLA rows, {s['variants']} variants. Tier counts (variants): "
                f"{esc(s['tier_counts_variants'])}. Selected for design: <b>{s['selected_for_design']}</b>.</p>"
                f"<p class='small'>Prediction algorithms: {esc(algorithms)}. Tier A = binds and passes all evaluable quality filters; "
                f"B = binds but fails a filter; C = weak binding.</p>")
    f = f"{d}/{patient}.top_candidates.tsv"
    if os.path.exists(f):
        df = pd.read_csv(f, sep="\t")
        cols = [c for c in ["rank_variant", "gene", "mutation", "hla_allele", "mt_peptide", "wt_peptide", "ic50_mt_nM",
                            "percentile_mt", "fold_change", "tumor_dna_vaf", "tumor_rna_vaf", "expression_tpm", "score", "tier",
                            "fail_reasons"] if c in df.columns]
        out += table(df[cols], 25)
    return out or "<p class='small'>no neoantigen results</p>"


def construct(res, patient):
    d = f"{res}/mrna/{patient}"
    rj = f"{d}/{patient}.design_report.json"
    if not os.path.exists(rj):
        return "<p class='small'>no construct designed</p>"
    r = json.load(open(rj))
    keys = ["n_units", "protein_length_aa", "cds_length_nt", "mrna_length_nt_incl_polyA", "gc_cds", "gc_mrna", "cai",
            "mfe_5utr_plus_60nt_kcal_mol", "mfe_full_mrna_kcal_mol", "internal_stop_codons", "translation_verified"]
    out = table(pd.DataFrame([{"metric": k, "value": r.get(k)} for k in keys]))
    if r.get("warnings"):
        out += "<div class='warn'><b>Warnings</b><ul>" + "".join(f"<li>{esc(w)}</li>" for w in r["warnings"]) + "</ul></div>"
    mp = f"{d}/{patient}.construct.map.tsv"
    if os.path.exists(mp):
        m = pd.read_csv(mp, sep="\t").drop(columns=["sequence"], errors="ignore")
        out += table(m)
    pf = f"{d}/{patient}.construct.protein.fa"
    if os.path.exists(pf):
        out += "<p>Protein:</p><pre>" + esc("".join(open(pf).read().split("\n")[1:])) + "</pre>"
    out += "<div class='warn'><b>Caveats</b><ul>" + "".join(f"<li>{esc(c)}</li>" for c in r.get("caveats", [])) + "</ul></div>"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--multiqc", default="")
    ap.add_argument("--aligner", default="")
    ap.add_argument("--gencode", default="")
    ap.add_argument("--ensembl", default="")
    ap.add_argument("--algorithms", default="")
    a = ap.parse_args()
    res, p = a.results_dir, a.patient
    h = [f"<html><head><meta charset='utf-8'><title>{esc(p)} neoantigen report</title><style>{CSS}</style></head><body>",
         f"<h1>Melanoma neoantigen pipeline — {esc(p)}</h1>",
         f"<p class='small'>Generated {datetime.datetime.now():%Y-%m-%d %H:%M} · aligner {esc(a.aligner)} · GRCh38 · "
         f"GENCODE {esc(a.gencode)} / Ensembl {esc(a.ensembl)}</p>",
         "<div class='warn'><b>Research use only.</b> Computational predictions; not validated and not for clinical use.</div>"]
    h.append(section("Read QC and trimming (fastp)", safe(fastp_summary, res, p)))
    h.append(section("Alignment and coverage", safe(flagstat, res, p)))
    h.append(section("HLA class I genotype", safe(hla, res, p)))
    h.append(section("Somatic variants", safe(variants, res, p)))
    h.append(section("Ranked neoantigen candidates", safe(neo, res, p, a.algorithms)))
    h.append(section("Candidate mRNA construct", safe(construct, res, p)))
    if a.multiqc and os.path.exists(a.multiqc):
        h.append(section("Full QC", f"<p><a href='{esc(os.path.relpath(a.multiqc, os.path.dirname(a.out)))}'>MultiQC report</a></p>"))
    h.append("</body></html>")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write("\n".join(h))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
