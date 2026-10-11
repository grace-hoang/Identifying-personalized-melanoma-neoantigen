#!/usr/bin/env python3
"""In-silico design of a candidate multi-epitope neoantigen mRNA construct. RESEARCH USE ONLY.

Protein   : [signal peptide] - unit1 - linker - unit2 - ... - [PADRE] - [MITD]
            units are mutant 25-mers (missense, built from VEP's WildtypeProtein) or bare MT epitopes.
CDS       : codon-optimised by simulated annealing (human codon usage, GC window, forbidden sites,
            homopolymers, repeats), translation verified.
mRNA      : 5'UTR + CDS + stop + 3'UTR + poly(A); DNA template adds a T7 promoter.
QC        : GC, CAI, motif/homopolymer scan, ViennaRNA MFE (if RNAfold is installed).
Outputs   : <P>.construct.{protein,cds,mrna,template}.fa, <P>.construct.map.tsv, <P>.design_report.json
"""
import argparse
import csv
import gzip
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
from collections import defaultdict

CODONS = {
    "F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "I": ["ATT", "ATC", "ATA"], "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"], "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "A": ["GCT", "GCC", "GCA", "GCG"], "Y": ["TAT", "TAC"], "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"], "N": ["AAT", "AAC"], "K": ["AAA", "AAG"], "D": ["GAT", "GAC"], "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"], "W": ["TGG"], "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "G": ["GGT", "GGC", "GGA", "GGG"],
}
# Homo sapiens codon usage, per thousand (Kazusa-style table; approximate).
USAGE = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0, "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4, "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1, "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 12.2, "TAC": 15.3, "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2, "AAT": 17.0, "AAC": 19.1,
    "AAA": 24.4, "AAG": 31.9, "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6, "TGT": 10.6, "TGC": 12.6,
    "TGG": 13.2, "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4, "AGT": 12.1, "AGC": 19.5, "AGA": 12.2,
    "AGG": 12.0, "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}
CODON2AA = {c: aa for aa, cs in CODONS.items() for c in cs}
W = {}                                    # relative adaptiveness
for aa, cs in CODONS.items():
    m = max(USAGE[c] for c in cs)
    for c in cs:
        W[c] = USAGE[c] / m
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def gc(s):
    return (s.count("G") + s.count("C")) / max(1, len(s))


def translate(cds):
    return "".join(CODON2AA.get(cds[i:i + 3], "*") for i in range(0, len(cds) - len(cds) % 3, 3))


def cai(cds):
    ws = [W[cds[i:i + 3]] for i in range(0, len(cds), 3) if cds[i:i + 3] in W]
    return math.exp(sum(math.log(w) for w in ws) / len(ws)) if ws else 0.0


def read_seq(path, required=True):
    if not os.path.exists(path):
        if required:
            sys.exit(f"missing design element file: {path}")
        return ""
    lines = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    return "".join(lines).upper().replace(" ", "")


def to_dna(s):
    return s.upper().replace("U", "T")


# ----------------------------------------------------------------------------- flanks
def load_vep_proteins(vcf):
    """(chrom, pos, transcript_no_version) -> (protein_position, 'R/H', WildtypeProtein)"""
    out = {}
    if not vcf or not os.path.exists(vcf):
        return out
    op = gzip.open if vcf.endswith(".gz") else open
    fmt = None
    with op(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format:")[1].strip().rstrip('">').split("|")
                continue
            if line.startswith("#") or fmt is None:
                continue
            f = line.rstrip("\n").split("\t")
            m = re.search(r"(?:^|;)CSQ=([^;]+)", f[7])
            if not m:
                continue
            for ent in m.group(1).split(","):
                d = dict(zip(fmt, ent.split("|")))
                tx = d.get("Feature", "").split(".")[0]
                if tx and d.get("WildtypeProtein"):
                    out[(f[0], int(f[1]), tx)] = (d.get("Protein_position", ""), d.get("Amino_acids", ""), d["WildtypeProtein"])
    return out


def build_unit(row, prots, flank):
    """Return (sequence, kind, note). Falls back to the bare MT epitope."""
    mt = row["mt_peptide"]
    base = (mt, "epitope", "bare MT epitope")
    if flank <= 0 or "missense" not in str(row.get("variant_type", "")).lower():
        return base
    try:
        key = (row["chrom"], int(float(row["start"])), str(row["transcript"]).split(".")[0])
        pp, aa, wt = prots[key]
        if "-" in pp or "/" not in aa:
            return base
        ref_aa, alt_aa = aa.split("/")
        if len(ref_aa) != 1 or len(alt_aa) != 1:
            return base
        p0 = int(pp) - 1
        if wt[p0] != ref_aa:
            return base
        s, e = max(0, p0 - flank), min(len(wt), p0 + flank + 1)
        win = wt[s:p0] + alt_aa + wt[p0 + 1:e]
        if mt not in win:
            return base
        return (win, "minigene", f"{len(win)}-mer around p.{ref_aa}{pp}{alt_aa}")
    except (KeyError, ValueError, IndexError):
        return base


# ----------------------------------------------------------------------------- ordering
def max_window_gravy(seq, w=9):
    if len(seq) < w:
        return sum(KD.get(a, 0) for a in seq) / max(1, len(seq))
    vals = [KD.get(a, 0) for a in seq[:w]]
    cur = sum(vals); best = cur
    for i in range(w, len(seq)):
        cur += KD.get(seq[i], 0) - KD.get(seq[i - w], 0)
        best = max(best, cur)
    return best / w


def order_units(units, linker_for, rng, iters=3000):
    """Keep high-ranked units early but avoid hydrophobic stretches at junctions (aggregation / TM-like)."""
    def cost(order):
        s = ""
        for i, u in enumerate(order):
            s += (linker_for(u) if i else "") + u["seq"]
        pos_pen = sum(i * u["rank"] for i, u in enumerate(order)) / max(1, len(order)) * 0.002
        return max_window_gravy(s) + pos_pen
    order = sorted(units, key=lambda u: u["rank"])
    best, best_c = order[:], cost(order)
    for _ in range(iters):
        if len(order) < 2:
            break
        i, j = rng.sample(range(len(order)), 2)
        cand = order[:]
        cand[i], cand[j] = cand[j], cand[i]
        c = cost(cand)
        if c < best_c:
            best, best_c, order = cand, c, cand
    return best


# ----------------------------------------------------------------------------- codons
def penalty(cds, gc_min, gc_max, sites):
    p = 0.0
    g = gc(cds)
    if g < gc_min:
        p += 80 * (gc_min - g)
    if g > gc_max:
        p += 80 * (g - gc_max)
    for s in sites:
        p += 10 * (cds.count(s) + cds.count(revcomp(s)))
    p += 6 * len(re.findall(r"A{6,}|T{6,}|G{6,}|C{7,}", cds))
    p += 8 * len(re.findall(r"AATAAA|ATTAAA", cds))
    bad = 0
    for i in range(0, max(0, len(cds) - 40), 20):
        wg = gc(cds[i:i + 40])
        if wg < 0.30 or wg > 0.80:
            bad += 1
    p += 0.5 * bad
    ks = [cds[i:i + 14] for i in range(0, len(cds) - 13)]
    p += 1.0 * (len(ks) - len(set(ks))) / 10
    return p


def optimise(protein, rng, iters, gc_min, gc_max, sites):
    aa_list = list(protein)
    cod = []
    for aa in aa_list:
        cs = CODONS[aa]
        cod.append(rng.choices(cs, weights=[USAGE[c] ** 2 for c in cs])[0])
    cur = "".join(cod)

    def J(seq):
        return penalty(seq, gc_min, gc_max, sites) + 25 * (1 - cai(seq))

    cur_j = J(cur)
    best, best_j = cur, cur_j
    mutable = [i for i, aa in enumerate(aa_list) if len(CODONS[aa]) > 1]
    for it in range(iters):
        T = 2.0 * (1 - it / iters) + 0.01
        i = rng.choice(mutable)
        old = cod[i]
        new = rng.choice([c for c in CODONS[aa_list[i]] if c != old])
        cod[i] = new
        cand = "".join(cod)
        j = J(cand)
        if j <= cur_j or rng.random() < math.exp((cur_j - j) / T):
            cur, cur_j = cand, j
            if j < best_j:
                best, best_j = cand, j
        else:
            cod[i] = old
    return best


# ----------------------------------------------------------------------------- RNAfold
def rnafold(seq):
    exe = shutil.which("RNAfold")
    if not exe:
        return None
    try:
        r = subprocess.run([exe, "--noPS"], input=seq.replace("T", "U") + "\n", capture_output=True, text=True, timeout=600)
        m = re.search(r"\(\s*(-?\d+\.\d+)\)\s*$", r.stdout.strip().splitlines()[-1])
        return float(m.group(1)) if m else None
    except Exception:
        return None


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", required=True)
    ap.add_argument("--design-input", required=True)
    ap.add_argument("--vep-vcf", default="")
    ap.add_argument("--design-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-epitopes", type=int, default=20)
    ap.add_argument("--flank", type=int, default=12)
    ap.add_argument("--linker-short", default="AAY")
    ap.add_argument("--linker-long", default="GGSGGGGSGG")
    ap.add_argument("--signal", type=int, default=1)
    ap.add_argument("--mitd", type=int, default=1)
    ap.add_argument("--padre", type=int, default=0)
    ap.add_argument("--polya-len", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--iters", type=int, default=6000)
    ap.add_argument("--gc-min", type=float, default=0.48)
    ap.add_argument("--gc-max", type=float, default=0.62)
    ap.add_argument("--avoid-sites", default="GGTCTC GAGACC CGTCTC GAGACG GAATTC AAGCTT")
    ap.add_argument("--t7-promoter", default="TAATACGACTCACTATAG")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(a.out_dir, exist_ok=True)
    warnings = []

    with open(a.design_input, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows = rows[: a.max_epitopes]
    if not rows:
        sys.exit("design input is empty: no neoantigen passed the ranking filters")

    prots = load_vep_proteins(a.vep_vcf) if a.flank > 0 else {}
    if a.flank > 0 and not prots:
        warnings.append("no WildtypeProtein found in VEP VCF: using bare epitopes")
    units = []
    for r in rows:
        seq, kind, note = build_unit(r, prots, a.flank)
        if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", seq):
            continue
        units.append({"seq": seq, "kind": kind, "note": note, "rank": int(float(r.get("design_rank") or 1)),
                      "gene": r.get("gene", ""), "mutation": r.get("mutation", "") or r.get("hgvsp", ""),
                      "hla": r.get("hla_allele", ""), "epitope": r["mt_peptide"], "score": r.get("score", ""),
                      "tier": r.get("tier", "")})
    if not units:
        sys.exit("no usable units after sequence validation")

    linker_for = lambda u: a.linker_long if len(u["seq"]) >= 20 else a.linker_short
    units = order_units(units, linker_for, rng)

    signal = read_seq(f"{a.design_dir}/signal_peptide.txt") if a.signal else ""
    mitd = read_seq(f"{a.design_dir}/mitd.txt") if a.mitd else ""
    padre = read_seq(f"{a.design_dir}/padre.txt") if a.padre else ""
    utr5 = to_dna(read_seq(f"{a.design_dir}/five_prime_utr.txt"))
    utr3 = to_dna(read_seq(f"{a.design_dir}/three_prime_utr.txt"))

    parts = []          # (label, aa sequence, meta)
    if signal:
        parts.append(("signal_peptide", signal, None))
    elif not units[0]["seq"].startswith("M"):
        parts.append(("start_Met", "M", None))
    for i, u in enumerate(units):
        if parts and (i > 0):
            parts.append(("linker", linker_for(u), None))
        elif parts and i == 0:
            parts.append(("linker", "GS", None))
        parts.append((f"unit{i + 1}", u["seq"], u))
    if padre:
        parts.append(("linker", a.linker_long, None)); parts.append(("PADRE", padre, None))
    if mitd:
        parts.append(("linker", a.linker_long, None)); parts.append(("MITD", mitd, None))
    protein = "".join(p[1] for p in parts)
    if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", protein):
        sys.exit("protein contains non-standard residues")

    sites = a.avoid_sites.split()
    cds_body = optimise(protein, rng, a.iters, a.gc_min, a.gc_max, sites)
    cds = cds_body + "TGA"
    if translate(cds_body) != protein:
        sys.exit("internal error: optimised CDS does not translate back to the protein")

    polya = "A" * a.polya_len
    mrna_dna = utr5 + cds + utr3 + polya
    template = a.t7_promoter + mrna_dna

    # map file
    pos = 0
    with open(f"{a.out_dir}/{a.patient}.construct.map.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["element", "aa_start", "aa_end", "gene", "mutation", "hla", "mt_epitope", "unit_type", "score", "tier", "sequence"])
        for lab, seq, meta in parts:
            w.writerow([lab, pos + 1, pos + len(seq)] + ([meta["gene"], meta["mutation"], meta["hla"], meta["epitope"], meta["note"], meta["score"], meta["tier"]] if meta else ["", "", "", "", "", "", ""]) + [seq])
            pos += len(seq)

    for name, seq, hdr in [
        ("protein", protein, f"{a.patient}_neoantigen_polyepitope_protein len={len(protein)}"),
        ("cds", cds, f"{a.patient}_cds len={len(cds)} RESEARCH_USE_ONLY"),
        ("mrna", mrna_dna.replace("T", "U"), f"{a.patient}_candidate_mRNA len={len(mrna_dna)} (cap/modified nucleosides not represented) RESEARCH_USE_ONLY"),
        ("template", template, f"{a.patient}_DNA_template_with_T7 len={len(template)}"),
    ]:
        with open(f"{a.out_dir}/{a.patient}.construct.{name}.fa", "w") as fh:
            fh.write(f">{hdr}\n{wrap(seq)}\n")

    # QC
    full_sites = {s: mrna_dna.count(s) + mrna_dna.count(revcomp(s)) for s in sites}
    cds_hits = {s: v for s, v in full_sites.items() if v}
    internal_stop = "*" in translate(cds_body)
    homo = re.findall(r"A{7,}|T{7,}|G{6,}|C{7,}", cds)
    mfe5 = rnafold((utr5 + cds[:60]))
    mfe_all = rnafold(mrna_dna) if len(mrna_dna) <= 6000 else None
    if mfe5 is None:
        warnings.append("RNAfold unavailable: MFE not computed")
    if gc(cds) < a.gc_min or gc(cds) > a.gc_max:
        warnings.append(f"CDS GC {gc(cds):.2f} outside [{a.gc_min},{a.gc_max}]")
    if cds_hits:
        warnings.append(f"forbidden sites remain in transcript (may be in UTRs): {cds_hits}")
    if "ATG" in utr5[:-3] and "ATG" in utr5:
        warnings.append("upstream ATG in 5'UTR — check uORF risk")
    report = {
        "patient": a.patient, "research_use_only": True,
        "n_units": len(units), "unit_types": {k: sum(1 for u in units if u["kind"] == k) for k in {u["kind"] for u in units}},
        "protein_length_aa": len(protein), "cds_length_nt": len(cds), "mrna_length_nt_incl_polyA": len(mrna_dna),
        "gc_cds": round(gc(cds), 4), "gc_mrna": round(gc(mrna_dna), 4), "cai": round(cai(cds_body), 4),
        "max_window_gravy": round(max_window_gravy(protein), 3), "internal_stop_codons": internal_stop,
        "translation_verified": True, "homopolymers_ge6_7": homo, "forbidden_sites_found": cds_hits,
        "mfe_5utr_plus_60nt_kcal_mol": mfe5, "mfe_full_mrna_kcal_mol": mfe_all,
        "elements": {"signal": bool(signal), "mitd": bool(mitd), "padre": bool(padre), "polyA": a.polya_len},
        "warnings": warnings,
        "caveats": ["Junctional neoepitopes between units were NOT screened; run pVACvector or NetMHCpan on junctions.",
                    "UTR/signal/MITD sequences are configurable examples; validate before any synthesis.",
                    "In-silico candidate only; not validated for expression, immunogenicity or safety."],
    }
    with open(f"{a.out_dir}/{a.patient}.design_report.json", "w") as fh:
        json.dump(report, fh, indent=2)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
