#!/usr/bin/env python3
"""Generate a small synthetic pVACseq all_epitopes.tsv + matching VEP-style VCF to smoke-test
stages 11 and 12 without running the whole pipeline.   python3 tests/make_synthetic.py OUTDIR"""
import gzip, os, random, sys

out = sys.argv[1] if len(sys.argv) > 1 else "/tmp/melsyn"
os.makedirs(out, exist_ok=True)
random.seed(1)
AA = "ACDEFGHIKLMNPQRSTVWY"
hla = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
genes = [f"GENE{i}" for i in range(1, 31)]
hdr = ["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "Ensembl Gene ID", "Variant Type", "Mutation",
       "Protein Position", "Gene Name", "HGVSc", "HGVSp", "HLA Allele", "Peptide Length", "Sub-peptide Position",
       "Mutation Position", "MT Epitope Seq", "WT Epitope Seq", "Median MT IC50 Score", "Median MT Percentile",
       "Median WT IC50 Score", "Median Fold Change", "Tumor DNA Depth", "Tumor DNA VAF", "Tumor RNA Depth",
       "Tumor RNA VAF", "Normal Depth", "Normal VAF", "Gene Expression", "Transcript Expression"]
rows, vcf = [], []
for i, g in enumerate(genes):
    prot = "M" + "".join(random.choice(AA) for _ in range(120))
    pos = random.randint(30, 90)               # 1-based protein position
    wt_aa = prot[pos - 1]
    mt_aa = random.choice([x for x in AA if x != wt_aa])
    mprot = prot[: pos - 1] + mt_aa + prot[pos:]
    tx = f"ENST{100000 + i:011d}"
    chrom, gpos = f"chr{1 + i % 20}", 1000000 + i * 5000
    for h in hla:
        L = random.choice([8, 9, 10, 11]); off = random.randint(0, L - 1)
        s = pos - 1 - off
        mt, wt = mprot[s:s + L], prot[s:s + L]
        ic = random.choice([8, 30, 120, 400, 900, 3000]) * random.uniform(0.7, 1.3)
        wic = ic * random.choice([0.5, 1, 3, 10, 40])
        rows.append([chrom, gpos, gpos, "A", "T", tx, f"ENSG{i:011d}", "missense_variant", f"{wt_aa}{pos}{mt_aa}",
                     str(pos), g, "c.1A>T", f"p.{wt_aa}{pos}{mt_aa}", h, L, s + 1, off + 1, mt, wt, f"{ic:.1f}",
                     f"{random.uniform(0.05, 4):.2f}", f"{wic:.1f}", f"{wic / ic:.2f}", random.randint(15, 120),
                     f"{random.uniform(0.04, 0.55):.2f}", random.randint(5, 80), f"{random.uniform(0, 0.6):.2f}", 60,
                     "0.00", f"{random.choice([0, 0.5, 3, 20, 90]):.2f}", f"{random.choice([0, 0.4, 2, 15, 80]):.2f}"])
    csq = f"T|missense_variant|MODERATE|{g}|ENSG{i:011d}|Transcript|{tx}|protein_coding|{pos}|{wt_aa}/{mt_aa}|{prot}"
    vcf.append(f"{chrom}\t{gpos}\t.\tA\tT\t.\tPASS\tCSQ={csq}\tGT\t0/1")
# one frameshift row
rows.append(["chr5", 777, 777, "-", "A", "ENST00000999999", "ENSG00000999999", "frameshift_variant", "FS", "10", "FSGENE", "", "",
             "HLA-A*02:01", 9, 1, 1, "KLVQWERTY", "NA", "35.0", "0.4", "NA", "NA", 50, "0.40", 30, "0.30", 40, "0.00", "12.0", "10.0"])
with open(f"{out}/syn.all_epitopes.tsv", "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
with gzip.open(f"{out}/syn.vep.vcf.gz", "wt") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    fh.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|Protein_position|Amino_acids|WildtypeProtein">\n')
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
    fh.write("\n".join(vcf) + "\n")
print("wrote", out)
