#!/usr/bin/env python3
"""Convert an OptiType result table (or a manual override line) into pVACseq allele names.

  hla_parse.py --patient P --optitype P_result.tsv --out-txt P.hla_alleles.txt --out-tsv P.hla_alleles.tsv
  hla_parse.py --patient P --override hla_override.tsv ...
"""
import argparse
import csv
import re
import sys

ALLELE_RX = re.compile(r"^HLA-[ABC]\*\d{2,3}:\d{2,3}[A-Z]?$")


def norm(a: str) -> str:
    a = a.strip()
    if not a or a.lower() in {"none", "nan", "na"}:
        return ""
    if not a.upper().startswith("HLA-"):
        a = "HLA-" + a
    return "HLA-" + a[4:]


def from_optitype(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"empty OptiType result: {path}")
    r = rows[0]
    out = []
    for k in ("A1", "A2", "B1", "B2", "C1", "C2"):
        v = norm(r.get(k, ""))
        if v:
            out.append(v)
    return out


def from_override(path, patient):
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t") if "\t" in line else line.split(None, 1)
            if parts[0].strip() == patient and len(parts) > 1:
                return [norm(x) for x in parts[1].split(",") if norm(x)]
    sys.exit(f"patient {patient} not found in {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", required=True)
    ap.add_argument("--optitype")
    ap.add_argument("--override")
    ap.add_argument("--out-txt", required=True)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()
    if bool(a.optitype) == bool(a.override):
        sys.exit("give exactly one of --optitype / --override")
    alleles = from_optitype(a.optitype) if a.optitype else from_override(a.override, a.patient)
    bad = [x for x in alleles if not ALLELE_RX.match(x)]
    if bad:
        sys.exit(f"unrecognised allele format: {bad}")
    uniq = list(dict.fromkeys(alleles))      # homozygous loci collapse to one allele
    if len(uniq) < 3:
        print(f"WARNING: only {len(uniq)} distinct class-I alleles for {a.patient}", file=sys.stderr)
    with open(a.out_txt, "w") as fh:
        fh.write(",".join(uniq) + "\n")
    with open(a.out_tsv, "w") as fh:
        fh.write("patient\tallele\tsource\n")
        for x in uniq:
            fh.write(f"{a.patient}\t{x}\t{'optitype' if a.optitype else 'override'}\n")
    print(",".join(uniq))


if __name__ == "__main__":
    main()
