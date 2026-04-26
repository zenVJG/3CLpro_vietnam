#!/usr/bin/env python3

# ============================================================
# Filter BLAST results by ambiguous characters in btop column
# Usage:
#   python filter_sequences.py
#   python filter_sequences.py --nuc my_blast.txt --prot my_protein.txt
# ============================================================

import argparse

# Ký tự cần loại trong nucleotide (ngoài ATGC)
AMBIGUOUS_NUC = set("NYRMKWSBDHVnyrmkwsbdhv")

# Ký tự cần loại trong protein
AMBIGUOUS_PROT = {"X", "x"}


def filter_blast(input_file, output_file):
    kept, removed = 0, 0
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            cols = line.strip().split("\t")
            btop = cols[4] if len(cols) >= 5 else ""

            # Loại dòng nếu có bất kỳ ký tự ambiguous nào trong btop
            if any(char in AMBIGUOUS_NUC for char in btop):
                removed += 1
            else:
                fout.write(line)
                kept += 1

    print(f"[nucleotide] Input : {input_file}")
    print(f"[nucleotide] Output: {output_file}")
    print(f"[nucleotide] Giữ lại: {kept} | Loại bỏ: {removed} (có N/Y/R/M/K...)")
    print()


def filter_protein(input_file, output_file):
    kept, removed = 0, 0
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            cols = line.strip().split("\t")
            btop = cols[4] if len(cols) >= 5 else ""

            # Loại dòng nếu có X trong btop
            if any(char in AMBIGUOUS_PROT for char in btop):
                removed += 1
            else:
                fout.write(line)
                kept += 1

    print(f"[protein]    Input : {input_file}")
    print(f"[protein]    Output: {output_file}")
    print(f"[protein]    Giữ lại: {kept} | Loại bỏ: {removed} (có X)")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BLAST results by ambiguous characters")
    parser.add_argument("--nuc",      default="blast_detail.txt",            help="Input nucleotide blast file")
    parser.add_argument("--prot",     default="protein_detail.txt",          help="Input protein blast file")
    parser.add_argument("--nuc_out",  default="blast_detail_filtered.txt",   help="Output nucleotide file")
    parser.add_argument("--prot_out", default="protein_detail_filtered.txt", help="Output protein file")
    args = parser.parse_args()

    filter_blast(args.nuc, args.nuc_out)
    filter_protein(args.prot, args.prot_out)