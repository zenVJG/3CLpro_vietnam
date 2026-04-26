#!/usr/bin/env python3

# ============================================================
# Thống kê đột biến từ kết quả BLAST (btop format)
# - blast_detail_filtered.txt : nucleotide (offset genome = 10054)
# - protein_data.txt          : protein (có header dòng đầu)
# ============================================================

import re
import argparse
from collections import defaultdict


# ── Hằng số ────────────────────────────────────────────────
NUC_FILE     = "blast_detail_filtered.txt"
PROT_FILE    = "protein_data.txt"
NUC_OUT      = "mutation_stats_nuc.txt"
PROT_OUT     = "mutation_stats_prot.txt"


# ── Hàm parse btop ─────────────────────────────────────────
def parse_btop(btop: str):
    """
    Phân tích chuỗi btop, trả về list các đột biến:
    [(position_in_query, ref_char, alt_char), ...]

    Ví dụ btop: '143CT165GT82'
      → 143 match → pos 144: C→T → 165 match → pos 310: G→T → 82 match
    """
    mutations = []
    pos = 0  # vị trí hiện tại trên query (0-based)

    # Tách thành tokens: số hoặc cặp chữ
    tokens = re.findall(r'\d+|[A-Za-z]{2}', btop)

    for token in tokens:
        if token.isdigit():
            pos += int(token)
        else:
            # Cặp chữ: ref (query) → alt (subject)
            ref_char = token[0]
            alt_char = token[1]
            pos += 1  # vị trí của mismatch này (1-based trong query)
            mutations.append((pos, ref_char, alt_char))

    return mutations


# ── Hàm thống kê chung ──────────────────────────────────────
def count_mutations(input_file, has_header):
    """
    Đọc file, parse btop, thống kê:
      - mutation_count[pos][ref>alt] = số lần xuất hiện
      - seq_with_mut[seq_id]         = list đột biến
    """
    mutation_count = defaultdict(lambda: defaultdict(int))  # {pos: {mut: count}}
    seq_with_mut   = {}                                      # {seq_id: [mut_str]}
    total_seqs     = 0
    zero_mut_seqs  = 0

    with open(input_file, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            if has_header and i == 0:
                continue  # bỏ qua header

            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue

            seq_id   = cols[0]
            mismatch = int(cols[3])
            btop     = cols[4]
            total_seqs += 1

            if mismatch == 0:
                zero_mut_seqs += 1
                continue

            mutations = parse_btop(btop)
            mut_list  = []

            for (qpos, ref, alt) in mutations:
                genome_pos = qpos
                mut_label  = f"{ref}{genome_pos}{alt}"
                mutation_count[genome_pos][mut_label] += 1
                mut_list.append(mut_label)

            seq_with_mut[seq_id] = mut_list

    return mutation_count, seq_with_mut, total_seqs, zero_mut_seqs


# ── Hàm ghi kết quả ─────────────────────────────────────────
def write_stats(mutation_count, seq_with_mut, total_seqs, zero_mut_seqs,
                output_file, label):

    total_mut_seqs = len(seq_with_mut)

    with open(output_file, "w", encoding="utf-8") as f:

        # ── Phần 1: Tóm tắt tổng quan ──
        f.write("=" * 60 + "\n")
        f.write(f"  THỐNG KÊ ĐỘT BIẾN — {label}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Tổng số trình tự       : {total_seqs}\n")
        f.write(f"Trình tự không đột biến: {zero_mut_seqs}\n")
        f.write(f"Trình tự có đột biến   : {total_mut_seqs}\n")
        f.write(f"Số vị trí đột biến     : {len(mutation_count)}\n\n")

        # ── Phần 2: Tần suất mỗi đột biến ──
        f.write("-" * 60 + "\n")
        f.write(f"{'VỊ TRÍ':<10} {'ĐỘT BIẾN':<12} {'SỐ TRÌNH TỰ':>12} {'TẦN SUẤT (%)':>14}\n")
        f.write("-" * 60 + "\n")

        # Sort theo vị trí
        for pos in sorted(mutation_count.keys()):
            for mut_label, count in sorted(mutation_count[pos].items(),
                                           key=lambda x: -x[1]):
                freq = count / total_seqs * 100
                f.write(f"{pos:<10} {mut_label:<12} {count:>12} {freq:>13.2f}%\n")
        f.write("\n")

        # ── Phần 3: Chi tiết từng trình tự ──
        f.write("-" * 60 + "\n")
        f.write("CHI TIẾT TỪNG TRÌNH TỰ\n")
        f.write("-" * 60 + "\n")
        for seq_id, muts in sorted(seq_with_mut.items()):
            f.write(f"{seq_id}\t{', '.join(muts)}\n")

    print(f"[{label}] Tổng: {total_seqs} | Có đột biến: {total_mut_seqs} | "
          f"Vị trí đột biến: {len(mutation_count)} → {output_file}")


# ── Main ────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Thống kê đột biến từ BLAST btop")
    parser.add_argument("--nuc",        default=NUC_FILE,  help="File nucleotide blast")
    parser.add_argument("--prot",       default=PROT_FILE, help="File protein blast")
    parser.add_argument("--nuc_out",    default=NUC_OUT,   help="Output nucleotide stats")
    parser.add_argument("--prot_out",   default=PROT_OUT,  help="Output protein stats")
    args = parser.parse_args()

    # Nucleotide
    mc, sm, total, zero = count_mutations(args.nuc, has_header=False)
    write_stats(mc, sm, total, zero, args.nuc_out, label="NUCLEOTIDE")

    # Protein
    mc, sm, total, zero = count_mutations(args.prot, has_header=True)
    write_stats(mc, sm, total, zero, args.prot_out, label="PROTEIN")