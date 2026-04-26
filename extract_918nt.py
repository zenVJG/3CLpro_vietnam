#!/usr/bin/env python3
"""
Extract 918nt regions using blastdbcmd (Windows compatible)
Usage:
    python extract_918nt.py

Chỉnh 3 biến DB_PATH, BLAST_FILE, OUT_FILE bên dưới rồi chạy.
"""

import subprocess

# ============================================================
# CHỈNH 3 BIẾN NÀY TRƯỚC KHI CHẠY
DB_PATH    = "nucl"       # đường dẫn tới BLAST database
BLAST_FILE = "filtered_nucl.txt"  # file kết quả BLAST
OUT_FILE   = "extracted_918nt.fasta"  # output
# ============================================================

def main():
    found   = 0
    missing = []

    with open(BLAST_FILE, encoding="utf-8") as f_in, \
         open(OUT_FILE, "w", encoding="utf-8") as f_out:

        for line in f_in: # Tách dòng
            line = line.strip()
            if not line:
                continue

            cols   = line.split("\t")
            seqid  = cols[0]
            sstart = int(cols[5])
            send   = int(cols[6])

            # Xác định strand
            if sstart <= send:
                coord  = f"{sstart}-{send}"
                strand = "plus"
            else:
                coord  = f"{send}-{sstart}"
                strand = "minus"

            # Gọi blastdbcmd
            cmd = [
                "blastdbcmd",
                "-db",     DB_PATH,
                "-entry",  seqid,
                "-range",  coord,
                "-strand", strand,
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0 or not result.stdout.strip(): # Loại bỏ dữ liệu trống
                missing.append(seqid)
                continue

            # Thay header gốc bằng header có tọa độ
            lines = result.stdout.strip().split("\n")
            seq   = "".join(lines[1:])  # bỏ dòng header gốc
            f_out.write(f">{seqid}\n{seq}\n")
            found += 1

            if found % 100 == 0:
                print(f"  Processed: {found} sequences...")

    print(f"\nDone!")
    print(f"  Extracted : {found} sequences")
    print(f"  Not found : {len(missing)} sequences")

    if missing:
        missing_file = OUT_FILE.replace(".fasta", "_missing.txt")
        with open(missing_file, "w", encoding="utf-8") as mf:
            mf.writelines(f"{m}\n" for m in missing)
        print(f"  -> Missing IDs saved to: {missing_file}")

if __name__ == "__main__":
    main()