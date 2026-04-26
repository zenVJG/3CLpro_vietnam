import os
import subprocess
import csv
import re
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

NUM_RUNS = 3  # ← Thay số lần chạy tại đây

def run_vina(receptor, ligand, config, output_path):
    vina_path = r'D:\Vina\vina.exe'
    cmd = [
        vina_path,
        "--receptor", receptor,
        "--ligand",   ligand,
        "--config",   config,
        "--out",      output_path,
    ]
    
    subprocess.run(cmd, check=True, capture_output=True, text=True)

def extract_best_affinity(pdbqt_file):
    best_affinity = None
    with open(pdbqt_file, 'r', encoding='utf-8') as f:
        for line in f:
            match = re.search(r'REMARK VINA RESULT:\s+(-?\d+\.\d+)', line)
            if match:
                affinity = float(match.group(1))
                if best_affinity is None or affinity < best_affinity:
                    best_affinity = affinity
    return best_affinity if best_affinity is not None else "No Result"

def clean_empty_outputs(output_folder):
    for f in os.listdir(output_folder):
        path = os.path.join(output_folder, f)
        if os.path.getsize(path) == 0:
            os.remove(path)

def process_ligand(args):
    receptor, ligand, config_file, output_folder, num_runs = args
    ligand_name = os.path.basename(ligand)
    stem = os.path.splitext(ligand_name)[0]  # tên không có đuôi .pdbqt

    row = [ligand_name]
    for run in range(1, num_runs + 1):
        # Mỗi run lưu file riêng: baicalin_run1.pdbqt, baicalin_run2.pdbqt...
        output_name = f"{stem}_run{run}.pdbqt"
        output_path = os.path.join(output_folder, output_name)

        if os.path.exists(output_path):
            affinity = extract_best_affinity(output_path)
            print(f"  ⏭️  Skipping {output_name}: {affinity} kcal/mol")
            row.append(affinity)
            continue

        try:
            run_vina(receptor, ligand, config_file, output_path)
            affinity = extract_best_affinity(output_path)
            print(f"  ✅ {output_name}: {affinity} kcal/mol")
            row.append(affinity)

        except subprocess.CalledProcessError as e:
            print(f"  ❌ Failed: {output_name}")
            print(f"     STDERR: {e.stderr.strip()}")
            row.append("Failed")

    return row  # [ligand_name, run1, run2, run3]

def dock_protein(receptor_path, config_path, ligands_folder, output_folder, csv_output, num_workers, num_runs):
    os.makedirs(output_folder, exist_ok=True)
    clean_empty_outputs(output_folder)

    ligands = [
        os.path.join(ligands_folder, f)
        for f in os.listdir(ligands_folder)
        if f.endswith(".pdbqt")
    ]

    if not ligands:
        print("  ❌ Không tìm thấy ligand nào!")
        return

    print(f"  🔬 {len(ligands)} ligands × {num_runs} runs = {len(ligands) * num_runs} jobs\n")

    args = [(receptor_path, lig, config_path, output_folder, num_runs) for lig in ligands]

    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for row in executor.map(process_ligand, args):
            results.append(row)

    if results:
        write_header = not os.path.exists(csv_output) or os.path.getsize(csv_output) == 0
        with open(csv_output, 'a', newline='') as f:
            writer = csv.writer(f)
            if write_header:
                # Header: Ligand | Run 1 | Run 2 | Run 3 | Best
                header = ["Ligand"] + [f"Run {i} (kcal/mol)" for i in range(1, num_runs + 1)] + ["Best (kcal/mol)"]
                writer.writerow(header)
            for row in results:
                # Tính best affinity trong tất cả các run
                affinities = [r for r in row[1:] if isinstance(r, float) or (isinstance(r, str) and r.replace('.','').replace('-','').isdigit())]
                best = min(float(a) for a in affinities) if affinities else "N/A"
                writer.writerow(row + [best])
        print(f"\n  📄 Saved: {csv_output}")

def run_all(proteins_folder, configs_folder, ligands_folder, results_folder, num_workers, num_runs):
    proteins = sorted([f for f in os.listdir(proteins_folder) if f.endswith(".pdbqt")])

    if not proteins:
        print("❌ Không tìm thấy protein nào!")
        return

    print(f"🧬 {len(proteins)} proteins | 🔁 {num_runs} runs/ligand\n{'='*45}")

    for i, protein_file in enumerate(proteins, 1):
        protein_name  = os.path.splitext(protein_file)[0]
        receptor_path = os.path.join(proteins_folder, protein_file)
        config_path   = os.path.join(configs_folder,  f"{protein_name}.txt")
        output_folder = os.path.join(results_folder,  protein_name)
        csv_output    = os.path.join(results_folder,  f"{protein_name}.csv")

        print(f"\n[{i}/{len(proteins)}] 🧪 {protein_name}")

        if not os.path.exists(config_path):
            print(f"  ⚠️  Thiếu config: {config_path} → Bỏ qua")
            continue

        dock_protein(receptor_path, config_path, ligands_folder, output_folder, csv_output, num_workers, num_runs)

    print(f"\n{'='*45}\n🎉 Hoàn tất!")

# -------- Cấu hình --------
proteins_folder = r"D:\auto_docking\proteins"
configs_folder  = r"D:\auto_docking\configs"
ligands_folder  = r"D:\auto_docking\ligands"
results_folder  = r"D:\auto_docking\results"
num_workers     = max(4, os.cpu_count() // 2)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    run_all(proteins_folder, configs_folder, ligands_folder, results_folder, num_workers, NUM_RUNS)