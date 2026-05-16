#!/c/Users/Ws/AppData/Local/Programs/Python/Python312/python.exe
"""Batch docking: run Vina for missing ligands only"""
import os, glob, subprocess, time

BASE = 'C:/Users/Ws/Desktop/vina_batch_package'
LIG_DIR = os.path.join(BASE, 'ligands')
OUT_DIR = os.path.join(BASE, 'results')
VINA = os.path.join(BASE, 'vina_1.2.5_win.exe')
RECEPTOR = os.path.join(BASE, '7PZC_NACHT.pdb')

# Find all ligand files
all_ligands = sorted(glob.glob(os.path.join(LIG_DIR, 'lig_*.pdbqt')))
existing_out = set(os.path.basename(f) for f in glob.glob(os.path.join(OUT_DIR, 'out_lig_*.pdbqt')))

# Find missing
missing = []
for lig_path in all_ligands:
    out_name = 'out_' + os.path.basename(lig_path)
    if out_name not in existing_out:
        missing.append(lig_path)

print(f'Total ligands: {len(all_ligands)}')
print(f'Existing results: {len(existing_out)}')
print(f'Missing: {len(missing)}')
if missing:
    print(f'First missing: {os.path.basename(missing[0])}')
    print(f'Last missing: {os.path.basename(missing[-1])}')
    
    # Run Vina for each missing
    total = len(missing)
    for i, lig_path in enumerate(missing, 1):
        lig_name = os.path.basename(lig_path)
        out_path = os.path.join(OUT_DIR, 'out_' + lig_name)
        
        cmd = [VINA, '--receptor', RECEPTOR, '--ligand', lig_path,
               '--center_x', '192.9', '--center_y', '204.7', '--center_z', '119.7',
               '--size_x', '20', '--size_y', '20', '--size_z', '20',
               '--out', out_path, '--cpu', '4', '--exhaustiveness', '8']
        
        print(f'[{i}/{total}] {lig_name} ... ', end='', flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print(f'OK ({i}/{total})')
        else:
            print(f'FAILED')
else:
    print('All done!')
