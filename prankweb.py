"""
PrankWeb (Local P2Rank) and Fpocket binding site prediction module
"""

import os
import shutil
import subprocess
import glob
import re
import pandas as pd
import gradio as gr
from pathlib import Path
from config import current_pdb_info, PRANKWEB_OUTPUT_DIR, PREPARED_PROTEIN_DIR

# --- CONFIGURATION ---
P2RANK_HOME = "p2rank_2.5.1" 
P2RANK_CMD = "prank.bat" if os.name == 'nt' else "prank"
WSL_FPOCKET_BIN = "/home/drugproj/projects/fpocket/bin/fpocket"

def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    try:
        subprocess.run(['obabel', pdbqt_file, '-O', pdb_file], 
                       check=True, capture_output=True, shell=True)
        return True
    except Exception as e:
        print(f"OpenBabel conversion failed: {e}")
        return False

def to_wsl_path(win_path):
    try:
        drive = win_path.drive.lower().replace(':', '')
        path_str = win_path.relative_to(win_path.anchor).as_posix()
        return f"/mnt/{drive}/{path_str}"
    except Exception:
        return str(win_path).replace('\\', '/')

def extract_fpocket_coords(output_dir):
    """Scans 'pockets' subdir for coords."""
    pockets_dir = output_dir / "pockets"
    search_dir = pockets_dir if pockets_dir.exists() else output_dir
    
    coords_dict = {}
    report_lines = ["", "="*50, "📍 CALCULATED POCKET CENTERS (XYZ)", f"Source: {search_dir}", "="*50]
    
    pocket_files = list(search_dir.glob("pocket*_atm.pdb"))
    
    def natural_sort_key(p):
        match = re.search(r'pocket(\d+)_atm', p.name)
        return int(match.group(1)) if match else 0
    pocket_files.sort(key=natural_sort_key)

    if not pocket_files:
        return {}, "No pocket files found."

    for p_file in pocket_files:
        match = re.search(r'pocket(\d+)_atm', p_file.name)
        if not match: continue
        
        pocket_num = int(match.group(1))
        x_sum, y_sum, z_sum, count = 0.0, 0.0, 0.0, 0
        
        try:
            with open(p_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                            x_sum += x; y_sum += y; z_sum += z
                            count += 1
                        except ValueError: continue
            
            if count > 0:
                cx, cy, cz = x_sum/count, y_sum/count, z_sum/count
                coords_dict[pocket_num] = {'x': cx, 'y': cy, 'z': cz}
                report_lines.append(f"Pocket {pocket_num:<3}:  X={cx:8.3f}  Y={cy:8.3f}  Z={cz:8.3f}")
        except Exception: pass

    return coords_dict, "\n".join(report_lines)

def parse_fpocket_results(info_file, coords_dict):
    """Parses _info.txt and merges with coords. Applies 'fpocket_' prefix."""
    data = []
    if not info_file.exists(): return pd.DataFrame()

    with open(info_file, 'r') as f: content = f.read()

    pocket_blocks = re.split(r'(Pocket\s+\d+\s+:)', content)
    current_pocket_num = None
    
    for block in pocket_blocks:
        header_match = re.match(r'Pocket\s+(\d+)\s+:', block.strip())
        if header_match:
            current_pocket_num = int(header_match.group(1))
            continue
            
        if current_pocket_num is not None:
            score_match = re.search(r'Score\s+:\s+([-\d\.]+)', block)
            drug_match = re.search(r'Druggability Score\s+:\s+([-\d\.]+)', block)
            
            score = float(score_match.group(1)) if score_match else 0.0
            drug_score = float(drug_match.group(1)) if drug_match else 0.0
            coords = coords_dict.get(current_pocket_num, {'x': 0.0, 'y': 0.0, 'z': 0.0})
            
            # --- NAMING CONVENTION APPLIED HERE ---
            data.append({
                'name': f"fpocket_pocket{current_pocket_num}", 
                'rank': current_pocket_num, 
                'score': score,
                'probability': drug_score,
                'center_x': coords['x'],
                'center_y': coords['y'],
                'center_z': coords['z'],
                'method': 'Fpocket'
            })
            current_pocket_num = None

    return pd.DataFrame(data)

def run_fpocket_wsl(abs_input_path):
    abs_input = Path(abs_input_path).resolve()
    if not abs_input.exists(): return "Error: Input not found", pd.DataFrame()

    expected_output_dir = abs_input.parent / (abs_input.stem + "_out")
    wsl_input_path = to_wsl_path(abs_input)
    cmd = ["wsl", WSL_FPOCKET_BIN, "-f", wsl_input_path]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        if expected_output_dir.exists():
            info_file = expected_output_dir / f"{abs_input.stem}_info.txt"
            info_text = info_file.read_text() if info_file.exists() else "Info file missing."
            coords_dict, coords_report = extract_fpocket_coords(expected_output_dir)
            fpocket_df = parse_fpocket_results(info_file, coords_dict)
            return info_text + "\n" + coords_report, fpocket_df
        return "Output dir not found", pd.DataFrame()
    except Exception as e:
        return f"Fpocket Failed: {e}", pd.DataFrame()

def run_prankweb_prediction():
    """Run Binding Site prediction, MERGE results, and SAVE combined CSV."""
    
    prepared_pdbqt_path = current_pdb_info.get("prepared_pdbqt")
    if not prepared_pdbqt_path:
         potential_path = os.path.join(PREPARED_PROTEIN_DIR, "prepared_protein.pdbqt")
         if os.path.exists(potential_path): prepared_pdbqt_path = potential_path

    if not prepared_pdbqt_path or not os.path.exists(prepared_pdbqt_path):
        return (gr.update(value="❌ No prepared protein found.", visible=True), None, "")
    
    if os.path.exists(PRANKWEB_OUTPUT_DIR):
        try: shutil.rmtree(PRANKWEB_OUTPUT_DIR)
        except: pass
    os.makedirs(PRANKWEB_OUTPUT_DIR, exist_ok=True)

    yield (gr.update(value="⚙️ Converting PDBQT -> PDB..."), None, "Waiting...")
    
    temp_pdb_filename = "prepared_for_p2rank.pdb"
    temp_pdb_path = os.path.abspath(os.path.join(PRANKWEB_OUTPUT_DIR, temp_pdb_filename))
    if not convert_pdbqt_to_pdb(prepared_pdbqt_path, temp_pdb_path):
         return (gr.update(value="❌ Conversion Failed."), None, "Error")

    # --- P2RANK ---
    yield (gr.update(value="⚡ Running P2Rank..."), None, "P2Rank Running...")
    p2rank_df = pd.DataFrame()
    p2rank_err = ""
    try:
        subprocess.run([P2RANK_CMD, "predict", "-f", temp_pdb_path, "-o", os.path.abspath(PRANKWEB_OUTPUT_DIR)],
                       cwd=P2RANK_HOME, capture_output=True, check=True, shell=True)
        
        csv_files = glob.glob(os.path.join(PRANKWEB_OUTPUT_DIR, "*_predictions.csv"))
        if csv_files:
            df = pd.read_csv(csv_files[0])
            df.columns = df.columns.str.strip()
            
            # --- NAMING CONVENTION APPLIED HERE ---
            # Standard P2Rank names are 'pocket1', 'pocket2'. Replace 'pocket' with 'p2rank_pocket'
            df['name'] = df['name'].apply(lambda x: x.replace('pocket', 'p2rank_pocket') if 'pocket' in x else f"p2rank_{x}")
            
            df['method'] = 'P2Rank'
            p2rank_df = df
    except Exception as e:
        p2rank_err = str(e)

    # --- FPOCKET ---
    yield (gr.update(value="⚡ Running Fpocket..."), p2rank_df, "P2Rank Done. Fpocket Running...")
    fpocket_text, fpocket_df = run_fpocket_wsl(temp_pdb_path)

    # --- MERGE & SAVE ---
    final_df = pd.DataFrame()
    target_cols = ['name', 'rank', 'score', 'probability', 'center_x', 'center_y', 'center_z', 'method']
    
    if not p2rank_df.empty:
        if 'probability' not in p2rank_df.columns: p2rank_df['probability'] = 0.0
        p2rank_subset = p2rank_df[p2rank_df.columns.intersection(target_cols)].copy()
        final_df = pd.concat([final_df, p2rank_subset], ignore_index=True)

    if not fpocket_df.empty:
        final_df = pd.concat([final_df, fpocket_df], ignore_index=True)

    # *** CRITICAL STEP: Save the combined CSV for Docking ***
    combined_csv_path = os.path.join(PRANKWEB_OUTPUT_DIR, "combined_pockets.csv")
    final_df.to_csv(combined_csv_path, index=False)
    
    # Update Global Config so Docking knows where to look
    current_pdb_info["combined_csv"] = combined_csv_path 

    status_msg = f"✅ Predictions Complete. Merged {len(final_df)} pockets (Saved to combined_pockets.csv)."
    if p2rank_err: status_msg = f"❌ P2Rank Error: {p2rank_err}. Fpocket results attached."

    yield (
        gr.update(value=f"<div style='padding:15px; background:#d4edda; color:#155724;'>{status_msg}</div>", visible=True),
        gr.update(value=final_df, visible=True),
        gr.update(value=fpocket_text, visible=True)
    )