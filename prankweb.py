"""
PrankWeb (Local P2Rank) and Fpocket binding site prediction module
"""

import os
import shutil  # Added for folder deletion
import subprocess
import glob
import re
import pandas as pd
import gradio as gr
from pathlib import Path
from config import current_pdb_info, PRANKWEB_OUTPUT_DIR, PREPARED_PROTEIN_DIR

# --- CONFIGURATION ---
# Path to your P2Rank folder (relative or absolute)
P2RANK_HOME = "p2rank_2.5.1" 
# On Windows use 'prank.bat', on Linux/Mac use 'prank'
P2RANK_CMD = "prank.bat" if os.name == 'nt' else "prank"

# --- FPOCKET CONFIGURATION (WSL) ---
WSL_FPOCKET_BIN = "/home/drugproj/projects/fpocket/bin/fpocket"

def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    """
    Converts PDBQT to PDB using OpenBabel (obabel).
    """
    try:
        # Using obabel to convert. -h adds hydrogens
        subprocess.run(['obabel', pdbqt_file, '-O', pdb_file], 
                       check=True, capture_output=True, shell=True)
        return True
    except Exception as e:
        print(f"OpenBabel conversion failed: {e}")
        return False

def to_wsl_path(win_path):
    """
    Converts a Windows Path object to a WSL path string.
    """
    try:
        drive = win_path.drive.lower().replace(':', '')
        # Use as_posix() to get forward slashes, relative_to anchor removes the drive letter/root
        path_str = win_path.relative_to(win_path.anchor).as_posix()
        return f"/mnt/{drive}/{path_str}"
    except Exception:
        # Fallback if path parsing fails, though Path objects usually handle this well
        return str(win_path).replace('\\', '/')

def calculate_pocket_centers(output_dir):
    """
    Scans the fpocket output directory (specifically the 'pockets' subdir)
    for 'pocket*_atm.pdb' files, calculates the geometric center (centroid),
    and returns a formatted string report.
    """
    # FIX: Fpocket creates a 'pockets' subdirectory inside the output folder.
    # We must look there first.
    pockets_dir = output_dir / "pockets"
    
    if pockets_dir.exists():
        search_dir = pockets_dir
        location_msg = f"Scanning directory: {pockets_dir}"
    else:
        # Fallback to root if 'pockets' doesn't exist (older versions)
        search_dir = output_dir
        location_msg = f"Scanning directory: {output_dir} ('pockets' subdir not found)"

    report_lines = ["", "="*50, "📍 CALCULATED POCKET CENTERS (XYZ)", location_msg, "="*50]
    
    # Find all pocket pdb files (e.g., pocket1_atm.pdb)
    pocket_files = list(search_dir.glob("pocket*_atm.pdb"))
    
    # Sort them naturally (1, 2, 10 instead of 1, 10, 2)
    def natural_sort_key(p):
        # Extract the number from filename 'pocket12_atm.pdb'
        match = re.search(r'pocket(\d+)_atm', p.name)
        return int(match.group(1)) if match else 0

    pocket_files.sort(key=natural_sort_key)

    if not pocket_files:
        report_lines.append("❌ No 'pocket*_atm.pdb' files found in the scan directory.")
        report_lines.append(f"Contents of {search_dir.name}: {os.listdir(search_dir) if search_dir.exists() else 'Directory does not exist'}")
        return "\n".join(report_lines)

    for p_file in pocket_files:
        x_sum, y_sum, z_sum, count = 0.0, 0.0, 0.0, 0
        try:
            with open(p_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        # PDB format: X is 30-38, Y is 38-46, Z is 46-54
                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            x_sum += x
                            y_sum += y
                            z_sum += z
                            count += 1
                        except ValueError:
                            continue
            
            if count > 0:
                cx = x_sum / count
                cy = y_sum / count
                cz = z_sum / count
                # Extract simple name like "Pocket 1"
                match = re.search(r'pocket(\d+)_atm', p_file.name)
                pocket_num = match.group(1) if match else "Unknown"
                
                # Format: Pocket 1:  X= 12.345  Y= 45.678  Z= -9.123
                report_lines.append(f"Pocket {pocket_num:<3}:  X={cx:8.3f}  Y={cy:8.3f}  Z={cz:8.3f}")
            else:
                report_lines.append(f"{p_file.name}: Empty or invalid atoms.")

        except Exception as e:
            report_lines.append(f"{p_file.name}: Error reading file ({str(e)})")

    return "\n".join(report_lines)

def run_fpocket_wsl(abs_input_path):
    """
    Runs fpocket via WSL on the specific PDB file.
    Returns the content of the _info.txt file PLUS calculated XYZ coordinates.
    """
    abs_input = Path(abs_input_path).resolve()
    
    if not abs_input.exists():
        return f"❌ Fpocket Error: Input file not found at {abs_input}"

    # Expected output directory logic from fpocket (filename_out)
    expected_output_dir = abs_input.parent / (abs_input.stem + "_out")
    
    # Convert to WSL path
    wsl_input_path = to_wsl_path(abs_input)

    command = [
        "wsl",
        WSL_FPOCKET_BIN,
        "-f", wsl_input_path
    ]

    print(f"🔄 Fpocket Processing: {abs_input.name}...")

    try:
        result = subprocess.run(
            command, 
            capture_output=True, 
            text=True, 
            check=True
        )

        if expected_output_dir.exists():
            # 1. Read the standard Info file
            # This file is usually in the root of the output folder (not in pockets subdir)
            info_file_name = f"{abs_input.stem}_info.txt"
            info_file_path = expected_output_dir / info_file_name
            
            main_content = ""
            if info_file_path.exists():
                with open(info_file_path, 'r') as f:
                    main_content = f.read()
            else:
                main_content = f"⚠️ Fpocket ran, but info file not found at: {info_file_path}\n"
            
            # 2. Calculate XYZ Coordinates from generated PDBs (inside 'pockets' folder)
            coordinates_report = calculate_pocket_centers(expected_output_dir)
            
            # Combine them
            full_report = main_content + "\n" + coordinates_report
            return full_report

        else:
            return f"⚠️ Fpocket command finished, but output folder not found.\nExpected: {expected_output_dir}\nWSL Stdout: {result.stdout}"

    except subprocess.CalledProcessError as e:
        return f"❌ Fpocket Execution Failed.\nError details: {e.stderr}"
    except Exception as ex:
        return f"❌ System Error running Fpocket: {str(ex)}"


def run_prankweb_prediction():
    """Run Binding Site prediction using LOCAL P2Rank AND Fpocket (WSL)."""
    
    # 1. Validate Input (Must come from Protein Prep step)
    prepared_pdbqt_path = current_pdb_info.get("prepared_pdbqt")
    
    # Fallback checks
    if not prepared_pdbqt_path:
         potential_path = os.path.join(PREPARED_PROTEIN_DIR, "prepared_protein.pdbqt")
         if os.path.exists(potential_path):
             prepared_pdbqt_path = potential_path

    if not prepared_pdbqt_path or not os.path.exists(prepared_pdbqt_path):
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No prepared protein found. Please run <b>Protein Preparation</b> first.</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(value="", visible=False)
        )
    
    # Check if P2Rank exists
    if not os.path.exists(P2RANK_HOME):
         return (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ P2Rank folder not found at: {P2RANK_HOME}. Please check the path.</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(value="", visible=False)
        )

    # --- CLEANUP LOGIC ---
    # Delete the output directory if it exists to ensure a fresh run
    if os.path.exists(PRANKWEB_OUTPUT_DIR):
        try:
            shutil.rmtree(PRANKWEB_OUTPUT_DIR)
            print(f"🗑️ Deleted previous results folder: {PRANKWEB_OUTPUT_DIR}")
        except Exception as e:
            print(f"⚠️ Warning: Could not delete output directory: {e}")

    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>⚙️ Converting PDBQT to PDB...</div>", visible=True),
        gr.update(value=None, visible=False),
        gr.update(value="Waiting...", visible=True)
    )
    
    output_dir = PRANKWEB_OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)
    
    # 2. Convert PDBQT -> PDB (Required for both P2Rank and Fpocket)
    temp_pdb_filename = "prepared_for_p2rank.pdb"
    # Use ABSOLUTE paths because we are changing CWD later
    temp_pdb_path = os.path.abspath(os.path.join(output_dir, temp_pdb_filename))
    abs_output_dir = os.path.abspath(output_dir)
    
    success = convert_pdbqt_to_pdb(prepared_pdbqt_path, temp_pdb_path)
    
    if not success:
         return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Failed to convert PDBQT to PDB. Please ensure <b>OpenBabel</b> is installed.</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(value="Conversion Failed", visible=True)
        )

    # --- RUNNING P2RANK ---
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>⚡ Running local P2Rank prediction...</div>", visible=True),
        gr.update(value=None, visible=False),
        gr.update(value="P2Rank Running...", visible=True)
    )

    p2rank_df = None
    p2rank_error = None

    try:
        # Command: prank predict -f "file.pdb" -o "folder"
        # We use cwd=P2RANK_HOME to simulate "cd into folder"
        cmd = [
            P2RANK_CMD,
            "predict",
            "-f", temp_pdb_path,
            "-o", abs_output_dir
        ]
        
        # Run subprocess inside the p2rank folder
        process = subprocess.run(
            cmd, 
            cwd=P2RANK_HOME, 
            capture_output=True, 
            text=True, 
            shell=True 
        )

        if process.returncode != 0:
            raise Exception(f"P2Rank exited with error: {process.stderr}")

        # Process Results
        expected_csv_name = f"{os.path.splitext(temp_pdb_filename)[0]}_predictions.csv"
        csv_path = os.path.join(output_dir, expected_csv_name)

        if not os.path.exists(csv_path):
             found_csvs = glob.glob(os.path.join(output_dir, "*_predictions.csv"))
             if found_csvs:
                 csv_path = found_csvs[0]
             else:
                 raise FileNotFoundError("Could not find _predictions.csv in output folder.")

        # Read and filter CSV
        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip()
        columns_to_drop = ['residue_ids', 'surf_atom_ids']
        p2rank_df = df.drop(columns=[col for col in columns_to_drop if col in df.columns], errors='ignore')
        
        # Store CSV path globally for later use in docking
        current_pdb_info["prankweb_csv"] = csv_path
        
    except Exception as e:
        p2rank_error = str(e)

    # --- RUNNING FPOCKET (WSL) ---
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>⚡ P2Rank Done. Running Fpocket (WSL)...</div>", visible=True),
        gr.update(value=p2rank_df, visible=(p2rank_df is not None)),
        gr.update(value="P2Rank Done. Starting Fpocket...", visible=True)
    )

    fpocket_text = run_fpocket_wsl(temp_pdb_path)

    # --- FINAL RETURN ---
    if p2rank_error:
        final_status = f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error executing P2Rank: {p2rank_error}. Fpocket ran (see below).</div>"
    else:
        final_status = f"<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>✅ P2Rank & Fpocket predictions completed!</div>"

    yield (
        gr.update(value=final_status, visible=True),
        gr.update(value=p2rank_df, visible=(p2rank_df is not None)),
        gr.update(value=fpocket_text, visible=True)
    )