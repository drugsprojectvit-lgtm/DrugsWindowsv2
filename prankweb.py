"""
PrankWeb (Local P2Rank) and Fpocket binding site prediction module
"""

import os
import subprocess
import glob
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

def run_fpocket_wsl(abs_input_path):
    """
    Runs fpocket via WSL on the specific PDB file.
    Returns the content of the _info.txt file if successful, or an error message.
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
            # Look for the info.txt file
            # usually named: [filename]_info.txt inside the [filename]_out folder
            info_file_name = f"{abs_input.stem}_info.txt"
            info_file_path = expected_output_dir / info_file_name
            
            if info_file_path.exists():
                with open(info_file_path, 'r') as f:
                    content = f.read()
                return content
            else:
                return f"⚠️ Fpocket ran, but info file not found at: {info_file_path}\nOutput dir content: {os.listdir(expected_output_dir)}"
        else:
            return f"⚠️ Fpocket command finished, but output folder not found.\nWSL Stdout: {result.stdout}"

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