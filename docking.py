# docking.py
"""
Molecular docking module using AutoDock Vina Executable
Supports splitting receptor into chains.
"""

import os
import glob
import subprocess
import re
import shutil
import pandas as pd
import gradio as gr
from config import current_pdb_info, DOCKING_RESULTS_DIR, LIGAND_DIR, PRANKWEB_OUTPUT_DIR

# Path to your Vina executable
VINA_EXE = "vina.exe" 

def split_pdbqt_chains(pdbqt_path, output_base_dir):
    """
    Splits a prepared PDBQT file into separate PDBQT files based on TER records.
    Returns a dictionary: {'Chain_A': path, 'Chain_B': path, ...}
    """
    chains = {}
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    chain_idx = 0
    
    current_lines = []
    
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if line.startswith("TER") or line.startswith("ENDMDL"):
            if current_lines:
                # We have a chain
                chain_name = f"Chain_{chain_letters[chain_idx % len(chain_letters)]}"
                if chain_idx >= len(chain_letters):
                    chain_name += str(chain_idx // len(chain_letters))
                
                chain_dir = os.path.join(output_base_dir, chain_name)
                os.makedirs(chain_dir, exist_ok=True)
                
                chain_file_path = os.path.join(chain_dir, "receptor.pdbqt")
                
                with open(chain_file_path, 'w') as out:
                    out.writelines(current_lines)
                
                chains[chain_name] = chain_file_path
                current_lines = []
                chain_idx += 1
        else:
            current_lines.append(line)
            
    # Handle last chain if no trailing TER/ENDMDL
    if current_lines:
        chain_name = f"Chain_{chain_letters[chain_idx % len(chain_letters)]}"
        chain_dir = os.path.join(output_base_dir, chain_name)
        os.makedirs(chain_dir, exist_ok=True)
        chain_file_path = os.path.join(chain_dir, "receptor.pdbqt")
        with open(chain_file_path, 'w') as out:
            out.writelines(current_lines)
        chains[chain_name] = chain_file_path
        
    return chains

def run_molecular_docking():
    """Run molecular docking on EACH CHAIN individually."""
    
    if not current_pdb_info.get("prepared_pdbqt"):
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No prepared protein found. Please prepare protein first.</div>", visible=True),
            gr.update(value=None, visible=False),
            gr.update(choices=[], visible=False), # Chain selector
            gr.update(choices=[], visible=False)  # Pose selector
        )
    
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>🔬 Splitting chains and running docking (this may take several minutes)...</div>", visible=True),
        gr.update(value=None, visible=False),
        gr.update(choices=[], visible=False),
        gr.update(choices=[], visible=False)
    )
    
    try:
        # --- CLEANUP OLD RESULTS ---
        if os.path.exists(DOCKING_RESULTS_DIR):
            try:
                shutil.rmtree(DOCKING_RESULTS_DIR)
            except Exception as e:
                yield (gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error deleting old results: {str(e)}</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))
                return
        os.makedirs(DOCKING_RESULTS_DIR, exist_ok=True)

        # Input files
        protein_pdbqt = current_pdb_info["prepared_pdbqt"]
        ligand_folder = LIGAND_DIR
        
        # --- SPLIT CHAINS ---
        chain_map = split_pdbqt_chains(protein_pdbqt, DOCKING_RESULTS_DIR)
        if not chain_map:
             yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Could not split protein into chains.</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))
             return
        
        # Find PrankWeb CSV
        output_dir = PRANKWEB_OUTPUT_DIR
        csv_file = None
        if current_pdb_info.get("prankweb_csv") and os.path.exists(current_pdb_info["prankweb_csv"]):
            csv_file = current_pdb_info["prankweb_csv"]
        else:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    if file.endswith("_predictions.csv"):
                        csv_file = os.path.join(root, file)
                        break
                if csv_file: break
        
        if not csv_file or not os.path.exists(csv_file):
            yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ PrankWeb results not found.</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))
            return
        
        if not os.path.exists(ligand_folder):
            yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Ligand folder not found.</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))
            return
        
        # Load CSV
        df = pd.read_csv(csv_file)
        df.columns = df.columns.str.strip()
        required_cols = ['name', 'center_x', 'center_y', 'center_z']
        if any(col not in df.columns for col in required_cols):
            yield (gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ CSV format error. Missing columns.</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))
            return

        ligand_files = glob.glob(os.path.join(ligand_folder, "*.pdbqt"))
        if not ligand_files:
            yield (gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No ligand files found.</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))
            return
        
        summary_data = []

        # --- DOCKING LOOP PER CHAIN ---
        for chain_id, chain_receptor_pdbqt in chain_map.items():
            
            # Create subdirs for this chain
            chain_base_dir = os.path.dirname(chain_receptor_pdbqt)
            output_dir_pdbqt = os.path.join(chain_base_dir, "docked_pdbqt")
            output_dir_pdb = os.path.join(chain_base_dir, "docked_pdb")
            os.makedirs(output_dir_pdbqt, exist_ok=True)
            os.makedirs(output_dir_pdb, exist_ok=True)

            # Convert Chain Receptor PDBQT to PDB for visualization later
            chain_receptor_pdb = os.path.join(output_dir_pdb, "receptor_ref.pdb")
            subprocess.run(['obabel', chain_receptor_pdbqt, '-O', chain_receptor_pdb], check=False, capture_output=True)

            for ligand_pdbqt in ligand_files:
                ligand_name = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
                ligand_best_poses = []
                
                for index, row in df.iterrows():
                    pocket_name = str(row['name']).strip()
                    cx, cy, cz = float(row['center_x']), float(row['center_y']), float(row['center_z'])
                    
                    output_pdbqt_file = os.path.join(output_dir_pdbqt, f"{ligand_name}_{pocket_name}_poses.pdbqt")

                    cmd = [
                        VINA_EXE, "--receptor", chain_receptor_pdbqt, "--ligand", ligand_pdbqt,
                        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                        "--size_x", "25", "--size_y", "25", "--size_z", "25",
                        "--exhaustiveness", "8", "--num_modes", "10", "--out", output_pdbqt_file
                    ]

                    try:
                        result = subprocess.run(cmd, capture_output=True, text=True, check=True, shell=True)
                        
                        # Parse Energies
                        for line in result.stdout.splitlines():
                            if re.match(r'^\s*\d+', line):
                                parts = line.split()
                                if len(parts) >= 2:
                                    try:
                                        mode_num = int(parts[0])
                                        affinity = float(parts[1])
                                        
                                        if affinity <= 0.0:
                                            ligand_best_poses.append({
                                                'chain': chain_id,
                                                'ligand': ligand_name,
                                                'pocket': pocket_name,
                                                'pose_number': mode_num,
                                                'binding_energy': affinity,
                                                'pdb_file': os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_ligand.pdb"),
                                                'receptor_pdb_file': chain_receptor_pdb, # Store ref to specific chain PDB
                                                'interaction_image': "N/A" 
                                            })
                                    except ValueError: continue

                        # Post-processing (Images/Complexes)
                        if os.path.exists(output_pdbqt_file):
                            pdb_ligand_only = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_ligand.pdb")
                            subprocess.run(['obabel', output_pdbqt_file, '-O', pdb_ligand_only, '-h'], check=False, capture_output=True)

                            # Generate Complex for PANDAMAP
                            complex_file = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_complex.pdb")
                            
                            # Read Receptor (The specific Chain)
                            rec_lines = []
                            if os.path.exists(chain_receptor_pdb):
                                with open(chain_receptor_pdb, 'r') as f:
                                    rec_lines = [l for l in f if l.startswith(("ATOM", "HETATM", "TER"))]

                            # Read Ligand (Model 1)
                            lig_lines = []
                            with open(output_pdbqt_file, 'r') as lf:
                                in_model_1 = False
                                for line in lf:
                                    if line.startswith("MODEL 1"): in_model_1 = True; continue
                                    if line.startswith("ENDMDL"): break
                                    if in_model_1 and line.startswith("ATOM"):
                                        lig_lines.append("HETATM" + line[6:])

                            with open(complex_file, 'w') as cf:
                                cf.writelines(rec_lines)
                                if rec_lines and not rec_lines[-1].strip() == "TER": cf.write("TER\n")
                                cf.writelines(lig_lines)
                                cf.write("END\n")

                            # Run Pandamap
                            if os.path.exists(complex_file):
                                interactions_png = os.path.join(output_dir_pdb, f"{ligand_name}_{pocket_name}_inter.png")
                                try:
                                    subprocess.run(f"pandamap {complex_file} --output {interactions_png}", shell=True, check=True, capture_output=True)
                                    if os.path.exists(interactions_png):
                                        for entry in ligand_best_poses:
                                            if entry['ligand'] == ligand_name and entry['pocket'] == pocket_name:
                                                entry['interaction_image'] = interactions_png
                                except: pass

                    except subprocess.CalledProcessError as e:
                        print(f"Failed {chain_id}/{ligand_name}: {e}")
                        continue
                
                if ligand_best_poses:
                    ligand_best_poses.sort(key=lambda x: x['binding_energy'])
                    summary_data.extend(ligand_best_poses[:3]) # Top 3 per ligand per chain

        # --- Output Logic ---
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            
            # Get list of unique chains
            unique_chains = sorted(summary_df['chain'].unique())
            
            # Default: Select first chain
            default_chain = unique_chains[0]
            
            # Filter poses for first chain to populate initial list
            chain_df = summary_df[summary_df['chain'] == default_chain]
            pose_choices = []
            for idx, row in chain_df.iterrows():
                label = f"{row['ligand']} - {row['pocket']} (Pose {row['pose_number']}) | {row['binding_energy']:.2f} kcal"
                value = f"{row['pdb_file']}::{row['pose_number']}::{row['chain']}" # Added chain to value for safety
                pose_choices.append((label, value))

            success_msg = f"""<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>
                ✅ Docking completed across {len(unique_chains)} chains! Found {len(summary_data)} total poses.
            </div>"""
            
            yield (
                gr.update(value=success_msg, visible=True),
                gr.update(value=summary_df, visible=True),
                gr.update(choices=unique_chains, visible=True, value=default_chain), # Update Chain Dropdown
                gr.update(choices=pose_choices, visible=True, value=pose_choices[0][1] if pose_choices else None) # Update Pose Dropdown
            )
        else:
            yield (
                gr.update(value="<div style='padding: 20px; background: #f8d7da; border-radius: 8px; color: #721c24;'>⚠️ No good poses found.</div>", visible=True),
                gr.update(value=None, visible=False),
                gr.update(choices=[], visible=False),
                gr.update(choices=[], visible=False)
            )
    
    except Exception as e:
        yield (gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>", visible=True), None, gr.update(choices=[]), gr.update(choices=[]))