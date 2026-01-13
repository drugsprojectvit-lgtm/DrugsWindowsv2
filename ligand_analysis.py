# ligand_analysis.py
"""
Ligand Sanity Check & Analysis Module
Updated to be 'Lenient'.
1. Uses OpenBabel to extract SMILES (fixes aromaticity/bond issues).
2. Uses RDKit on the SMILES for properties (much safer than loading 3D files directly).
3. Uses OpenBabel to generate PDBs for visualization.
"""

import os
import glob
import pandas as pd
import subprocess
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
import gradio as gr
from config import LIGAND_DIR

# Suppress RDKit warnings just in case
RDLogger.DisableLog('rdApp.*') 

# Define output directory
ANALYSIS_OUTPUT_DIR = "ligand_analysis_results"
os.makedirs(ANALYSIS_OUTPUT_DIR, exist_ok=True)

def get_smiles_from_obabel(pdbqt_path):
    """
    Uses OpenBabel to convert PDBQT -> SMILES.
    OpenBabel is lenient and will 'fix' broken aromaticity (Kekulization) 
    that causes RDKit to crash.
    """
    try:
        # -osmi: Output SMILES
        # -h: Keep hydrogens (optional, but usually better to let RDKit re-add them)
        cmd = ['obabel', pdbqt_path, '-osmi']
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0 and result.stdout.strip():
            # Output format is usually "SMILES\tFilename"
            # We split by whitespace and take the first part
            smiles = result.stdout.strip().split()[0]
            return smiles
    except Exception as e:
        print(f"Error getting SMILES for {pdbqt_path}: {e}")
    return None

def pdbqt_to_pdb_obabel(pdbqt_path, output_pdb_path):
    """
    Uses OpenBabel to convert PDBQT -> PDB for the 3D viewer.
    """
    try:
        cmd = ['obabel', pdbqt_path, '-O', output_pdb_path]
        subprocess.run(cmd, check=False, capture_output=True)
        return os.path.exists(output_pdb_path)
    except:
        return False

def classify_molecule(mol, mw):
    """Classifies molecule based on MW and Ring info."""
    tags = []
    
    # MW Classification
    if mw < 250:
        tags.append("Fragment")
    elif 250 <= mw <= 500:
        tags.append("Small Molecule")
    else:
        tags.append("Large Molecule")

    # Macrocycle Detection
    try:
        ring_info = mol.GetRingInfo()
        is_macrocycle = any(len(ring) >= 12 for ring in ring_info.AtomRings())
        if is_macrocycle:
            tags.append("Macrocycle")
    except:
        pass 
        
    return ", ".join(tags)

def run_ligand_analysis():
    """Main execution function."""
    
    if not os.path.exists(LIGAND_DIR):
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Ligand directory not found.</div>", visible=True),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False)
        )

    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>🔬 extracting SMILES (Lenient Mode)...</div>", visible=True),
        gr.update(visible=False),
        gr.update(visible=False),
        gr.update(visible=False)
    )

    pdbqt_files = glob.glob(os.path.join(LIGAND_DIR, "*.pdbqt"))
    
    if not pdbqt_files:
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No PDBQT files found.</div>", visible=True),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(visible=False)
        )

    data = []
    seen_scaffolds = {} 
    seen_smiles = set()

    for file_path in pdbqt_files:
        filename = os.path.basename(file_path)
        
        # 1. Get SMILES via OpenBabel (The "Strictness" Bypass)
        smiles = get_smiles_from_obabel(file_path)
        
        if not smiles:
            data.append({"Filename": filename, "Status": "OpenBabel Conversion Failed", "Class": "-", "MW": 0, "Path": "N/A"})
            continue

        # 2. Load RDKit from SMILES (Much safer than from PDB/PDBQT)
        mol = Chem.MolFromSmiles(smiles)
        
        if not mol:
            # If even SMILES fails, we fallback to just showing the file
            data.append({"Filename": filename, "Status": "Invalid SMILES", "Class": "-", "MW": 0, "Path": "N/A"})
            continue
            
        # 3. Calculate Properties
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            rot_bonds = Descriptors.NumRotatableBonds(mol)
            formal_charge = Chem.GetFormalCharge(mol)
        except:
            mw, logp, rot_bonds, formal_charge = 0, 0, 0, 0
        
        # 4. Classification
        classification = classify_molecule(mol, mw)
        
        # 5. Scaffold Analysis
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smi = Chem.MolToSmiles(core)
            if not scaffold_smi: scaffold_smi = "Acyclic"
        except:
            scaffold_smi = "N/A"

        # 6. Redundancy Check
        redundancy_status = "Unique"
        if smiles in seen_smiles:
            redundancy_status = "Duplicate Structure"
        elif scaffold_smi in seen_scaffolds and scaffold_smi not in ["N/A", "Acyclic"]:
            redundancy_status = f"Common Scaffold with {seen_scaffolds[scaffold_smi]}"
        
        seen_smiles.add(smiles)
        if scaffold_smi not in seen_scaffolds and scaffold_smi not in ["N/A", "Acyclic"]:
            seen_scaffolds[scaffold_smi] = filename

        # 7. Generate PDB for visualization
        # We use Obabel to go PDBQT -> PDB directly.
        pdb_vis_path = os.path.join(ANALYSIS_OUTPUT_DIR, f"{filename}.pdb")
        pdbqt_to_pdb_obabel(file_path, pdb_vis_path)

        data.append({
            "Filename": filename,
            "Status": redundancy_status,
            "Class": classification,
            "MW": round(mw, 2),
            "Charge": formal_charge,
            "Rot. Bonds": rot_bonds,
            "LogP": round(logp, 2),
            "Scaffold SMILES": scaffold_smi,
            "Path": pdb_vis_path 
        })

    df = pd.DataFrame(data)
    csv_path = os.path.join(ANALYSIS_OUTPUT_DIR, "ligand_sanity_check.csv")
    df.to_csv(csv_path, index=False)
    
    choices = []
    for idx, row in df.iterrows():
        if row['Path'] != "N/A" and os.path.exists(row['Path']):
            label = f"{row['Filename']} [{row['Class']}]"
            choices.append((label, row['Path']))

    success_msg = f"""<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>
        ✅ Analysis Complete! Processed {len(df)} ligands (Lenient Mode).<br>
    </div>"""

    return (
        gr.update(value=success_msg, visible=True),
        gr.update(value=df, visible=True),
        gr.update(choices=choices, value=choices[0][1] if choices else None, visible=True),
        gr.update(value=csv_path, visible=True)
    )