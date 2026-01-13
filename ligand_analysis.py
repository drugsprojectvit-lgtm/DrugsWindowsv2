# ligand_analysis.py
import os
import pandas as pd
import subprocess

def analyze_pdbqt_manual(pdbqt_path):
    """Parses PDBQT to calculate Molecular Weight and Rotatable Bonds."""
    atom_masses = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'P': 30.974, 
        'S': 32.065, 'F': 18.998, 'CL': 35.45, 'BR': 79.904, 'I': 126.90
    }
    total_mass = 0.0
    rotatable_bonds = 0
    found_torsdof = False

    with open(pdbqt_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    parts = line.split()
                    atom_type = parts[-1].upper()
                    element = atom_type if atom_type in atom_masses else atom_type[0]
                    if element in atom_masses:
                        total_mass += atom_masses[element]
                except: continue

            if line.startswith("REMARK") and "TORSDOF" in line:
                try:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "TORSDOF": 
                            rotatable_bonds = int(parts[i+1])
                    found_torsdof = True
                except: pass

            if not found_torsdof and line.startswith("BRANCH"):
                rotatable_bonds += 1

    return {"MW": total_mass, "RotBonds": rotatable_bonds}

def convert_pdbqt_to_pdb(pdbqt_path, output_folder):
    """Converts a PDBQT file to PDB format using obabel."""
    try:
        base_name = os.path.basename(pdbqt_path)
        pdb_name = base_name.replace(".pdbqt", ".pdb")
        output_path = os.path.join(output_folder, pdb_name)
        
        # subprocess command for obabel
        subprocess.run(
            ["obabel", "-ipdbqt", pdbqt_path, "-opdb", "-O", output_path],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        return pdb_name
    except Exception as e:
        print(f"Error converting {pdbqt_path}: {e}")
        return None

def run_ligand_classification(ligand_folder="pdbqt"):
    """Analyzes PDBQT files and converts them to PDB for visualization."""
    if not os.path.exists(ligand_folder):
        return None, "Folder not found."

    # Create folder for converted PDBs
    pdb_out_folder = "ligand_pdb"
    os.makedirs(pdb_out_folder, exist_ok=True)

    results = []
    pdbqt_files = [f for f in os.listdir(ligand_folder) if f.endswith(".pdbqt")]
    
    if not pdbqt_files:
        return None, "No PDBQT files found in folder."

    for file in pdbqt_files:
        path = os.path.join(ligand_folder, file)
        
        # 1. Analyze properties (using original PDBQT)
        props = analyze_pdbqt_manual(path)
        mw = props['MW']
        rb = props['RotBonds']
        
        # 2. Convert to PDB for visualization
        converted_pdb = convert_pdbqt_to_pdb(path, pdb_out_folder)
        
        # Logic
        if mw < 300:
            cat, guide = "Fragment", "Good for fragment-based screening."
        elif 300 <= mw <= 600:
            if rb < 12:
                cat, guide = "Small Molecule", "Optimal for Standard Docking."
            else:
                cat, guide = "Flexible Small Molecule", "High flexibility; caution with accuracy."
        else:
            if rb > 15:
                cat, guide = "Peptide / Biologic", "Use specialized Vina protocols."
            else:
                cat, guide = "Macrocycle", "Check ring conformation sampling."

        results.append({
            "Ligand": file,
            "Converted_PDB": converted_pdb,
            "MW (Da)": round(mw, 2),
            "Rotatable Bonds": rb,
            "Category": cat,
            "Guidance": guide
        })

    df = pd.DataFrame(results)
    output_csv = "ligand_classification_report.csv"
    df.to_csv(output_csv, index=False)
    
    return df, output_csv