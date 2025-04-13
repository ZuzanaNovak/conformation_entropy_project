import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

sdf_folder = "ligands"
csv_path = "CREST/results/tables/007_JAK1_rotatable_bonds.csv"

df = pd.read_csv(csv_path)

def analyze_molecule(mol):
    Chem.SanitizeMol(mol)

    # Size
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    mol_weight = Descriptors.MolWt(mol)

    # H-bond donors and acceptors
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    num_HBD = rdMolDescriptors.CalcNumHBD(mol)
    num_HBA = rdMolDescriptors.CalcNumHBA(mol)

    return num_heavy_atoms, mol_weight, num_HBD, num_HBA

results = []

for fname in os.listdir(sdf_folder):
    if not fname.startswith("007-JAK1") or not fname.endswith(".sdf"):
        continue

    sdf_path = os.path.join(sdf_folder, fname)
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = suppl[0]

    if mol is None:
        continue

    ligand_id = fname.split("__")[-1].replace(".sdf", "")
    props = analyze_molecule(mol)

    results.append({
        "ligand_pdb": ligand_id,
        "num_heavy_atoms": props[0],
        "mol_weight": props[1],
        "num_HBD": props[2],
        "num_HBA": props[3],
    })

df_props = pd.DataFrame(results)
df_merged = df.merge(df_props, on="ligand_pdb", how="left")
df_merged.to_csv(csv_path, index=False)
