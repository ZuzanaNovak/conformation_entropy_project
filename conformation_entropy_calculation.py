import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd

ligand_folder = "ligands"
rotatable_bonds = []
not_found = []

#rotatable bonds calculation
for filename in os.listdir(ligand_folder):
    if filename.endswith(".sdf"):
        filepath = os.path.join(ligand_folder, filename)
        suppl = Chem.SDMolSupplier(filepath)
        for mol in suppl:
            if mol is not None:
                bond_number = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True)
                #target protein and ligand pdb separation
                name = filename.replace(".sdf", "")
                if "__" in name:
                    target_protein, ligand_pdb = name.split("__", 1)
                else:
                    target_protein = name
                    ligand_pdb = "none"

                rotatable_bonds.append({
                    "target_protein": target_protein,
                    "ligand_pdb": ligand_pdb,
                    "rotatable_bonds": bond_number
                })
            else:
                print(f"{filename} error with this molecule")
                not_found.append(filename)
    else:
        print(f"{filename} is not SDF")

#saving results
rotatable_bonds_df = pd.DataFrame(rotatable_bonds)
rotatable_bonds_df = rotatable_bonds_df.sort_values(by=["target_protein", "ligand_pdb"])
rotatable_bonds_df = rotatable_bonds_df.reset_index(drop=True)
print(rotatable_bonds_df)
rotatable_bonds_df.to_csv("results/ligand_rotatable_bonds.csv", index=False)
rotatable_bonds_df.to_excel("results/ligand_rotatable_bonds.xlsx", index=False)


#min/max rotatable bonds for every protein
min_max_df = rotatable_bonds_df.groupby("target_protein")["rotatable_bonds"].agg(["min", "max"]).reset_index()
min_max_df.columns = ["target_protein", "min_rotatable_bonds", "max_rotatable_bonds"]
min_max_df.to_csv("results/min_max_rotatable_bonds_per_protein.csv", index=False)
print(min_max_df)

#plotting
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.hist(rotatable_bonds_df["rotatable_bonds"], bins=range(rotatable_bonds_df["rotatable_bonds"].max()+2), edgecolor='black', align='left')
plt.title("Rotatable Bond Distribution")
plt.xlabel("Number of Rotatable Bonds")
plt.ylabel("Count")
plt.xticks(range(rotatable_bonds_df["rotatable_bonds"].max()+1))
plt.grid(axis='y', linestyle='-', alpha=0.7)
plt.tight_layout()
plt.savefig("results/rotatable_bonds_distribution_hist.png", dpi =300)
plt.show()