import os
import re
import pandas as pd
# File paths
input_csv = "rotatable_bonds/results/ligand_rotatable_bonds.csv"
filtered_csv = "CREST/results/tables/007_JAK1_rotatable_bonds.csv"
crest_directory = "CREST/CREST_on_007-JAK1"

with open(input_csv, "r") as infile:
    lines = infile.readlines()

header, *data_lines = lines
matching_lines = [line for line in data_lines if line.startswith("007-JAK1")]

if matching_lines:
    with open(filtered_csv, "w") as outfile:
        outfile.write(header)
        outfile.writelines(matching_lines)

df = pd.read_csv(filtered_csv)

def extract_sconf(ligand_id):
    log_file = os.path.join(crest_directory, f"{ligand_id}_crest.log")
    
    if os.path.exists(log_file):
        with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            match = re.search(r"Sconf\s*=\s*([0-9.]+)", content)
            if match:
                sconf = float(match.group(1))
                print(f"Sconf for {ligand_id}: {sconf}")
                return sconf
            else:
                print(f"No Sconf found in {log_file}")
    else:
        print(f"File not found: {log_file}")
    
    return None


df["Sconf"] = df["ligand_pdb"].apply(extract_sconf)


df.to_csv(filtered_csv, index=False)
