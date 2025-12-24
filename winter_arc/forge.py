import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

print(">> IRON FORGE: ONLINE")
print(">> INITIATING 3D PROTOTYPING FOR ELITE CANDIDATES...")

# 1. Load the Elite Squad
df = pd.read_csv('top_candidates.csv')
print(f">> PROCESSING {len(df)} TARGETS...")

# 2. Prepare the Output Container
writer = Chem.SDWriter('elite_vanguard.sdf')

success_count = 0

print("-" * 50)
print(f"{'RANK':<4} | {'NAME':<15} | {'STATUS'}")
print("-" * 50)

for index, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    name = f"Rank_{index+1}_{row['Name']}"
    
    if mol:
        mol.SetProp("_Name", name)
        
        # STEP A: Add Hydrogens (Critical for mass)
        mol_3d = Chem.AddHs(mol)
        
        # STEP B: 3D Embedding (The Inflation)
        res = AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        
        if res == 0:
            # STEP C: Physics Optimization (MMFF Force Field)
            AllChem.MMFFOptimizeMolecule(mol_3d)
            
            # STEP D: Secure the Asset
            writer.write(mol_3d)
            success_count += 1
            print(f"{index+1:<4} | {name[:15]:<15} | 🧊 3D FORGED")
        else:
            print(f"{index+1:<4} | {name[:15]:<15} | ⚠️ STRUCTURE UNSTABLE")

writer.close()
print("-" * 50)
print(f">> FORGE COMPLETE. {success_count} PROTOTYPES BUILT.")
print(">> ASSETS SECURED: 'elite_vanguard.sdf'")
