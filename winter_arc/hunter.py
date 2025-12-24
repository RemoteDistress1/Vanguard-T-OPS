import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
from rdkit import DataStructs

print(">> SYSTEM ONLINE. LOADING LIBRARY...")

# 1. Define the Enemy (Acridine)
target_smiles = "C1=CC=C2C(=C1)C=C3C=CC=CC3=N2"
target_mol = Chem.MolFromSmiles(target_smiles)
target_fp = Chem.RDKFingerprint(target_mol)
print(f">> TARGET LOCKED: ACRIDINE (MW: {Descriptors.MolWt(target_mol):.1f})")

# 2. Load the Haystack (Delaney CSV)
# The file has a header: "Compound ID", "ESOL predicted log solubility...", "SMILES", etc.
# We just need the "SMILES" column.
df = pd.read_csv('delaney.csv')
print(f">> DATABASE LOADED: {len(df)} CANDIDATES")

print("-" * 65)
print(f"{'RANK':<4} | {'NAME':<15} | {'SIM%':<5} | {'QED':<5} | {'STATUS'}")
print("-" * 65)

# 3. The Screening Loop
results = []
for index, row in df.iterrows():
    smiles = row['SMILES']
    mol = Chem.MolFromSmiles(smiles)
    if not mol: continue

    # Comparison Metrics
    mol_fp = Chem.RDKFingerprint(mol)
    similarity = DataStructs.TanimotoSimilarity(target_fp, mol_fp)
    qed_score = QED.qed(mol)
    
    # Filter: > 35% Similar
    if similarity > 0.35:
        # Use the SMILES as the name since this dataset doesn't have common names
        short_name = smiles[:10] + "..." 
        results.append({
            'Name': short_name,
            'SMILES': smiles,
            'Similarity': similarity,
            'QED': qed_score
        })

# 4. Sort and Report
results.sort(key=lambda x: x['Similarity'], reverse=True)

for i, res in enumerate(results[:15]):
    sim_pct = res['Similarity'] * 100
    print(f"{i+1:<4} | {res['Name']:<15} | {sim_pct:<5.1f}% | {res['QED']:.2f} | 🟢 MATCH")

# 5. Export
pd.DataFrame(results).to_csv('top_candidates.csv', index=False)
print("-" * 65)
print(f">> SCAN COMPLETE. {len(results)} CANDIDATES FOUND.")
