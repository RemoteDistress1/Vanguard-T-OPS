
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

print("🔍 AI DIAGNOSTICS: Deconstructing the QED Score...")

# 1. Define our contestants
# The Original vs. The 'Best' Failure (Amino-Ethyl) vs. A Different Failure (Acetyl)
contestants = {
    "Original Lead": "Oc1ccc2ncccc2c1",
    "Amino-Ethyl": "NCCOc1ccc2ncccc2c1",  
    "Acetyl (Ester)": "CC(=O)Oc1ccc2ncccc2c1"
}

data = []

for name, smiles in contestants.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # Calculate the "Big 4" Drug Properties
        mw = Descriptors.MolWt(mol)       # Weight (< 500 is good)
        logp = Descriptors.MolLogP(mol)   # Greasiness (< 5 is good)
        hbd = Descriptors.NumHDonors(mol) # H-Donors (< 5 is good)
        tpsa = Descriptors.TPSA(mol)      # Polarity (Related to cell permeability)
        qed = Descriptors.qed(mol)        # The Master Score
        
        data.append([name, qed, mw, logp, hbd, tpsa])

# 2. Create the Comparison Table
df = pd.DataFrame(data, columns=["Molecule", "QED Score", "Weight (MW)", "Lipophilicity (LogP)", "H-Donors", "Polarity (TPSA)"])

# 3. Display
print("\n📊 COMPARATIVE ANALYSIS:")
print(df.to_string(index=False))

print("\n🔎 ANALYSIS GUIDE:")
print("   - If MW > 500: Too heavy.")
print("   - If LogP > 5: Too greasy (won't dissolve).")
print("   - If LogP < 0: Too polar (won't cross membranes).")
