
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from AI_DISCOVERY import train_qsar_model
import os

# 1. Reload the Brain
library_path = os.path.join("Data_Depot", "Compound_Library.sdf")
model = train_qsar_model(library_path)

if model:
    print("\n🧪 AI OPTIMIZER: Starting Virtual Screen...")
    
    # 2. Load the Best Lead
    suppl = Chem.SDMolSupplier(library_path)
    mols = [m for m in suppl if m is not None]
    best_mol = max(mols, key=lambda m: Descriptors.qed(m))
    best_score = Descriptors.qed(best_mol)
    print(f"   >> Current Best Candidate: QED {best_score:.4f}")

    # 3. Virtual Mutations
    smiles = Chem.MolToSmiles(best_mol)
    variations = []
    
    # Let's try different groups this time since the last ones failed
    modifications = {
        'Methylated': 'C', 
        'Fluorinated': 'F', 
        'Chlorinated': 'Cl', 
        'Amine': 'N',
        'Methoxy': 'OC'
    }

    print(f"   >> Generating {len(modifications)} virtual analogs...")
    
    for name, group in modifications.items():
        new_smiles = smiles + group 
        new_mol = Chem.MolFromSmiles(new_smiles)
        
        if new_mol:
            # Calculate features
            mw = Descriptors.MolWt(new_mol)
            logp = Descriptors.MolLogP(new_mol)
            hbd = Descriptors.NumHDonors(new_mol)
            hba = Descriptors.NumHAcceptors(new_mol)
            tpsa = Descriptors.TPSA(new_mol)
            
            # --- THE FIX IS HERE ---
            # We create a DataFrame with columns to match the training data
            features_df = pd.DataFrame([[mw, logp, hbd, hba, tpsa]], 
                                     columns=['MW', 'LogP', 'HBD', 'HBA', 'TPSA'])
            
            predicted_score = model.predict(features_df)[0]
            variations.append((name, predicted_score))

    # 4. Results
    print("\n📊 PREDICTION RESULTS (Higher is Better):")
    print("-" * 40)
    for name, score in sorted(variations, key=lambda x: x[1], reverse=True):
        diff = score - best_score
        symbol = "⬆️" if diff > 0 else "⬇️"
        print(f"   {name} Analog: {score:.4f} ({symbol} {diff:+.4f})")
    print("-" * 40)
