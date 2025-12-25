
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from AI_DISCOVERY import train_qsar_model
import os

# 1. Reload the Brain
library_path = os.path.join("Data_Depot", "Compound_Library.sdf")
model = train_qsar_model(library_path)

if model:
    print("\n🧪 AI OPTIMIZER v2: Bioisosteric Scanning...")
    
    # 2. Load the Best Lead
    suppl = Chem.SDMolSupplier(library_path)
    mols = [m for m in suppl if m is not None]
    best_mol = max(mols, key=lambda m: Descriptors.qed(m))
    best_score = Descriptors.qed(best_mol)
    
    # Get the SMART pattern for an aromatic Carbon (c)
    aromatic_carbon = Chem.MolFromSmarts('c')
    
    print(f"   >> Base Structure QED: {best_score:.4f}")
    print(f"   >> Strategy: Swapping Ring Carbons (C) -> Nitrogen (N)")

    # 3. The Swap Operation
    matches = best_mol.GetSubstructMatches(aromatic_carbon)
    variations = []
    
    print(f"   >> Found {len(matches)} possible substitution sites...")

    for i, match in enumerate(matches):
        # Create a writable version of the molecule
        rw_mol = Chem.RWMol(best_mol)
        idx = match[0]
        
        # Change atom at index 'idx' to Nitrogen (Atomic Num = 7)
        rw_mol.GetAtomWithIdx(idx).SetAtomicNum(7)
        
        # Sanitize to check if the chemistry is valid
        try:
            Chem.SanitizeMol(rw_mol)
            new_mol = rw_mol.GetMol()
            
            # Predict Score
            mw = Descriptors.MolWt(new_mol)
            logp = Descriptors.MolLogP(new_mol)
            hbd = Descriptors.NumHDonors(new_mol)
            hba = Descriptors.NumHAcceptors(new_mol)
            tpsa = Descriptors.TPSA(new_mol)
            
            features_df = pd.DataFrame([[mw, logp, hbd, hba, tpsa]], 
                                     columns=['MW', 'LogP', 'HBD', 'HBA', 'TPSA'])
            
            pred_score = model.predict(features_df)[0]
            variations.append((f"Pos_{idx}_Nitrogen_Swap", pred_score))
            
        except:
            continue

    # 4. Results
    print("\n📊 PREDICTION RESULTS (Higher is Better):")
    print("-" * 45)
    
    sorted_vars = sorted(variations, key=lambda x: x[1], reverse=True)
    
    for name, score in sorted_vars[:5]:
        diff = score - best_score
        symbol = "⬆️" if diff > 0 else "⬇️"
        print(f"   {name}: {score:.4f} ({symbol} {diff:+.4f})")
        
    if not sorted_vars:
        print("   No valid substitutions found.")
    print("-" * 45)
