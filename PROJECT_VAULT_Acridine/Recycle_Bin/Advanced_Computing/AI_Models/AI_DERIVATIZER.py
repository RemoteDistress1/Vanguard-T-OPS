
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from AI_DISCOVERY import train_qsar_model
import os

# 1. Reload the Brain
library_path = os.path.join("Data_Depot", "Compound_Library.sdf")
model = train_qsar_model(library_path)

if model:
    print("\n🧪 AI DERIVATIZER: Targeted Side-Chain Optimization...")
    
    # 2. Load the Best Lead
    suppl = Chem.SDMolSupplier(library_path)
    mols = [m for m in suppl if m is not None]
    best_mol = max(mols, key=lambda m: Descriptors.qed(m))
    best_score = Descriptors.qed(best_mol)
    
    print(f"   >> Base Structure QED: {best_score:.4f}")
    
    # 3. Define the 'Reaction' (Alkylation of Phenol/Alcohol)
    # This SMARTS pattern says: "Take an OH group and attach a Carbon chain to it"
    rxn_smarts = '[O:1][H]>>[O:1][C]' 
    # Note: For this demo, we will use a simpler approach: Re-building from SMILES
    # to ensure we don't get 'Kekulize' errors.
    
    base_smiles = Chem.MolToSmiles(best_mol)
    print(f"   >> Base SMILES: {base_smiles}")
    
    # List of groups to attach to the Oxygen (assuming it starts with 'O')
    # We strip the first 'O' and replace it with 'Group-O'
    # Or more safely: we look for the alcohol handle.
    
    tails = {
        'Methoxy': 'CO',          # -OCH3
        'Ethoxy': 'CCO',          # -OCH2CH3
        'Isopropoxy': 'CC(C)O',   # -OCH(CH3)2
        'Acetyl': 'CC(=O)O',      # -OCOCH3 (Ester)
        'Trifluoroethoxy': 'FC(F)(F)CO', # -OCH2CF3 (Fluorinated)
        'Amino-Ethyl': 'NCCO',    # -OCH2CH2NH2 (Solubility boost)
    }

    print(f"   >> Synthesizing {len(tails)} derivatives...")
    
    variations = []
    
    # Identify the core (Remove the OH or O- group to get the scaffold)
    # This is a heuristic for the demo.
    if base_smiles.startswith('O'):
        scaffold = base_smiles[1:] # Chop off the first O
    elif base_smiles.startswith('Oc'):
        scaffold = base_smiles[2:]
        tails = {k: v + 'c' for k, v in tails.items()} # Fix connection
    else:
        # Fallback: Just append to end (less precise)
        scaffold = base_smiles
    
    for name, tail_smiles in tails.items():
        try:
            # Stitch them together
            new_smiles = tail_smiles + scaffold
            new_mol = Chem.MolFromSmiles(new_smiles)
            
            if new_mol:
                # Score it
                mw = Descriptors.MolWt(new_mol)
                logp = Descriptors.MolLogP(new_mol)
                hbd = Descriptors.NumHDonors(new_mol)
                hba = Descriptors.NumHAcceptors(new_mol)
                tpsa = Descriptors.TPSA(new_mol)
                
                features_df = pd.DataFrame([[mw, logp, hbd, hba, tpsa]], 
                                         columns=['MW', 'LogP', 'HBD', 'HBA', 'TPSA'])
                
                pred_score = model.predict(features_df)[0]
                variations.append((name, pred_score))
        except:
            continue

    # 4. Results
    print("\n📊 DERIVATIVE RANKING (Higher is Better):")
    print("-" * 45)
    
    sorted_vars = sorted(variations, key=lambda x: x[1], reverse=True)
    
    for name, score in sorted_vars:
        diff = score - best_score
        symbol = "⬆️" if diff > 0 else "⬇️"
        print(f"   {name}: {score:.4f} ({symbol} {diff:+.4f})")
        
    print("-" * 45)
