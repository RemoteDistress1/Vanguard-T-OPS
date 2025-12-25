
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from rdkit import Chem
from rdkit.Chem import Descriptors
import os

def train_qsar_model(library_path):
    print(f"🧠 AI ENGINE: Loading data from {library_path}...")
    
    # 1. Load Data
    suppl = Chem.SDMolSupplier(library_path)
    mols = [m for m in suppl if m is not None]
    
    if not mols:
        print("❌ Error: No molecules found.")
        return None

    # 2. Featurize (Turn molecules into numbers)
    print(f"⚙️  Processing {len(mols)} molecules into features...")
    data = []
    for mol in mols:
        # Calculate descriptors (The 'Inputs' for the AI)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Calculate Target (What we want to predict: QED Drug Likeness)
        qed = Descriptors.qed(mol) 
        
        data.append([mw, logp, hbd, hba, tpsa, qed])

    df = pd.DataFrame(data, columns=['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'QED_Score'])

    # 3. Train Model
    print("🎓 Training Random Forest Regressor...")
    X = df[['MW', 'LogP', 'HBD', 'HBA', 'TPSA']]
    y = df['QED_Score']
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # 4. Evaluate
    score = model.score(X_test, y_test)
    print(f"✅ Model Trained! Accuracy (R² Score): {score:.4f}")
    print("   (1.0 = Perfect Prediction, 0.0 = Random Guessing)")
    
    return model

if __name__ == "__main__":
    # Auto-run if executed directly
    path = os.path.join("Data_Depot", "Compound_Library.sdf")
    if os.path.exists(path):
        train_qsar_model(path)
    else:
        print("❌ Library not found in Data_Depot.")
