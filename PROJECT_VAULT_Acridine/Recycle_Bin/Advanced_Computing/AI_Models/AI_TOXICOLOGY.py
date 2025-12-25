
import os
from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem import FragmentMatcher
import sys

# Append RDKit Contrib path for SA Score (Synthetic Accessibility)
# Note: This often lives in a specific folder. For this demo, we will use a heuristic 
# because standard installs sometimes hide the 'SAScore' module.
# We will focus on "Structural Alerts" which are native.

def run_safety_check(smiles):
    print(f"🛡️  SAFETY CHECK: Analyzing {smiles}...")
    mol = Chem.MolFromSmiles(smiles)
    
    issues = []
    
    # 1. Check for 'PAINS' (Pan-Assay Interference Compounds)
    # These are chemical groups that react with *everything* (False Positives)
    # We define a few common "Bad Actors" manually for this demo.
    
    bad_patterns = {
        "Nitro Group (Explosive/Toxic Risk)": "[N+](=O)[O-]",
        "Peroxide (Unstable)": "OO",
        "Aldehyde (Reactive)": "[CH]=O",
        "Michael Acceptor (Covalent Binder)": "C=CC(=O)"
    }
    
    print("   >> Scanning for Toxicophores...")
    for name, smarts in bad_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            issues.append(f"ALERT: Contains {name}")
            
    # 2. Fragment Analysis (Is it too weird?)
    num_rotatable = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable > 10:
        issues.append("WARNING: Too flexible (>10 Rotatable Bonds). Poor oral bioavailablity.")
        
    num_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 5:
        issues.append("WARNING: Too rigid/complex (>5 Rings). Hard to synthesize.")

    # 3. Report
    print("\n📋 TOXICOLOGY REPORT:")
    print("-" * 30)
    if not issues:
        print("✅ CLEAN BILL OF HEALTH.")
        print("   No obvious structural alerts found.")
        print("   Molecule appears stable and non-reactive.")
    else:
        for alert in issues:
            print(f"❌ {alert}")
    print("-" * 30)

if __name__ == "__main__":
    # Test on our Champion
    champion_smiles = "NCCOc1ccc2ncccc2c1"
    run_safety_check(champion_smiles)
