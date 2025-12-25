# CHEM_TOOLS V6.0 REPAIR
import sys
try:
    from rdkit import Chem, AllChem, Descriptors
    import py3Dmol
except ImportError:
    pass # Handled in main code

def calculate_dou(s, n="M"):
    m = Chem.MolFromSmiles(s)
    if not m: return -1
    m = Chem.AddHs(m)
    atoms = [a.GetSymbol() for a in m.GetAtoms()]
    dou = atoms.count('C') - (atoms.count('H')/2) - (sum(atoms.count(x) for x in 'FClBrI')/2) + (atoms.count('N')/2) + 1
    print(f"DoU: {dou}")
    return dou

def check_lipinski(s, n="M"):
    m = Chem.MolFromSmiles(s)
    if not m: return {"passed": False}
    mw = Descriptors.MolWt(m)
    logp = Descriptors.MolLogP(m)
    passed = mw < 500 and logp < 5
    print(f"MW: {mw:.1f} | LogP: {logp:.1f} | {'PASS' if passed else 'FAIL'}")
    return {"passed": passed}

# ... (Simplified for verification) ...
