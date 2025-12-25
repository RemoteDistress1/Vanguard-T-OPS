# FILE: ORGANIC_I_TOOLS.py
# CONTEXT: Undergraduate Organic Chemistry (CHEM 2423)
# DESCRIPTION: Computational aids for homework and concept verification.

from rdkit import Chem
from rdkit.Chem import AllChem

# ==============================================================================
# CHAPTER 1: BONDING & MOLECULAR STRUCTURE
# ==============================================================================
def calculate_dou(smiles, name="Molecule"):
    '''
    Calculates Degree of Unsaturation (DoU).
    Formula: DoU = C - (H/2) - (X/2) + (N/2) + 1
    Usage: Quick check for rings/pi-bonds in unknown structures.
    '''
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return -1
    
    # Add implicit Hydrogens to get accurate count
    mol = Chem.AddHs(mol)
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    
    C = atoms.count('C')
    H = atoms.count('H')
    N = atoms.count('N')
    # Halogens count as 1, Oxygen/Sulfur count as 0
    X = sum(atoms.count(hal) for hal in ['F', 'Cl', 'Br', 'I'])
    
    dou = C - (H/2) - (X/2) + (N/2) + 1
    
    print(f"--- {name} DoU Analysis ---")
    print(f"Formula inputs: C={C}, H={H}, N={N}, X={X}")
    print(f"Degree of Unsaturation: {dou}")
    print("-" * 30)
    return dou

# ==============================================================================
# CHAPTER 2: RESONANCE & STABILITY
# ==============================================================================
def check_conjugation(smiles):
    '''
    Estimates resonance capability by counting conjugated bonds.
    Useful for comparing stability between two isomers.
    '''
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return 0
    
    conjugated_bonds = 0
    total_bonds = mol.GetNumBonds()
    
    for bond in mol.GetBonds():
        if bond.GetIsConjugated():
            conjugated_bonds += 1
            
    percentage = (conjugated_bonds / total_bonds) * 100 if total_bonds > 0 else 0
    
    print(f"--- Resonance Report ---")
    print(f"Conjugated Bonds: {conjugated_bonds} / {total_bonds}")
    print(f"Conjugation Coverage: {percentage:.1f}%")
    return percentage

# ==============================================================================
# CHAPTER 5: STEREOCHEMISTRY
# ==============================================================================
def check_chirality(smiles):
    '''
    Identifies Chiral Centers and assigns R/S configuration.
    Note: Requires 3D embedding for some complex cases, but works for most 2D.
    '''
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return
    
    try:
        # Force RDKit to calculate stereo
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        centers = Chem.FindMolChiralCenters(mol, includeStereo=True)
        
        print(f"--- Stereochemistry Report ---")
        if not centers: 
            print("Result: Achiral (Meso or no centers)")
        else:
            for idx, stereo in centers:
                atom_symbol = mol.GetAtomWithIdx(idx).GetSymbol()
                print(f"Atom {idx} ({atom_symbol}): {stereo}")
    except Exception as e:
        print(f"Error determining chirality: {e}")

# ==============================================================================
# FUTURE MODULES (Placeholder)
# ==============================================================================
# - Alkenes/Alkynes Reactions
# - Substitution (Sn1/Sn2) Logic
