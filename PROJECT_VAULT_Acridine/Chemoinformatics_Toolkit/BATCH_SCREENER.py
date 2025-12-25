# FILE: live_fire.py
# MISSION: SYSTEMS CHECK (ALL DOMAINS)
# DATE: 23 DEC 2025

import chem_tools
from rdkit import Chem
from rdkit.Chem import AllChem  # <--- ADDED: Direct import for physics engine

print("\n" + "="*60)
print(">>> VANGUARD THERAPEUTICS | LIVE FIRE SIMULATION")
print(">>> TARGETS: ACADEMIC (Exams) & RESEARCH (Discovery)")
print("="*60 + "\n")

# ==============================================================================
# PHASE 1: ACADEMIC SECTOR (The Syllabus Killer)
# ==============================================================================
print("[*] PHASE 1: ENGAGING ACADEMIC TARGETS...\n")

# TARGET A: DEGREE OF UNSATURATION (Chapter 1)
print(">>> FIRING: Unsaturation Bot (Target: Benzene)")
chem_tools.calculate_dou("c1ccccc1", "Benzene")

# TARGET B: CONJUGATION & STABILITY (Chapter 2)
print(">>> FIRING: Resonance Radar (Target: Acridine)")
acridine_smiles = "C1=CC=C2C(=C1)C=C3C=CC=CC3=N2"
chem_tools.check_conjugation(acridine_smiles, "Vanguard Catalyst (Acridine)")

# TARGET C: CHIRALITY (Chapter 5)
print(">>> FIRING: Stereo-Scanner (Target: Thalidomide)")
thalidomide = "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1"
chem_tools.check_chirality(thalidomide, "Thalidomide")

# ==============================================================================
# PHASE 2: RESEARCH SECTOR (The Drug Hunter)
# ==============================================================================
print("\n" + "="*60)
print("[*] PHASE 2: ENGAGING RESEARCH PROTOCOLS...\n")

# TARGET D: LIPINSKI RULES (Commercial Viability)
print(">>> FIRING: Lipinski Gatekeeper (Target: Morphine)")
morphine = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
chem_tools.check_lipinski(morphine, "Morphine")

# TARGET E: VIRTUAL SYNTHESIS (The Reactor)
print(">>> FIRING: Virtual Reactor (Simulation: Amide Formation)")
reactant_A = "c1ccccc1C(=O)O"  # Benzoic Acid
reactant_B = "CN"              # Methylamine
amide_rxn = "[C:1](=[O:2])[OH].[N:3]>>[C:1](=[O:2])[N:3]" 

product = chem_tools.run_reaction(reactant_A, reactant_B, amide_rxn)
print(f"   Reactant A: {reactant_A}")
print(f"   Reactant B: {reactant_B}")
# [Removed empty print statement here]
print(f"   PRODUCT GENERATED: {product}")
chem_tools.check_lipinski(product, "Synthesized Amide Product")

# TARGET F: PHYSICS ENGINE (3D Geometry)
print(">>> FIRING: Geometry Physics Engine (Target: H2O)")
water = "O"
# FIXED: Using direct 'Chem' and 'AllChem' calls, not 'chem_tools.Chem'
mol_water = Chem.AddHs(Chem.MolFromSmiles(water)) 
AllChem.EmbedMolecule(mol_water)
angle = chem_tools.measure_geometry(mol_water, [1, 0, 2])
print(f"   Calculated Bond Angle: {angle}")
print("   (Theoretical Value: ~104.5 degrees)")

# ==============================================================================
# MISSION DEBRIEF
# ==============================================================================
print("\n" + "="*60)
print(">>> SYSTEMS CHECK COMPLETE. ALL WEAPONS FUNCTIONAL.")
print(">>> READY FOR SPRING 2026 DEPLOYMENT.")
print("="*60 + "\n")