
import csv
import sys
# SAFETY CHECK: Ensure we can find our tools
sys.path.append('.') 

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors # <--- THE CRITICAL FIX
    import chem_tools
except ImportError:
    print("!!! CRITICAL: RDKit not found in this environment.")
    sys.exit(1)

print(f"{'ID':<6} | {'NAME':<15} | {'DoU':<5} | {'LIPINSKI':<10} | {'STATUS'}")
print("-" * 65)

try:
    with open('targets.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mol = Chem.MolFromSmiles(row['SMILES'])
            if not mol:
                print(f"{row['ID']:<6} | {row['Name']:<15} | ERR   | INVALID    | ❌ DEAD")
                continue

            # 1. CALCULATE DATA (Using the fixed import)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            passed = (mw < 500) and (logp < 5)
            
            # 2. CALCULATE DoU (Using our library)
            # We suppress the library's print output for a clean table
            original_stdout = sys.stdout
            sys.stdout = None # Mute output
            dou = chem_tools.calculate_dou(row['SMILES'])
            sys.stdout = original_stdout # Unmute

            # 3. PRINT ROW
            status = "✅ DEPLOY" if passed else "⚠️ REJECT"
            print(f"{row['ID']:<6} | {row['Name']:<15} | {dou:<5.1f} | {'PASS' if passed else 'FAIL':<10} | {status}")

except FileNotFoundError:
    print("!!! ERROR: targets.csv not found. Run the Ammo Fabricator.")
    
print("-" * 65)
print(">>> BATCH FIRE COMPLETE.")
