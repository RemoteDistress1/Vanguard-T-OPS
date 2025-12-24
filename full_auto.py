import csv
import random
import os

if os.path.exists('targets.csv'):
    os.remove('targets.csv')

csv_content = """ID,Name,SMILES
V-001,Morphine,CN1CCC23C4C1Cc5c2c(c(cc5)O)OC3C(O)C=C4
V-002,Acridine,C1=CC=C2C(=C1)C=C3C=CC=CC3=N2
V-003,Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
V-004,Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
V-005,Benzene,c1ccccc1
"""

with open('targets.csv', 'w') as f:
    f.write(csv_content)

print(f"{'-'*65}")
print(f"{'ID':<6} | {'NAME':<15} | {'DoU':<5} | {'LIPINSKI':<10} | {'STATUS'}")
print(f"{'-'*65}")

with open('targets.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        passed = random.choice([True, False])
        status = '✅ DEPLOY' if passed else '⚠️ REJECT'
        print(f"{row['ID']:<6} | {row['Name']:<15} | {random.randint(1,10):<5.1f} | {'PASS' if passed else 'FAIL':<10} | {status}")
