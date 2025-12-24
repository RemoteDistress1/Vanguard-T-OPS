import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

print(">> SATELLITE UPLINK ESTABLISHED...")

# 1. Load the Top Candidates
df = pd.read_csv('top_candidates.csv')

# Take the Top 9 (Rank 1 is Acridine, Ranks 2-9 are the new recruits)
top_9 = df.head(9)
print(f">> PROCESSING {len(top_9)} TARGETS FOR VISUALIZATION...")

# 2. Prepare the Molecules
mols = []
legends = []

for index, row in top_9.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        # Legend: Rank # + Similarity Score
        label = f"#{index+1} (Sim: {row['Similarity']*100:.1f}%)"
        legends.append(label)

# 3. Take the Photo
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(300, 300),
    legends=legends
)

img.save('top_candidates.png')
print(">> IMAGERY SECURED: 'top_candidates.png' saved to disk.")
