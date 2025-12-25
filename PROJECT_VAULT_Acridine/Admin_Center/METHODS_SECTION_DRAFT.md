
## 2. Materials and Methods (Expanded)

### 2.1 Computational Libraries & Environment
All chemoinformatic analyses were performed in a Python 3.10 environment. 
Structure handling and descriptor calculation were conducted using **RDKit** (v2023.09.3) [1]. 
Machine learning models were implemented using **scikit-learn** (v1.3.0) [2].

### 2.2 Machine Learning Pipeline
A Random Forest Regressor (n_estimators=100, random_state=42) was trained on a curated library of acridine derivatives. 
The model utilized five key physicochemical descriptors: Molecular Weight (MW), LogP (Lipophilicity), 
Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA), and Topological Polar Surface Area (TPSA). 
Model performance was evaluated via the Coefficient of Determination ($R^2$), achieving a score of 0.72 on the validation set.

### 2.3 Structure Optimization & Generation
Candidate structures were generated via a custom bioisosteric replacement algorithm and side-chain derivatization script (`AI_DERIVATIZER.py`). 
Lipinski's Rule of Five was applied as a post-generation filter to ensure oral bioavailability.

### 2.4 Molecular Docking
Molecular docking studies were performed using **AutoDock Vina** (v1.2.3) [3]. 
The crystal structure of Human Topoisomerase II Beta (PDB ID: 3QX3) was retrieved from the RCSB Protein Data Bank.
Ligand preparation (PDBQT conversion) was handled via **Meeko** [4]. 
A search space of 20x20x20 Å was defined around the active site (Center: 10, 40, 20). 
Binding affinity was calculated in kcal/mol, with an exhaustiveness setting of 8.

### 2.5 Software Availability
The custom code used for this study has been consolidated into the `DrugDiscovery_Toolkit` Python library 
and is available in the supplementary material to ensure reproducibility.
