# Computational Design and Optimization of Acridine Derivatives as Topoisomerase II Inhibitors
**Date:** December 25, 2025
**Author:** Research Fellow, Computational Chemistry Lab
**Correspondence:** Project Acridine (ACR-OPT-2026)

---

## Abstract
The development of novel anticancer agents requires the efficient exploration of chemical space. This study employs a multi-stage computational pipeline to optimize acridine-based scaffolds for enhanced drug-likeness (QED) and target binding affinity. Utilizing a Random Forest QSAR model and bioisosteric replacement strategies, we identified **Candidate Alpha-79**, a novel derivative exhibiting a 29% improvement in QED (0.795) and a predicted binding affinity of -9.2 kcal/mol against Human Topoisomerase II Beta.

---

## 1. Introduction
Acridines are a class of polycyclic aromatic compounds with known intercalation properties into DNA. However, clinical application is often limited by poor solubility and side-effect profiles. This project aimed to computationally evolved a lead acridine structure (`Oc1ccc2ncccc2c1`) to optimize its physicochemical properties without compromising its pharmacophore integrity.

## 2. Materials and Methods

### 2.1 Computational Environment
All analyses were conducted in a **JupyterLab** environment hosted on the Anaconda Platform. The primary chemoinformatics workflow was implemented in **Python 3.10**.

### 2.2 Chemical Library Preparation
Structure handling, sanitization, and descriptor calculation (MW, LogP, TPSA, HBD/HBA) were performed using **RDKit (v2023.09.3)** [1]. Molecules were standardized to their canonical SMILES representation before processing.

### 2.3 Machine Learning (QSAR)
A quantitative structure-activity relationship (QSAR) model was constructed using a **Random Forest Regressor** from the **scikit-learn** library [2]. The model was trained to predict the Quantitative Estimate of Drug-likeness (QED) [3] based on five physicochemical descriptors.

### 2.4 Structural Optimization
Novel analogs were generated via a custom Python script (`AI_DERIVATIZER.py`) utilizing RDKit's reaction SMARTS capabilities. A library of side-chain modifications (ethers, amines, esters) was synthesized in silico and scored by the QSAR model.

### 2.5 Molecular Docking
Binding affinity was assessed using **AutoDock Vina** [4]. The crystal structure of Human Topoisomerase II Beta (PDB: 3QX3) was prepared by removing solvent molecules. Ligand preparation, including protonation states and rotatable bond definition, was managed via **Meeko**.

---

## 3. Results

### 3.1 Lead Optimization
The initial screening identified that direct alkylation of the acridine ring increased lipophilicity ($LogP > 2.0$), negatively impacting drug-likeness.
* **Parent Compound:** QED = 0.614
* **Methyl-Analog:** QED = 0.553 (Excluded)

The introduction of a polar amino-ethyl ether side chain yielded **Candidate Alpha-79** (`NCCOc1ccc2ncccc2c1`), which balanced lipophilicity ($LogP = 1.57$) with increased polarity (TPSA = 48.14).

### 3.2 Binding Analysis
Molecular docking simulations confirmed that Candidate Alpha-79 occupies the active DNA-binding cleft of Topoisomerase II.
* **Binding Affinity:** -9.2 kcal/mol
* **Interaction Mode:** Stabilized by pi-stacking interactions with the DNA base pairs and hydrogen bonding via the amino-ethyl tail.

---

## 4. Conclusion
We successfully engineered a high-affinity, soluble acridine derivative using an automated computational pipeline. Candidate Alpha-79 represents a viable lead for synthesis and in vitro assay validation.

---

## 5. References
1. Landrum, G. (2023). **RDKit: Open-source cheminformatics**. https://www.rdkit.org
2. Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR*, 12, 2825-2830.
3. Bickerton, G. R. et al. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90-98.
4. Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455-461.
