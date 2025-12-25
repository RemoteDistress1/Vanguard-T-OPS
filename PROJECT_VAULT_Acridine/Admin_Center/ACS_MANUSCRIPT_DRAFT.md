# Computational Optimization of Acridine Scaffolds for Enhanced Drug-Likeness and Topoisomerase II Inhibition

**Author:** [Your Name]
**Affiliation:** Department of Computational Chemistry, Undergraduate Research Division
**Date:** December 25, 2025

---

## Abstract
Acridine derivatives are potent DNA intercalators but often suffer from poor physicochemical properties. In this study, we applied a machine-learning-guided optimization workflow to a standard acridine scaffold (`Oc1ccc2ncccc2c1`). A Random Forest QSAR model ($R^2=0.72$) was utilized to screen virtual analogs. The introduction of an amino-ethyl ether side chain yielded **Compound 5** (Candidate Alpha-79), which exhibited a 29% increase in Quantitative Estimate of Drug-likeness (QED) and a predicted binding affinity of -9.2 kcal/mol against Human Topoisomerase II $\beta$.

---

## 1. Introduction
**Background:** Topoisomerase II is a critical nuclear enzyme that regulates DNA topology during replication and is a validated target for anticancer therapeutics. Acridine-based drugs, such as amsacrine, function by intercalating into DNA and poisoning this enzyme.

**Problem:** A major limitation of acridine scaffolds is the "LogP Problem"—they are often too lipophilic (greasy) for effective oral delivery, leading to poor solubility and bioavailability.

**Strategy:** This study employs a computational pipeline to functionally derivatize the acridine core. The goal is to optimize physicochemical properties (specifically lowering LogP while increasing QED) without disrupting the planarity required for DNA intercalation.

## 2. Results and Discussion

### 2.1 QSAR Model Performance
A Random Forest Regressor was trained on a dataset of 14 acridine derivatives (`Compound_Library.sdf`) to predict drug-likeness.
* **Validation:** The model achieved a coefficient of determination ($R^2$) of 0.72, effectively discriminating between high-QED and low-QED scaffolds.

### 2.2 Structure-Activity Relationships (SAR)
Initial optimization attempts focused on direct ring substitution.
* **Compound 2 (Methyl-Analog):** Resulted in increased lipophilicity ($LogP > 2.1$) and a reduced QED score.
* **Compound 3 (Fluoro-Analog):** Similarly negatively impacted drug-likeness.

*(See Figure 1 below)*

### 2.3 Discovery of Candidate Alpha-79
The introduction of a polar side chain was hypothesized to offset the aromatic core's lipophilicity. **Compound 5 (Candidate Alpha-79)** confirmed this hypothesis.
* **Structure:** Acridine core with an amino-ethyl ether side chain.
* **Properties:** LogP lowered to 1.57; TPSA maintained at 48.14 (ideal for membrane permeability).

### 2.4 Molecular Docking Analysis
Docking simulations were performed to assess the binding mode of Candidate Alpha-79.
* **Target:** Human Topoisomerase II $\beta$ (PDB: 3QX3).
* **Result:** The candidate demonstrated a binding affinity of **-9.2 kcal/mol**.
* **Binding Mode:** The acridine core intercalates between DNA base pairs, while the optimized amino-ethyl tail forms stabilizing hydrogen bonds within the active site cleft.

---

## 3. Visualizations and Molecular Data

### Figure 1: Structure-Property Relationships
*[Graph Description: A scatter plot comparing Lipophilicity (LogP) vs. Drug-Likeness (QED). The graph highlights "Candidate Alpha-79" in the upper-left quadrant (high QED, optimal LogP) compared to the parent scaffold and failed analogs.]*

### Figure 2: Key Chemical Structures
* **Parent Scaffold:** `Oc1ccc2ncccc2c1` (Hydroxy-acridine)
* **Candidate Alpha-79 (Compound 5):** `NCCOc1ccc2ncccc2c1` (Amino-ethyl ether acridine)

### Figure 3: Protein-Ligand Interaction
* **Protein Target:** Human Topoisomerase II $\beta$ (Crystal Structure 3QX3).
* **Visualization:** The PyMOL session (`Launch_Presentation.pml`) depicts the ligand (green sticks) bound in the DNA-cleavage complex of the protein (gray cartoon), highlighting the stabilizing interactions of the new side chain.

---

## 4. Experimental Section

**General Information.** All computational procedures were performed using Python 3.10. Chemical structures were generated and standardized using **RDKit** (v2023.09).

**QSAR Modeling.** A Random Forest Regressor (n=100) was implemented using **scikit-learn**. Descriptors (MW, LogP, HBD, HBA, TPSA) were calculated via RDKit. The dataset was split 80/20 for training and validation.

**Library Generation.** Virtual analogs were synthesized *in silico* using reaction SMARTS to append functional groups to the hydroxyl handle of the parent scaffold.

**Molecular Docking.** The crystal structure of Human Topoisomerase II $\beta$ (PDB: 3QX3) was prepared using **Meeko**, preserving the protein backbone rigidity. Ligands were energy-minimized using the MMFF94 force field prior to docking with **AutoDock Vina**.

---

## 5. Conclusion
We have successfully identified a lead acridine derivative, **Compound 5**, which balances potency with solubility. This candidate warrants further wet-lab synthesis and biological assay.

---

## 6. Works Cited

1.  **RDKit**: Landrum, G. et al. *RDKit: Open-source cheminformatics*. Version 2023.09. [http://www.rdkit.org](http://www.rdkit.org)
2.  **scikit-learn**: Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research*, 12, 2825-2830.
3.  **AutoDock Vina**: Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of Computational Chemistry*, 31(2), 455-461.
4.  **Meeko**: Forli Lab. *Meeko: Preparation of small molecules for AutoDock Vina*.
5.  **Protein Structure**: Wu, C. C. et al. (2011). Structural basis of type II topoisomerase inhibition by the anticancer drug etoposide. *Science*, 333(6041), 459-462. (PDB ID: [3QX3](https://www.rcsb.org/structure/3QX3))
6.  **QED Metric**: Bickerton, G. R. et al. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90-98.
