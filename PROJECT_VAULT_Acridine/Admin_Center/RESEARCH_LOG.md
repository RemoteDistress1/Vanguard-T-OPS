# 🧪 Computational Discovery Log: Project Acridine
**Date:** 2025-12-25
**Principal Investigator:** Undergraduate Research Lead
**Subject:** Optimization of Acridine Scaffolds for Drug-Likeness (QED)

## 1. Executive Summary
We successfully identified and optimized a lead candidate, **Candidate_Alpha_79**, derived from the Acridine core. The new candidate exhibits a **29% increase in QED** (Drug-Likeness) and a significantly improved solubility profile compared to the starting material.

## 2. Methodology
* **Initial Screening:** Processed 'Compound_Library.sdf' (14 candidates).
* **AI Model:** Random Forest Regressor trained on physicochemical descriptors (MW, LogP, HBD, HBA, TPSA).
* **Optimization Strategy:**
    * *Attempt 1 (Direct Addition):* Failed. Fluorination/Methylation increased lipophilicity (LogP > 2.0).
    * *Attempt 2 (Bioisosteric Swap):* Failed. Nitrogen swaps disrupted aromaticity.
    * *Attempt 3 (Derivatization):* **Success.** Amino-ethyl ether linkage lowered LogP (1.94 -> 1.57) and boosted TPSA.

## 3. Key Findings (The "Table 1")
| Compound ID | Structure (SMILES) | QED Score | LogP (Lipophilicity) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Original Lead** | `Oc1ccc2ncccc2c1` | 0.6141 | 1.94 | Reference |
| *Methyl-Analog* | `COc1ccc2ncccc2c1` | 0.5534 | 2.12 | Rejected (Too Greasy) |
| *Fluoro-Analog* | `Fc1ccc2ncccc2c1` | 0.5606 | 2.35 | Rejected (Toxic Risk) |
| **Candidate Alpha-79** | `NCCOc1ccc2ncccc2c1` | **0.7949** | **1.57** | **Lead Candidate** |

## 4. Toxicology Profile
* **PAINS Filter:** Pass (No interference structures).
* **Structural Alerts:** None.
* **Synthesizability:** High (Standard Ether Synthesis from Acridin-ol).

## 5. Conclusion
Candidate_Alpha_79 is recommended for wet-lab synthesis.
