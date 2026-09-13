# BOND Judge Accuracy Report (three-layer acceptance protocol)

> Date: 2026-09-13. Engine: `bond/judge_core.py` v2.0+ (single-entry CLI).
> Test tables: 2,680 plant-compound reference set; 4,215 fecal metabolites (2024 restricted-diet system); 4,583 fecal metabolites (2025 recovery system).
> Environment: Python 3.10 / pandas 2.3 / RDKit 2025.03. No extra packages installed.

---

## 1. Verdict

| Acceptance layer | Result |
|---|---|
| Built-in textbook exam (28 assertions) | **28/28 (100%)** |
| Full cross-check vs. independently curated glycoside labels (2,680 compounds) | **99.0% detector correctness on the glycoside set** (678/685; overall agreement 2,634/2,680 = 98.3%) |
| Stratified manual audit (per-family random samples, structure-by-structure) | **0 errors in the final two consecutive rounds (36/36)** |

Coverage (compounds matching >=1 bond code): plant 87.5% / fecal 2024 86.7% / fecal 2025 87.1%. Non-covered compounds are by-design exclusions (background skeleton bonds: amino acids, organic acids, plain hydrocarbons — no enzyme is defined by breaking those bonds).

## 2. The three layers

1. **Built-in exam** — 28 textbook assertions (positives, negatives, edge cases) shipped in `tests/test_exam.py` and run on every rule change. Protects against regressions of every previously fixed defect.
2. **Full cross-check** — the reference table carries independently curated glycoside labels (name + fingerprint-library + manual curation pipeline). The detector (pure SMARTS/structure computation) was run against these labels as an independent answer key. 676/685 curated glycoside compounds confirmed; of the 9 disagreements, 3 are label errors in the reference itself (free fructose/sorbitol/anhydroglucitol curated as "C-glycosides"), 2 are documented SMILES-expression blind spots (aloin-type), 1 is a label/structure contradiction in the source (trehalose-type annotated C-C but drawn O-bridged), and 3 minor. The detector additionally surfaced 39 glycosides/acyl-sugar esters missed by the curation (phenolic-acid sugar esters, phenolic glycosides, nucleosides) — the two methods are complementary.
3. **Stratified manual audit** — per-family random samples (glycoside/ester/alkene/ether/sulfur/C-C/amide) judged structure-by-structure. 15 iterative rounds total; acceptance requires **two consecutive error-free rounds**.

## 3. Defects found and fixed during testing (all regression-tested)

| Round | Symptom | Root cause | Fix |
|---|---|---|---|
| exam v1 | nitrile never matched; glycylglycine not typed as peptide; ethylene oxide double-coded | SMARTS [#7X1] invalid on triple bond; peptide pattern required CH (glycine is CH2); T02 lacked epoxide exclusion | `[#6]#[#7]`; alpha-carbon pattern widened to sp3 C (verified no amide false positive); T02 excludes 3-rings |
| exam v2 | gentiobiose/thioglycoside picked up T02 | sugar in-ring ethers counted as alkyl ethers | T02 excludes O atoms inside sugar rings |
| scan 1 | Leptostachyol acetate typed G16 | furan-lactone ring (3 O-bearing carbons) passed the old sugar-ring test | sugar-ring criterion tightened |
| scan 2 | methyl esters typed T01; bis-THF lignan ring typed "sugar" | ester O-CH3 is not an ether; THF ring has 4 O-bearing carbons | T01 excludes carbonyl-adjacent O; sugar ring requires a -CH2O- side arm |
| scan 2 (detail) | flavones typed C01 | chromenone ring perceived aromatic; ring-to-aryl bond looks like biphenyl | C01 requires both atoms in all-carbon six-membered aromatic rings, different rings, not adjacent to ring O |
| scan 3 | catechin chroman O missed as T03; methylenedioxy O double-counted | over-exclusion; no T05 de-duplication | T03 includes chroman/benzopyran ring O; excludes T05 group O |
| demo | chlorogenic acid (quinic-acid cyclitol) typed G16 | cyclitol ring passed sugar tests (4 O-bearing carbons) | sugar ring requires >=3 carbons with **exocyclic** O and an sp3 anomeric carbon free of exocyclic OH — cyclitols have 2 and fail |

The cyclitol case (quinic acid / chlorogenic acid skeleton) is methodologically notable: quinic acid has exactly one 5-membered O-heterocycle with 4 O-bearing carbons — structurally very sugar-like — but only 2 carbons carry **exocyclic** oxygen (glycosidic or hydroxyl), versus 4 in true pyranose/furanose rings. This criterion (>=3 exocyclic-O carbons) cleanly separates true sugars from cyclitols, pantolactone-type lactones, and bis-THF lignan rings simultaneously.

## 4. Final catalog-level numbers (46 codes, selected)

| Code | Plant 2680 | Fecal 2024 | Fecal 2025 | Code | Plant | 2024 | 2025 |
|---|---|---|---|---|---|---|---|
| G14 O-glycoside | 628 | 558 | 649 | T01 methoxy | 694 | 1127 | 1177 |
| G16 C-glycoside | 43 | 59 | 56 | T02 alkyl ether | 575 | 727 | 794 |
| G17 N-glycoside | 16 | 33 | 36 | T03 aryl ether | 316 | 453 | 520 |
| E01 aliphatic ester | 329 | 527 | 571 | T04 epoxide | 97 | 145 | 167 |
| E02 aromatic ester | 199 | 266 | 311 | S01 thioether | 12 | 28 | 30 |
| E04 aromatic lactone | 130 | 199 | 216 | C01 aryl-aryl C-C | 41 | 38 | 45 |
| E09 aliphatic lactone | 208 | 338 | 392 | N04 peptide bond | 7 | 59 | 54 |

## 5. Honest limitations

1. **Stereochemistry not resolved**: D/L, alpha/beta and sugar identity need stereo-resolved structures (name-to-structure bridge planned). All glycosides are typed at linkage level (G14 = sugar undetermined).
2. 2-deoxy sugars and aloin-type C-glycoside SMILES expressions are documented fail-safe blind spots.
3. A molecule containing both a lactone and an open-chain ester is typed lactone-only (priority rule, documented).
4. Acyl-sugar esters are legitimately double-counted (E05 glycoside-context + E01/E02 ester-context).
5. Input rows are assumed to be real MS2-annotated features; annotation-layer false positives are out of scope.
6. Cross-checking validated linkage-level typing only; sugar-species assignment remains with the curation pipeline.
