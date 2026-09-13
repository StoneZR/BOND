# BOND Bond-Type Classification: Closure Check Report (v1.4)

> Date: 2026-09-13 (English release version of the working report v1.4)
> Environment: Python 3.10 / pandas 2.3 / RDKit 2025.03 (no extra packages)
> Companion files: `data/bond_catalog_v1.4.tsv` (46-code catalog), `docs/accuracy_report.md`, `docs/evidence_report.md`

---

## 1. Summary

1. **7 bond families × 46 bond codes** (two-layer classification): Glycosidic 18, Ester 9, Amide-Nitrile 4, Alkene-Alkyne 5, Ether-Epoxide 5, Sulfur 2, C-C skeleton 3.
2. **Compound-side coverage (measured)**: 2,680 plant compounds — 2,374 (88.6%) matched at least one bond code; 4,215 fecal metabolites — 3,678 (87.3%); 4,583 recovery-phase metabolites — 87.1%.
3. **Enzyme-side closure**: every official IUBMB "acting on ... bonds" subclass (3.1–3.8, plus EC 1.3/1.8/2.5.1/4.2.1 for non-hydrolytic bonds) maps to a catalog row or is documented as out-of-domain.
4. **Positive-control regression**: 22/22 textbook enzyme-substrate pairs map correctly.
5. **Evidence grading**: A (dual official anchors) 28 codes, B (single anchor + literature) 15, C (marker-only) 3.

## 2. Why these families (the dual-anchor rule)

Each family must hold **two official credentials**:

- **Chemistry anchor**: an IUPAC functional-group term (Gold Book / Blue Book: ether, ester, acetal, amide, nitrile, ...) placed on the carbon oxidation-state ladder;
- **Enzymology anchor**: an IUBMB EC subclass officially defined by the bond it acts upon (EC 3.1 esters, 3.2 glycosylases, 3.3 ethers, 3.4 peptide bonds, 3.5 C-N bonds incl. 3.5.5 nitriles; non-hydrolytic bonds covered by EC 1.3, 1.8, 2.5.1, 4.2.1).

Why both: chemistry alone cannot justify "why these bonds matter" (no enzyme landing point for omics association); enzymology alone cannot ground the annotation in structure (no SMARTS realizability). The intersection of the two systems equals the space of bonds that are **structurally definable and enzymologically meaningful**.

| Family | Chemistry anchor (IUPAC) | Enzymology anchor (EC; entries measured in local IUBMB 2026-06 release) |
|---|---|---|
| Glycosidic | acetal / glycoside | EC 3.2 (3.2.1: 232; 3.2.2: 32; 3.2.3: 1) |
| Ester | ester (incl. phosphate/sulfate/thioester) | EC 3.1 (3.1.1: 126; 3.1.2: 35; 3.1.3/3.1.4: 181; 3.1.6: 22) |
| Amide-Nitrile | amide; nitrile | EC 3.4 (557) + 3.5 (3.5.5 nitrile: 8) |
| Alkene-Alkyne | alkene (+ conjugation criteria) | EC 1.3 (131) + 4.2.1 hydratases + 1.13.11 dioxygenases |
| Ether-Epoxide | ether | EC 3.3 (18; main: 3.3.2 epoxide hydrolase) + 1.14 CYP oxidative cleavage |
| Sulfur | thioether; disulfide | EC 1.8 (55) + 3.3.1 (rare) |
| C-C skeleton | biphenyl-type | EC 1.10.3 laccase / 1.11.1 peroxidase (radical); synthesis anchor EC 2.5.1 |

**Why exactly 7**: fewer families would force alkene/sulfur/C-C into other families (violating the dual-anchor rule); more families (e.g. promoting epoxide or nitrile) is impossible — they are sub-subclasses, not top EC positions. 7 is derived, not chosen.

**Chemical corroboration (oxidation-state ladder)**: ether, ester and glycosidic bonds are all C-O bonds; the split follows the oxidation state of the O-bonded carbon and the leaving-group ability (ether: pKa-16 alkoxide leaving group — almost no hydrolases exist, CYP oxidative cleavage instead; ester: pKa-5 carboxylate — EC 3.1 hydrolases are ubiquitous; glycoside: acid-labile acetal — EC 3.2 glycosidases use general-acid catalysis). Chemistry and enzymology point to the same boundary.

**Admission rule for a new family**: must present both ① an IUPAC official term and ② at least one EC bond-acting subclass. Missing either → the item stays a code within a family or a marker-only type.

## 3. Why the codes are unevenly distributed across families

Measured family resolution of the IUBMB database:

| Family (EC anchor) | 4th-level entries | Entries naming stereo/config | Named-context examples |
|---|---|---|---|
| O-glycosidases (3.2.1) | 232 | **166 (72%)** | 118 entries named by sugar species (glucosidase, rhamnosidase, xylosidase, ...) |
| N-glycosidases (3.2.2) | 32 | 8 (25%) | all nucleosidase variants — one bond context |
| S-glycosidases (3.2.3) | **1** | 0 | myrosinase only |
| Carboxylesterases (3.1.1) | 126 | 38 (30%) | galloyl/feruloyl/coumaroyl/tannin/lactonase contexts (64) |
| Amides (3.5.1/2) | 162+20 | 49 (30%) | named by substrate, one bond context |
| Peptidases (3.4) | 557 | 144 (26%) | named by protein substrate |
| Ene reductases (1.3.1) | 131 | 49 (37%) | named by hydride donor, not alkene context |

The imbalance is the real structure of enzymology databases: plants decorate metabolites with "sugar × configuration × position" diversity and microbes co-evolved matching specialist glycosidase arsenals, whereas ethers/thioethers are chemically inert and never spawned enzyme-family diversity. **Stereochemistry in non-glycosidic families exists chemically (alkene E/Z, epoxide R/S) but fails three gates**: no EC entry is defined by that configuration (E/Z isomerases belong to EC 5.2, outside "bond-acting"; epoxide enantioselectivity is a within-entry property), and the data lack stereo-resolved structures. These are documented as attribute-layer upgrades pending a name-to-structure bridge.

## 4. Natural-product coverage

Measured on the 2,680-compound reference table (superclass-level):

| Superclass | Total | Hit | Coverage |
|---|---|---|---|
| Phenylpropanoids (flavonoids/lignans/coumarins/tannins) | 945 | 898 | 95% |
| Terpenoids | 654 | 616 | 94% |
| Lipids | 315 | 252 | 80% |
| Alkaloids | 278 | 229 | 82% |
| Polyketides | 178 | 151 | 85% |
| Steroids | 71 | 60 | 85% |
| Carbohydrates | 61 | 50 | 82% |
| Organic acids / Amino acids | 99 | 51 | 52% (background bonds, by design) |
| **Total** | **2680** | **2374** | **88.6%** |

All 15 major natural-product classes (flavonoids, lignans, coumarins, tannins, terpenes, biflavonoids, alkaloids, saponins, cyanogenic glycosides, glucosinolates, polyacetylenes, polyketide macrolides, phosphorylated metabolites, organosulfur compounds, prenylated derivatives) have their hallmark bonds covered by at least one catalog code or explicitly excluded with documented reasons (background C-C/C-H/C-O skeleton bonds: no enzyme is defined by breaking them; their enzymatic chemistry proceeds by site-level oxidation/hydroxylation, EC 1.14 — 845 entries, all named by substrate/position, none by bond).

## 5. Closure checks (four gates, all passed)

1. **Compound-side closure**: every NP-class hallmark bond is either a catalog code or has an exclusion record.
2. **Enzyme-side closure**: every EC 3.1–3.8 bond-acting subclass maps to ≥1 code or carries a documented domain exemption (3.1.5/3.1.7/3.1.8, 3.1.11–3.1.31, 3.6, 3.7.1, 3.8 halides — halogenated NPs are a marine/fungal domain, near-zero frequency in terrestrial plant studies).
3. **Positive-control regression**: 22/22 textbook pairs (beta-glucosidase→G01, tannase→E03, epoxide hydrolase→T04, peptidase→N04, gluconolactonase→E09, C-glycoside oxidase 1.1.3.50→G16, laccase→C01, ...).
4. **No dangling codes**: 43 official codes all carry both a chemical definition and an enzymatic anchor; 3 marker-only codes (isolated alkene, alkyne, prenyl C-C) are explicitly excluded from enzyme association.

## 6. Feasibility & credibility

- **Feasibility measured**: full scan of 2,680+4,215+4,583 compounds runs in ~2 minutes in RDKit; no extra packages.
- **Credibility grades**: A 28 codes (dual official anchors), B 15 (single anchor + literature), C 3 (marker-only). Grade A requires a dedicated EC entry defined by the bond context (beta-glucosidase EC 3.2.1.21 is "born for" code G01, whereas CYP O-demethylation for T01 is a broad-specificity passenger, hence B).
- **Known limitations**: sugar species and alpha/beta D/L configuration are not resolved (name-to-structure bridging planned); aromatic-fused lactone perception required element-wildcard SMARTS; the T02 count includes sugar in-ring ethers (documented); annotation-layer false positives are out of scope.

## 7. Versioning

v1.0 (44 codes) → v1.1 (peptide bond split; lactone aromatic/aliphatic split; aromatic-perception fix) → v1.2 (essence/context columns; EC 3.8 exemption) → v1.3 (C-glycoside anchors corrected to EC 1.1.3.50/4.1.99.28) → v1.4 (linkage-atom typing; five-rule sugar-ring sharpening: all-aliphatic ring carbons, >=4 carbons with O, sp3 anomeric carbon, single-bond linkage, carbonyl-checked ester routing). Each change is logged with its chemical basis and enzymatic split point in the versioned curation records (44 entries).
