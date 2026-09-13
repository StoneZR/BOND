# BOND: Bond-centric Omics Nexus Decipherer

**Bond-type annotation engine that turns metabolite SMILES into a 46-code chemical-bond profile — the bridge between metabolomics and metagenomics.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![RDKit](https://img.shields.io/badge/RDKit-required-green)

---

## Why BOND

Metabolomics and metagenomics speak different languages: compounds on one side, enzymes on the other. Existing integration tools (MIMOSA, MelonnPan, MMvec) correlate abundances statistically but cannot explain **which chemical bond of which compound is acted on by which enzyme**. BOND makes the **chemical bond** the shared unit of both omics layers:

```
metabolite structure ──► bond-type profile (46 codes)     [this repo]
metagenome (KO/CAZy) ──► enzyme→bond action spectrum       [next layer]
                         ▼
        compound × bond-type × enzyme association
```

No single existing tool provides bond-level, traceable, enzyme-anchored annotation (see `docs/evidence_report.md` for a four-database analysis: IUBMB/KEGG/BRENDA/CAZy).

## What it does

For every input SMILES, BOND reports:

1. **Glycosidic linkage typing** — O- / C- / N- / S-glycoside, acyl-sugar ester, with blind-spot guards (see Accuracy);
2. **25 non-glycosidic bond codes** — carboxylic esters (aliphatic/aromatic/galloyl), lactones (aromatic/aliphatic), phosphate/sulfate/thio esters, amides/lactams/peptide bonds, nitriles, activated/cinnamoyl/conjugated/isolated alkenes, alkynes, methoxy/alkyl/aryl ethers, epoxides, methylenedioxy bridges, thioethers, disulfides, aryl–aryl C–C bonds;
3. **Full traceability** — every rule is a named SMARTS with a documented chemical basis and EC/CAZy anchor in the catalog.

## Quick start

```bash
pip install -r requirements.txt   # pandas + rdkit

# CLI: annotate a table (any CSV/TSV with a SMILES column)
python -m bond.judge_core examples/example_input.csv output.tsv --smiles smiles --name name

# or as a library
```

```python
from bond import BondJudge

judge = BondJudge()
res = judge.judge("OC[C@H]1O[C@@H](Oc2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O")
print(sorted(res.bond_codes))        # ['G14']  -> O-glycoside (sugar species undetermined)
print(res.glycosidic_linkage)        # 'O'
```

Run the built-in acceptance exam (28 textbook assertions):

```bash
python tests/test_exam.py            # built-in exam: 28/28 passed
```

## The bond catalog

`data/bond_catalog_v1.4.tsv` — 46 bond codes in 7 families, each with:

- structural definition + SMARTS,
- **chemical anchor** (IUPAC functional-group terminology, carbon oxidation-state ladder),
- **enzymatic anchor** (IUBMB EC sub-subclass + CAZy family + KEGG KO),
- reaction direction (degradation/synthesis/oxidation/reduction),
- evidence grade A/B/C (A = dedicated EC entries exist for this bond context).

Bond families are the intersection of two official systems acting on **the bond being broken**: IUPAC functional groups (chemistry side) × EC "acting on ... bonds" subclasses (enzymology side). New family admission requires both anchors. Details: `docs/classification_report.md`.

## Accuracy (three-layer acceptance)

| Layer | Method | Result |
|---|---|---|
| Built-in exam | 28 textbook assertions (positive/negative/edge), run on every change | **28/28** |
| Full cross-check | 2,680-compound reference table with independently curated glycoside labels | **99.0%** detector correctness on the glycoside set (678/685); 3/6 remaining misses are label errors in the reference, each documented |
| Stratified manual audit | Random samples per bond family, structure-by-structure verification over 15 iterative rounds | **0 errors in the final two consecutive rounds (36/36)** |

Iterative testing surfaced and fixed 10 systematic rule defects (aromatic-perceived lactones, flavone pyranone rings mimicking biphenyls, cyclitol rings mimicking sugars, methyl esters masquerading as methoxy ethers, and more) — every fix is regression-tested. Full log: `docs/accuracy_report.md`.

## Scope & honest limitations

- **Stereochemistry is not resolved**: D/L, α/β and sugar identity require stereo-resolved structures (name→structure bridging via OPSIN is planned). All glycosides are currently typed at linkage level (`G14` = O-glycoside, sugar undefined).
- 2-deoxy sugars and aloin-type C-glycoside SMILES expressions are documented blind spots (fail-safe: omitted rather than mis-typed).
- Input rows are assumed to be real MS2-annotated features; annotation-layer false positives are out of scope.

## Repo layout

```
bond/           core package (judge_core.py — bond judge + table runner)
data/           bond_catalog_v1.4.tsv (46 codes, dual anchors, evidence grades)
tests/          test_exam.py (28 built-in assertions)
scripts/        batch verification helpers
docs/           classification & accuracy & enzyme-evidence reports
examples/       demo input/output
```

## Citation

If you use BOND, please cite the bond catalog (dual-anchor classification of 46 bond codes) and the accuracy protocol (three-layer acceptance). A methods paper is in preparation.

## License

MIT — see [LICENSE](LICENSE).
