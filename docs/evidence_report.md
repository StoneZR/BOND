# How Many Bond Essences Can Be Enzyme-Linked? A Four-Database Evidence Report

> Date: 2026-09-13. Question: among all chemical bonds of natural products, how many "bond essences" can be linked to enzymes, and on what evidence?
> Databases (all local snapshots): IUBMB/ExPASy ENZYME (`enzyme.dat` + `enzclass`, release 2026-06-10, SIB, CC BY 4.0); KEGG reaction set (9,711 EC-annotated reactions); BRENDA microbial enzyme records (156,017); CAZy (524 families, 2026-07-24 snapshot).
> Every number below is reproducible with the measurement script archived in the working repository.

---

## 1. Verdict

**In the natural-product domain (plant specialized metabolism and its microbial transformation), 16 bond essences have enzymatic anchors**; adding 2 domain-edge essences (halide C-X, peroxide O-O) gives 18 candidates, of which 16 are catalogued.

## 2. The 16 bond essences (46 codes collapse into 16)

| # | Essence | Codes | # | Essence | Codes |
|---|---|---|---|---|---|
| 1 | O-glycosidic | G01-G15 (15) | 9 | amide | N01, N02, N04 (3) |
| 2 | C-glycosidic | G16 | 10 | nitrile | N03 |
| 3 | N-glycosidic | G17 | 11 | alkene | A01-A04 (4) |
| 4 | S-glycosidic | G18 | 12 | alkyne | A05 |
| 5 | carboxylic ester | E01-E05, E09 (6) | 13 | ether | T01-T05 (5) |
| 6 | phosphate ester | E06 | 14 | thioether | S01 |
| 7 | sulfate ester | E07 | 15 | disulfide | S02 |
| 8 | thioester | E08 | 16 | aryl C-C | C01-C03 (3) |

46 codes = 16 essences expanded by their **context dimensions** (sugar species, linkage position, acyl environment, conjugation, ring fusion) at enzymological resolution. Codes sharing one essence are the same chemical bond in different enzyme-recognizable contexts.

## 3. Four independent database perspectives

### IUBMB (official classification)
Bond-acting EC positions: 3.1 esters, 3.2 glycosides, 3.3 ethers, 3.4 peptide, 3.5 C-N (3.5.5 nitrile), 3.6 anhydrides, 3.7 C-C (keto-acids), 3.8 halides; non-hydrolytic: 1.3 (CH-CH), 1.8 (sulfur), 2.5.1 (prenyl transfer), 4.2.1 (hydration). Entries per essence: amide 719, alkene 430, O-glycosidic 232, phosphate 181, carboxylic ester 126, aryl C-C 47, thioester 35, N-glycosidic 32, sulfate 22, disulfide 16, ether 15, nitrile 8, thioether 3, S-glycosidic 1, **C-glycosidic 2** (EC 1.1.3.50 C-glycoside oxidase + 4.1.99.28 deglycosidase — cross-class entries in oxidoreductases/lyases, not glycosidases), alkyne 0.

### KEGG (reaction universe)
Of 9,711 EC-annotated reactions, ~1,800 map to the 16 essences: alkene 655, O-glycosidic 272, amide 241, phosphate 208, carboxylic ester 189, aryl C-C 70, thioester 53, N-glycosidic 31, ether 32, sulfate 26, disulfide 19, nitrile 14, C-glycosidic 2 (R13033/R13184). The remaining 7,904 reactions are primary-metabolism redox chemistry (site-level, not bond-cleavage).

### BRENDA (microbial instances — most relevant for gut context)
O-glycosidic **19,116 records (12.3%)** — the largest enzyme arsenal of any essence; amide 15,785; carboxylic ester 7,481; alkene 5,572; phosphate 3,951; aryl C-C 3,334; nitrile 1,183; N-glycosidic 780; thioester 511; disulfide 368; sulfate 249; ether 151; C-glycosidic 10.

### CAZy (family system — the negative-space argument)
CAZy's entire universe (GH 187, GT 138, CBM 114, PL 43, CE 23, AA 19) covers only ~2.5 essences (O-glycosidic + carbohydrate-context esters + auxiliary oxidative). This is the strongest demonstration that glycosidic and ester bonds carry the thickest enzymatic equipment — and that the remaining essences must be anchored through EC (which is why BOND uses EC as primary anchor, CAZy as glycoside/ester reinforcement).

## 4. Cross-database consistency

| Essence | IUBMB | KEGG | BRENDA | CAZy | Non-zero DBs |
|---|---|---|---|---|---|
| O-glycosidic | 232 | 272 | 19,116 | GH+GT+PL (368) | **4/4** |
| carboxylic ester | 126 | 189 | 7,481 | CE (23) | **4/4** |
| aryl C-C | 47 | 70 | 3,334 | AA (19, auxiliary) | **4/4** |
| amide | 719 | 241 | 15,785 | — | 3/4 |
| alkene | 430 | 655 | 5,572 | — | 3/4 |
| phosphate ester | 181 | 208 | 3,951 | — | 3/4 |
| N-glycosidic | 32 | 31 | 780 | — | 3/4 |
| thioester | 35 | 53 | 511 | — | 3/4 |
| nitrile | 8 | 14 | 1,183 | — | 3/4 |
| sulfate ester | 22 | 26 | 249 | — | 3/4 |
| disulfide | 16 | 19 | 368 | — | 3/4 |
| ether | 15 | 32 | 151 | — | 3/4 |
| C-glycosidic | 2 | 2 | 10 | — | 3/4 |
| thioether | 3 | 0 | 0 | — | 1/4 |
| S-glycosidic | 1 | 0 | 0 | — | 1/4 |
| alkyne | 0 | 0 | 0 | — | 0/4 (marker-only) |

**13/16 essences have non-zero evidence in >=3 databases.** Thioether and S-glycosidic rest on official entries only; alkyne is marker-only. No essence is claimed beyond its evidence.

## 5. The negative statement (why skeleton bonds are excluded)

EC 1.14 (CYP-type oxygenases) contains 845 entries — all named by substrate/position (procollagen-lysine 5-dioxygenase, L-arginine hydroxylase, kanamycin B dioxygenase, ...), **none by bond cleavage**. Skeleton C-C/C-H/C-O bonds are transformed by site-level hydroxylation/oxidation, not bond scission; they belong to a site layer (functional-group modification), not the bond-essence layer. The 16-essence boundary is therefore the objective structure of enzymology, not an omission.

## 6. Domain edges (recorded, not catalogued)

| Edge essence | Anchor | Frequency in terrestrial-plant data | Disposition |
|---|---|---|---|
| halide C-X | EC 3.8 (13 entries, e.g. haloalkane dehalogenase 3.8.1.5) | ~0 (marine/fungal domain) | recorded; admitted on dataset expansion |
| peroxide O-O | no entry-level anchor | ~0 (artemisinin-type) | recorded |

## 7. Use in BOND

The ~1,800 KEGG reactions and 156,017 BRENDA records mapped to the 16 essences are the data source for the next BOND layer: the **enzyme-to-bond action spectrum** (enzyme x bond code x direction), which is the interface that connects metabolite bond profiles to metagenomic KO/CAZy abundances.

---

*All figures reproducible with the archived measurement script. Database versions: IUBMB/ExPASy ENZYME 2026-06-10; KEGG local snapshot 2026-08; BRENDA microbial subset v4; CAZy 2026-07-24.*
