# -*- coding: utf-8 -*-
"""Built-in exam: 28 textbook assertions (positive / negative / edge cases).

Run standalone:  python tests/test_exam.py
Run with pytest: pytest tests/test_exam.py -v

These are the acceptance assertions from iterative accuracy testing
(see docs/accuracy_report.md). Each case carries the lesson tag of the
rule it protects.
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from bond.judge_core import _judge  # noqa: E402

EXAM = [
    # (SMILES, expected bond codes, name, lesson tag)
    ('OC[C@H]1O[C@@H](Oc2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O', {'G14'},
     'phenyl beta-D-glucoside', '-'),
    ('CC1(COC(=O)C1O)C', {'E09'}, 'pantolactone (aliphatic lactone, not glycoside)', '#5'),
    ('O=C1OC2=CC=CC=C2C=C1', {'E04'}, 'coumarin (aromatic lactone)', '#1'),
    ('C1[C@H]([C@@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O', {'T03'},
     'catechin (chroman ring O = aryl-alkyl ether; no sugar/ester)', '#4'),
    ('O=C(O)C1=CC=C(O)C=C1', set(), 'p-hydroxybenzoic acid (carboxyl: not a cleavable bond)', 'excl'),
    ('COC1=CC(=CC=C1OC)C(=O)OC2C(C(C(C(O2)CO)O)O)O', {'E02', 'T01', 'E05'},
     'benzoic acid glucose ester (aromatic ester + methoxy + acyl-sugar)', '#6'),
    ('C/C=C/C', {'A04'}, '2-butene (isolated alkene)', '-'),
    ('O=C(/C=C/c1ccccc1)O', {'A01', 'A02'}, 'cinnamic acid (cinnamoyl alkene = activated alkene)', '-'),
    ('COc1ccccc1OC', {'T01'}, 'veratrole (methoxy ethers)', '-'),
    ('C1=CC=C2C(=C1)C=CO2', set(), 'benzofuran core (in-ring aromatic ether: not T03)', '-'),
    ('SC1C(C(C(C(O1)CO)O)O)O', {'G18'}, 'S-glucoside (thioglycoside, blind spot 17a)', '#17a'),
    ('CCCCCCCCCCCCCCC(=O)OCC1OC(=O)C(O)C(O)C1O', {'E01'},
     'fatty acid glucose diester (aliphatic ester)', '-'),
    ('NCC(=O)NCC(=O)O', {'N04'}, 'glycylglycine (peptide bond)', '-'),
    ('CC(C)CC1=CC=C(C=C1)C(C)C(=O)N', {'N01'}, 'ibuprofen amide (open-chain, non-peptide)', '-'),
    ('C1CC(=O)NC1=O', {'N02'}, '2-pyrrolidone (lactam)', '-'),
    ('CC(C)(C)c1ccc(cc1)C(C)C(=O)O', set(), 'ibuprofen (carboxylic acid)', 'excl'),
    ('C1=COC=C1', set(), 'furan (aromatic in-ring ether: not T02/T03)', '-'),
    ('C1OC1', {'T04'}, 'ethylene oxide (epoxide)', '-'),
    ('COc1ccc2c(c1)OCO2', {'T05', 'T01'}, 'methylenedioxybenzene + methoxy', '-'),
    ('c1ccc(-c2ccccc2)cc1', {'C01'}, 'biphenyl (aryl-aryl C-C)', '-'),
    ('OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@H]1O', {'G14'},
     'gentiobiose (sugar-sugar O-glycoside; in-ring ethers not T02)', '-'),
    ('C1CCC2C1CCC3C2CCC4C3CCC(C4)(C)C', set(), 'sterane skeleton (background bonds only)', 'excl'),
    ('OC(=O)C(=O)OCC1OC(=O)C(O)C(O)C1O', {'E01'}, 'malonyl glucose (aliphatic ester on sugar)', '-'),
    ('C1=CC(=CC=C1C=CC(=O)OCC2C(C(C(C(O2)CO)O)O)O)O', {'A01', 'A02', 'E02'},
     'caffeoyl glucose C6 ester (cinnamoyl aromatic ester; no anomeric glycoside)', '-'),
    ('C1CCCS1', {'S01'}, 'thiolane (thioether)', '-'),
    ('CC#N', {'N03'}, 'acetonitrile (nitrile)', '-'),
    ('C#CC', {'A05'}, 'propyne (alkyne)', '-'),
    ('C1CCC(C)CC1', set(), 'methylcyclohexane (background bonds)', 'excl'),
]


def run_exam():
    fails = []
    for smi, exp, name, tag in EXAM:
        hits, _, _ = _judge(smi)
        if hits != exp:
            fails.append((name, sorted(exp), sorted(hits), smi[:50]))
    return fails


if __name__ == '__main__':
    fails = run_exam()
    for name, exp, got, smi in fails:
        print(f'[FAIL] {name}: expected {exp} got {got} | {smi}')
    print(f'built-in exam: {len(EXAM) - len(fails)}/{len(EXAM)} passed')
    sys.exit(1 if fails else 0)


def test_exam():
    fails = run_exam()
    assert not fails, f'{len(fails)} exam failures: {fails}'
