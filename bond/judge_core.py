# -*- coding: utf-8 -*-
"""BOND bond-type annotation core.

Two-layer bond classification (7 families / 46 bond codes):
glycosidic linkage typing (O/C/N/S/ester, blind spots patched) and
26 non-glycosidic SMARTS rules (esters, lactones, alkenes, ethers,
epoxides, amides, nitriles, sulfur, aryl-aryl C-C).

All rules encode lessons learned from iterative accuracy testing
(see docs/accuracy_report.md). Each rule carries its lesson tag.
"""
import os, sys, re, argparse, dataclasses
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ==================== 1. Glycosidic linkage typing (lessons #3,4,5,6,7,17a,17b,17d) ====================
def detect_glycoside(mol):
    """Type the glycosidic linkage of a molecule.

    Returns dict(link=O/C/N/S/ester/free/None, n_rings=int, quinone_c=bool).
    Rules (each tagged with its accuracy-testing lesson):
      R1 sugar ring = 5/6-membered, exactly one O, rest C          (ring definition)
      R2 ring carbons must be all-aliphatic (lesson #4: flavone benzopyran ring mimics)
      R3 >=3 ring carbons carry an exocyclic O (true sugar hallmark; 2-deoxy sugars
         are a documented blind spot 17c, fail-safe excluded)
      R4 anomeric carbon must be sp3 (lesson #5: lactone carbonyl C mimics)
      R5 linkage must be a single bond (lesson #5: double-bond O mimics)
      R6 exocyclic O without H: carbonyl end -> acyl-sugar ester; else O-glycoside (lesson #6)
      R7 exocyclic O with H -> free reducing sugar, not counted as glycoside (lesson #7)
      R8 C linkage: aromatic C -> C-glycoside; quinone-ring sp2 C -> aloin-type (17b)
      R9 N/S linkage -> N- / S-glycoside (blind spot 17a)
    """
    ri = mol.GetRingInfo()
    n_sugar = 0
    link_found = None
    quinone_c = False
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6): continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        os_ = [a for a in atoms if a.GetSymbol() == 'O']
        cs_ = [a for a in atoms if a.GetSymbol() == 'C']
        if len(os_) != 1 or len(cs_) != len(ring) - 1: continue          # R1
        if any(a.GetIsAromatic() for a in cs_): continue                 # R2
        # R3 (v2.5, synced with is_sugar_ring): >=3 ring carbons carry exocyclic O
        # (anomeric glycosidic O or hydroxyl); excludes pantolactone rings, quinic-acid
        # cyclitols (2), and 2-deoxy sugars (blind spot 17c)
        n_exo_o = sum(1 for c in cs_ if any(nb.GetSymbol() == 'O' and nb.GetIdx() != os_[0].GetIdx()
                                            for nb in c.GetNeighbors()))
        if n_exo_o < 3: continue
        n_sugar += 1
        o_idx = os_[0].GetIdx()
        for c in cs_:
            if c.GetHybridization() != Chem.HybridizationType.SP3: continue   # R4
            if o_idx not in [nb.GetIdx() for nb in c.GetNeighbors()]: continue
            # v2.1 (lesson #17g, quinic-acid case): a true anomeric carbon carries no
            # exocyclic hydroxyl O; cyclitol alcohol carbons do, so their aromatic
            # neighbors do not imply a C-glycoside
            exo_oh = any(x.GetSymbol() == 'O' and x.GetIdx() != o_idx and x.GetTotalNumHs() >= 1
                         for x in c.GetNeighbors())
            if exo_oh:
                continue
            for nb in c.GetNeighbors():
                if nb.GetIdx() == o_idx: continue
                b = mol.GetBondBetweenAtoms(c.GetIdx(), nb.GetIdx())
                if b.GetBondType() != Chem.BondType.SINGLE: continue          # R5
                sym = nb.GetSymbol()
                if sym == 'O':
                    if nb.GetTotalNumHs() >= 1:
                        link_found = link_found or 'free'                     # R7
                    else:
                        for ob in nb.GetNeighbors():
                            if ob.GetIdx() == c.GetIdx(): continue
                            if ob.GetSymbol() == 'C':
                                for bb in ob.GetNeighbors():
                                    if bb.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(
                                            ob.GetIdx(), bb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                        link_found = link_found or 'ester'; break   # R6
                            elif ob.GetSymbol() == 'P' or ob.GetSymbol() == 'S':
                                link_found = link_found or 'ester'                # phospho/sulfate sugar ester
                        link_found = link_found or 'O'
                elif sym == 'C':
                    if nb.GetIsAromatic():
                        link_found = link_found or 'C'                        # R8 aryl C-glycoside
                    elif nb.GetHybridization() != Chem.HybridizationType.SP3:
                        # blind spot 17b: aloin-type C-glycoside (anomeric C bonded to non-aromatic sp2 carbon)
                        quinone_c = True
                        link_found = link_found or 'C'
                    else:
                        # blind spot 17f: trehalose-type sugar-sugar C-C bridge (outer C in a different sugar ring)
                        nb_rings = [r for r in mol.GetRingInfo().AtomRings() if nb.GetIdx() in r]
                        if any(is_sugar_ring(mol, r) for r in nb_rings) and \
                           not any(set(r) == set(ring) for r in nb_rings):
                            link_found = link_found or 'C'
                elif sym == 'N':
                    link_found = link_found or 'N'                            # R9
                elif sym == 'S':
                    link_found = link_found or 'S'                            # R9
    return {'link': link_found, 'n_rings': n_sugar, 'quinone_c': quinone_c}

# ==================== 2. Non-glycosidic SMARTS (26 rules; lessons #1 wildcard form, #2 tuple index) ====================
SP = {
 'E01_aliphatic_ester': r'[CX3](=[OX1])[OX2H0][CX4]', 'E02_aromatic_ester': r'c[CX3](=[OX1])[OX2H0]',
 'E03_galloyl': r'c1c(O)c(O)c(C(=[OX1]))c(O)c1',
 'E04_aromatic_lactone': r'[#6X3](=[#8X1])[#8X2;R]~[c]', 'E09_aliphatic_lactone': r'[#6X3](=[#8X1])[#8X2;R]',
 'E06_phosphate': r'[PX4](=[OX1])([OX2])[OX2]', 'E07_sulfate': r'[SX4](=[OX1])(=[OX1])([OX2])[OX2]',
 'E08_thioester': r'[#6X3](=[#8X1])[#16X2][#6]',
 'N01_openchain_amide': r'[#6X3](=[#8X1])[#7X3;!R]', 'N02_lactam': r'[#6X3](=[#8X1])[#7X3;R]',
 'N04_peptide': r'[#6X4][#6X3](=[#8X1])[#7X3][#6X4]',
 'N03_nitrile': r'[#6]#[#7]',
 'A01_activated_alkene': r'[CX3]=[CX3][CX3]=[OX1]', 'A02_cinnamoyl_alkene': r'[c][CX3]=[CX3][CX3]=[OX1]',
 'A03_conjugated_diene': r'[CX3]=[CX3][CX3]=[CX3]', 'A05_alkyne': r'[#6X2]#[#6X2]',
 'T01_methoxy': r'[OX2][CH3]', 'T04_epoxide': r'[CX4]1[OX2][CX4]1',
 'T05_methylenedioxy': r'[OX2][CH2][OX2]', 'S01_thioether': r'[#6][#16X2][#6]', 'S02_disulfide': r'[#16X2][#16X2]',
 'C03a_prenyl_a': r'[c][CX3]([CH3])=[CX3]', 'C03b_prenyl_b': r'[c][CH2][CX3]([CH3])=[CX3]',
}
SM = {k: Chem.MolFromSmarts(v) for k, v in SP.items()}


def is_sugar_ring(mol, ring):
    """Sugar-ring test: 5/6-membered, one O, all-aliphatic carbons, >=3 carbons with exocyclic O."""
    if len(ring) not in (5, 6):
        return False
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    os_ = [a for a in atoms if a.GetSymbol() == 'O']
    os_idx = [a.GetIdx() for a in atoms if a.GetSymbol() == 'O']
    cs_ = [a.GetIdx() for a in atoms if a.GetSymbol() == 'C']
    if len(os_) != 1 or len(cs_) != len(ring) - 1:
        return False
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in cs_):
        return False
    # v2.5 (lessons #43 + #17g): a true sugar ring has >=4 ring carbons carrying exocyclic O
    # (anomeric glycosidic O or hydroxyl). Excludes pantolactone-type lactone rings (3, ring O
    # counted twice) and quinic-acid cyclitols (2). 2-deoxy sugars = blind spot 17c, fail-safe excluded.
    o_ring = os_idx[0]
    n_exo_o = sum(1 for ci in cs_ if any(nb.GetSymbol() == 'O' and nb.GetIdx() != o_ring
                                         for nb in mol.GetAtomWithIdx(ci).GetNeighbors()))
    if n_exo_o < 3:
        return False
    # round-12 sharpening: a true sugar ring carries a -CH2-O- side arm (hydroxymethyl C5-C6);
    # distinguishes bis-tetrahydrofuran lignan rings (no CH2O side arm)
    for ci in cs_:
        c_atom = mol.GetAtomWithIdx(ci)
        for nb in c_atom.GetNeighbors():
            if nb.GetIdx() in os_idx:
                continue
            if nb.GetSymbol() == 'C' and nb.GetTotalNumHs() == 2 and \
               any(x.GetSymbol() == 'O' for x in nb.GetNeighbors()):
                return True
    return False

def scan_nongly(mol):
    """Return the set of non-glycosidic bond codes (context-aware typing for ethers/alkenes; program-level aryl C-C)."""
    hits = set()
    # --- Ester family split (E01/E02/E03/E09; E04 aromatic lactone = fused aromatic ring) ---
    ester_p = Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#8X2H0]')
    lactone_p = Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#8X2;R]')
    galloyl_p = SM['E03_galloyl']
    if mol.HasSubstructMatch(ester_p):
        ri = mol.GetRingInfo()
        is_lactone = False
        lactone_arom = False
        for match in mol.GetSubstructMatches(lactone_p):
            o_idx = match[2] if len(match) >= 3 else match[1]
            if any(is_sugar_ring(mol, r) for r in ri.AtomRings() if o_idx in r):
                continue  # sugar-hydroxyl esterification, not a lactone
            for ring in ri.AtomRings():
                if o_idx in ring:
                    is_lactone = True
                    rs = set(ring)
                    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                        lactone_arom = True
                    else:
                        for ring2 in ri.AtomRings():
                            if len(rs & set(ring2)) >= 2 and any(
                                    mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring2):
                                lactone_arom = True
                    break
            if is_lactone:
                break
        if is_lactone:
            hits.add('E04' if lactone_arom else 'E09')
        else:
            def _aryl_within2(mol, c_idx):
                # aryl reachable within 2 bonds of the carbonyl C (cinnamoyl/benzoyl type)
                for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors():
                    if nb.GetSymbol() != 'C':
                        continue
                    if nb.GetIsAromatic():
                        return True
                    for nb2 in nb.GetNeighbors():
                        if nb2.GetSymbol() == 'C' and nb2.GetIsAromatic():
                            return True
                return False
            aryl_acyl = any(nb.GetIsAromatic() or
                            (nb.GetHybridization() != Chem.HybridizationType.SP3 and _aryl_within2(mol, nb.GetIdx()))
                            for match in mol.GetSubstructMatches(ester_p)
                            for nb in mol.GetAtomWithIdx(match[0]).GetNeighbors()
                            if nb.GetSymbol() == 'C' and nb.GetIdx() != match[1])
            if aryl_acyl:
                hits.add('E02')
                if mol.HasSubstructMatch(galloyl_p):
                    hits.add('E03')
            else:
                hits.add('E01')
    # --- simple SMARTS codes ---
    simple = {'E06': 'E06_phosphate', 'E07': 'E07_sulfate', 'E08': 'E08_thioester', 'N03': 'N03_nitrile',
              'A05': 'A05_alkyne', 'T01': 'T01_methoxy', 'T04': 'T04_epoxide', 'T05': 'T05_methylenedioxy',
              'S01': 'S01_thioether', 'S02': 'S02_disulfide'}
    for code, patt in simple.items():
        if SM[patt] is not None and mol.HasSubstructMatch(SM[patt]): hits.add(code)
    # T01 fix (round 12): ester/acid O-CH3 is not a methoxy ether (O must not touch carbonyl C)
    if 'T01' in hits:
        ester_ome = Chem.MolFromSmarts(r'[CX3](=[OX1])[OX2][CH3]')
        if mol.HasSubstructMatch(ester_ome):
            me_p = SM['T01_methoxy']
            keep = False
            for match in mol.GetSubstructMatches(me_p):
                o_atom = mol.GetAtomWithIdx(match[0])
                if not any(nb.GetSymbol() == 'C' and any(
                        bb.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(nb.GetIdx(), bb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                        for bb in nb.GetNeighbors()) for nb in o_atom.GetNeighbors()):
                    keep = True; break
            if not keep:
                hits.discard('T01')
    # --- Amides: N01/N02/N04 (N02 lactam first, N04 peptide next, else N01) ---
    amide_p = Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#7X3]')
    if mol.HasSubstructMatch(amide_p):
        if mol.HasSubstructMatch(Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#7X3;R]')):
            hits.add('N02')
        elif mol.HasSubstructMatch(SM['N04_peptide']):
            hits.add('N04')
        else:
            hits.add('N01')
    # --- Alkene grading: A01 > A02 > A03 > A04 ---
    cc_p = Chem.MolFromSmarts(r'[CX3]=[CX3]')
    if mol.HasSubstructMatch(cc_p):
        if mol.HasSubstructMatch(SM['A01_activated_alkene']):
            hits.add('A01')
            if mol.HasSubstructMatch(SM['A02_cinnamoyl_alkene']): hits.add('A02')
        elif mol.HasSubstructMatch(SM['A03_conjugated_diene']):
            hits.add('A03')
        else:
            hits.add('A04')
    # --- Ether family: T01 methoxy / T02 alkyl-cyclic / T03 aryl ---
    # exclusion set: sugar-ring O + glycosidic bridge O (anomeric sp3 C-O-aryl);
    # these belong to the glycoside context, not standalone ethers
    sugar_others = set()
    glycosidic_others = set()
    for r in mol.GetRingInfo().AtomRings():
        if is_sugar_ring(mol, r):
            for i in r:
                if mol.GetAtomWithIdx(i).GetSymbol() == 'O':
                    sugar_others.add(i)
    for o in mol.GetAtoms():
        if o.GetSymbol() != 'O' or o.GetIdx() in sugar_others:
            continue
        nbs = list(o.GetNeighbors())
        if len(nbs) != 2:
            continue
        c_nbs = [nb for nb in nbs if nb.GetSymbol() == 'C']
        if len(c_nbs) != 2:
            continue
        # glycosidic-bridge O criteria (v13): ring test narrowed to true sugar rings
        # sugar-sugar bridge: O bonded to two sp3 carbons both in sugar rings;
        # aryl-glycoside bridge: O bonded to aromatic C + sp3 C in a sugar ring
        ri_ = mol.GetRingInfo()

        def _in_sugaring(ci):
            return any(is_sugar_ring(mol, r) for r in ri_.AtomRings() if ci in r)

        both_sp3_in_sugar = all(nb.GetHybridization() == Chem.HybridizationType.SP3 and _in_sugaring(nb.GetIdx())
                                for nb in c_nbs)
        aromC = [nb for nb in c_nbs if nb.GetIsAromatic()]
        sp3Cs = [nb for nb in c_nbs if nb.GetHybridization() == Chem.HybridizationType.SP3]
        arom_plus_sugar = bool(aromC) and len(sp3Cs) == 1 and _in_sugaring(sp3Cs[0].GetIdx())
        if both_sp3_in_sugar or arom_plus_sugar:
            glycosidic_others.add(o.GetIdx())
    alkyl_ether = Chem.MolFromSmarts(r'[CX4][OX2][CX4]')
    ri = mol.GetRingInfo()
    for match in mol.GetSubstructMatches(alkyl_ether):
        o_idx = match[1]
        if o_idx in sugar_others or o_idx in glycosidic_others:
            continue
        o_atom = mol.GetAtomWithIdx(o_idx)
        if any(nb.GetSymbol() == 'C' and not nb.GetIsAromatic() and nb.GetTotalNumHs() == 3
               for nb in o_atom.GetNeighbors()):
            continue  # methoxy -> T01
        if any(len(r) == 3 for r in ri.AtomRings() if o_idx in r):
            continue  # epoxide -> T04
        hits.add('T02')
        break
    aryl_d = Chem.MolFromSmarts(r'c[OX2]c')
    aryl_a = Chem.MolFromSmarts(r'c[OX2][CX4]')
    # round 14: aromatic O of a T05 methylenedioxy bridge belongs to T05, not double-counted as T03
    t05_os = set()
    for match in mol.GetSubstructMatches(SM['T05_methylenedioxy']):
        t05_os.update(match)
    if mol.HasSubstructMatch(aryl_d):
        for match in mol.GetSubstructMatches(aryl_d):
            if match[1] not in sugar_others and match[1] not in glycosidic_others and match[1] not in t05_os:
                hits.add('T03')
                break
    elif mol.HasSubstructMatch(aryl_a):
        for match in mol.GetSubstructMatches(aryl_a):
            o_idx = match[1]
            if o_idx in sugar_others or o_idx in glycosidic_others or o_idx in t05_os:
                continue
            o_atom = mol.GetAtomWithIdx(o_idx)
            if any(nb.GetSymbol() == 'C' and not nb.GetIsAromatic() and nb.GetTotalNumHs() == 3
                   for nb in o_atom.GetNeighbors()):
                continue  # methoxy -> T01
            hits.add('T03')
            break
    # --- Aryl C-C: C01 biphenyl / C02 lignan (C03 prenyl SMARTS) ---
    ri = mol.GetRingInfo()
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetIsAromatic() and a2.GetIsAromatic() and not ri.NumBondRings(bond.GetIdx()):
            # round 13: different SSSR rings; round 15: both atoms in all-carbon six-membered
            # aromatic rings, not adjacent to a ring O -- chromenone/benzopyran perceived-aromatic
            # ring-to-aryl bonds are not cleavable biphenyls
            def _benzene_ok(a):
                for r in ri.AtomRings():
                    if a.GetIdx() in r and len(r) == 6:
                        ratoms = [mol.GetAtomWithIdx(i) for i in r]
                        if all(x.GetSymbol() == 'C' and x.GetIsAromatic() for x in ratoms):
                            ring_o = [x for x in ratoms if any(
                                n.GetSymbol() == 'O' for n in x.GetNeighbors())]
                            # double guard: the carbon must not be adjacent to an aromatic O
                            return True
                return False
            r1 = [set(r) for r in ri.AtomRings() if a1.GetIdx() in r]
            r2 = [set(r) for r in ri.AtomRings() if a2.GetIdx() in r]
            adj_o = False
            for a in (a1, a2):
                for nb in a.GetNeighbors():
                    if nb.GetSymbol() == 'O' and nb.GetIsAromatic():
                        adj_o = True
            if r1 and r2 and not any(x & y for x in r1 for y in r2) \
                    and _benzene_ok(a1) and _benzene_ok(a2) and not adj_o:
                hits.add('C01')
                break
    if mol.HasSubstructMatch(Chem.MolFromSmarts(r'[c][CX4;H1,H2][CX4;H1,H2][c]')):
        hits.add('C02')
    if mol.HasSubstructMatch(SM['C03a_prenyl_a']) or mol.HasSubstructMatch(SM['C03b_prenyl_b']):
        hits.add('C03')
    return hits

class BondJudge:
    """Judge metabolite SMILES against the 46-code bond catalog.

    judge(smiles) -> JudgeResult(bond_codes, glycosidic_linkage,
                                 n_sugar_rings, status)
    """

    @staticmethod
    def judge(smiles):
        hits, g, status = _judge(smiles)
        g = g or {}
        return JudgeResult(
            bond_codes=set(hits),
            glycosidic_linkage=g.get('link'),
            n_sugar_rings=g.get('n_rings', 0),
            status=status,
        )


@dataclasses.dataclass
class JudgeResult:
    bond_codes: set
    glycosidic_linkage: str
    n_sugar_rings: int
    status: str


def _judge(smiles):
    """Single-compound entry point: returns (bond-code set, glycoside info, status)."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles.strip() else None
    if mol is None:
        return set(), None, 'parse_failed'
    g = detect_glycoside(mol)
    hits = scan_nongly(mol)
    # glycoside linkage -> G codes
    if g['link'] == 'O': hits.add('G14')
    elif g['link'] == 'C': hits.add('G16')
    elif g['link'] == 'N': hits.add('G17')
    elif g['link'] == 'S': hits.add('G18')
    elif g['link'] == 'ester': hits.add('E05')
    # 'free'/None are not counted as glycosides
    return hits, g, 'ok'

def run_table(path, out_path, smiles_col='smiles', name_col='name'):
    """Annotate every row of a metabolite table with its bond-type profile.

    Parameters
    ----------
    path : str  Input table (CSV/TSV; encoding auto-detected).
    out_path : str  Output TSV (one row per compound).
    smiles_col, name_col : str  Column names for structure and compound name.

    Returns (DataFrame, Counter, n_ok, n_parse_failed).
    """
    from collections import Counter
    df = None
    for enc in ('utf-8-sig', 'gbk', 'utf-8'):
        try:
            df = pd.read_csv(path, encoding=enc, low_memory=False)
            break
        except Exception:
            continue
    if df is None:
        raise IOError(f'cannot read {path}')
    rows = []
    hit_cnt = Counter()
    n_ok = n_fail = 0
    for _, row in df.iterrows():
        smi = row.get(smiles_col)
        hits, g, status = _judge(smi)
        if status != 'ok':
            n_fail += 1
        else:
            n_ok += 1
        for h in hits:
            hit_cnt[h] += 1
        rows.append({
            'compound_id': row.get('number', row.get('MS2name', '')),
            'name': str(row.get(name_col, ''))[:60],
            'superclass': str(row.get('superclass', ''))[:28],
            'bond_codes': ';'.join(sorted(hits)) if hits else '',
            'n_bond_types': len(hits),
            'glycosidic_linkage': (g['link'] or '') if g else '',
            'n_sugar_rings': (g['n_rings'] if g else 0),
            'status': status,
        })
    out = pd.DataFrame(rows)
    out.to_csv(out_path, sep='	', index=False, encoding='utf-8-sig')
    return out, hit_cnt, n_ok, n_fail


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        prog='bond-scan',
        description='BOND: annotate metabolite SMILES with 46 bond-type codes')
    ap.add_argument('input', help='input table (CSV/TSV with a SMILES column)')
    ap.add_argument('output', help='output TSV (per-compound bond profile)')
    ap.add_argument('--smiles', default='smiles', help='SMILES column name')
    ap.add_argument('--name', default='name', help='compound name column')
    a = ap.parse_args(argv)
    out, cnt, n_ok, n_fail = run_table(a.input, a.output, a.smiles, a.name)
    print(f'annotated {n_ok} rows, parse failed {n_fail} -> {a.output}')
    print('bond-code counts:', dict(sorted(cnt.items(), key=lambda x: -x[1])))


if __name__ == '__main__':
    main()
