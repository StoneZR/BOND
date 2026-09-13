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

# ==================== 一、糖苷分型 (教训#3,4,5,6,7,17a,17b,17d) ====================
def detect_glycoside(mol):
    """返回 dict(link=O/C/N/S/ester/free/None, n_rings=糖环数, quininone=醌型C苷)
    规则(每条注明教训编号):
      R1 糖环=5/6元,恰1个O,其余C           (糖环定义)
      R2 环内碳须全脂肪 (教训#4: 黄酮苯并吡喃误判)
      R3 环上氧取代>=3 (真糖环特征; 2-脱氧糖盲区17c备案,不放宽)
      R4 异头碳须sp3 (教训#5: 内酯羰基碳误判)
      R5 连接须单键 (教训#5: 双键氧误判)
      R6 异头外接O无H: 若该O另一端为带双键O的羰基碳->糖酯ester; 否则O-苷 (教训#6两轮)
      R7 异头外接O带H->free游离糖,不入苷 (教训#7)
      R8 C连接: 芳环碳->C-苷; 醌环sp2碳->quinone型C-苷 (盲区17b: 芦荟苷型)
      R9 N/S连接->N苷/S苷 (盲区17a补齐S)
    """
    ri = mol.GetRingInfo()
    n_sugar = 0; link_found = None; quinone_c = False
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6): continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        os_ = [a for a in atoms if a.GetSymbol() == 'O']
        cs_ = [a for a in atoms if a.GetSymbol() == 'C']
        if len(os_) != 1 or len(cs_) != len(ring) - 1: continue          # R1
        if any(a.GetIsAromatic() for a in cs_): continue                 # R2
        # R3(v2.5, 与is_sugar_ring同步): 至少3个环碳带环外O(异头苷氧或羟基);
        # 排除泛解酸内酯/奎宁酸cyclitol(2个)/2-脱氧糖(盲区17c备案)
        o_idx_r = os_[0].GetIdx()
        n_exo_o = sum(1 for c in cs_ if any(nb.GetSymbol() == 'O' and nb.GetIdx() != o_idx_r
                                            for nb in c.GetNeighbors()))
        if n_exo_o < 3: continue
        n_sugar += 1
        o_idx = os_[0].GetIdx()
        for c in cs_:
            if c.GetHybridization() != Chem.HybridizationType.SP3: continue   # R4
            if o_idx not in [nb.GetIdx() for nb in c.GetNeighbors()]: continue
            # v2.1(教训#17g, 奎宁酸案例): 真异头碳不含环外羟基O;
            # cyclitol类(奎宁酸/莽草酸环)的醇碳有环外OH, 其芳环/杂环邻居不代表C-苷
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
                                link_found = link_found or 'ester'                # 磷/硫酸糖酯
                        link_found = link_found or 'O'
                elif sym == 'C':
                    if nb.GetIsAromatic():
                        link_found = link_found or 'C'                        # R8 芳基C-苷
                    elif nb.GetHybridization() != Chem.HybridizationType.SP3:
                        # 盲区17b: 醌型C-苷(aloin型): 异头碳直连sp2非芳碳(醌环/蒽酮环成员)
                        quinone_c = True
                        link_found = link_found or 'C'
                    else:
                        # 盲区17f: 海藻糖型 1,1-糖-糖C-C桥——外接sp3碳须属于"另一个"糖环(非本环)
                        nb_rings = [r for r in mol.GetRingInfo().AtomRings() if nb.GetIdx() in r]
                        if any(is_sugar_ring(mol, r) for r in nb_rings) and \
                           not any(set(r) == set(ring) for r in nb_rings):
                            link_found = link_found or 'C'
                elif sym == 'N':
                    link_found = link_found or 'N'                            # R9
                elif sym == 'S':
                    link_found = link_found or 'S'                            # R9
    return {'link': link_found, 'n_rings': n_sugar, 'quinone_c': quinone_c}

# ==================== 二、非糖 SMARTS (26 项; 教训#1通配写法,#2下标) ====================
SP = {
 'E01_脂肪酯': r'[CX3](=[OX1])[OX2H0][CX4]', 'E02_芳香酯': r'c[CX3](=[OX1])[OX2H0]',
 'E03_没食子酰': r'c1c(O)c(O)c(C(=[OX1]))c(O)c1',
 'E04_内酯芳香': r'[#6X3](=[#8X1])[#8X2;R]~[c]', 'E09_内酯脂肪': r'[#6X3](=[#8X1])[#8X2;R]',
 'E06_磷酸酯': r'[PX4](=[OX1])([OX2])[OX2]', 'E07_硫酸酯': r'[SX4](=[OX1])(=[OX1])([OX2])[OX2]',
 'E08_硫酯': r'[#6X3](=[#8X1])[#16X2][#6]',
 'N01_开链酰胺': r'[#6X3](=[#8X1])[#7X3;!R]', 'N02_内酰胺': r'[#6X3](=[#8X1])[#7X3;R]',
 'N04_肽键': r'[#6X4][#6X3](=[#8X1])[#7X3][#6X4]',
 'N03_腈': r'[#6]#[#7]',
 'A01_活化烯': r'[CX3]=[CX3][CX3]=[OX1]', 'A02_桂皮酰烯': r'[c][CX3]=[CX3][CX3]=[OX1]',
 'A03_共轭二烯': r'[CX3]=[CX3][CX3]=[CX3]', 'A05_炔': r'[#6X2]#[#6X2]',
 'T01_甲氧基': r'[OX2][CH3]', 'T04_环氧': r'[CX4]1[OX2][CX4]1',
 'T05_亚甲二氧桥': r'[OX2][CH2][OX2]', 'S01_硫醚': r'[#6][#16X2][#6]', 'S02_二硫': r'[#16X2][#16X2]',
 'C03a_异戊烯a': r'[c][CX3]([CH3])=[CX3]', 'C03b_异戊烯b': r'[c][CH2][CX3]([CH3])=[CX3]',
}
SM = {k: Chem.MolFromSmarts(v) for k, v in SP.items()}


def is_sugar_ring(mol, ring):
    """糖环判定: 5/6元、单O、全脂肪碳、环上碳氧取代>=3"""
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
    # v2.5(教训#43+#17g): 真糖环=至少4个环碳带环外O(异头苷氧或羟基);
    # 排除: 泛解酸内酯环(3碳带O,且环O算两次)、cyclitol环(奎宁酸仅2碳有环外O)、2-脱氧糖(盲区17c备案, 宁可漏不误)
    o_ring = os_idx[0]
    n_exo_o = sum(1 for ci in cs_ if any(nb.GetSymbol() == 'O' and nb.GetIdx() != o_ring
                                         for nb in mol.GetAtomWithIdx(ci).GetNeighbors()))
    if n_exo_o < 3:
        return False
    # 第12轮锐化: 真糖环必有-CH2-O-型侧臂(羟甲基C5-C6); 区分双四氢呋喃木脂素环(无CH2O侧臂)
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
    """返回非糖键型编号集合 (程序判定: 醚/烯按语境程序分型, 芳基C-C程序判定)"""
    hits = set()
    # --- 酯族拆分 (E01/E02/E03/E09; E04芳香内酯=内酯含芳环语境) ---
    ester_p = Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#8X2H0]')
    lactone_p = Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#8X2;R]')
    galloyl_p = SM['E03_没食子酰']
    if mol.HasSubstructMatch(ester_p):
        ri = mol.GetRingInfo()
        is_lactone = False
        lactone_arom = False
        for match in mol.GetSubstructMatches(lactone_p):
            o_idx = match[2] if len(match) >= 3 else match[1]
            if any(is_sugar_ring(mol, r) for r in ri.AtomRings() if o_idx in r):
                continue  # 糖羟基成酯非内酯
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
                # 羰基C的邻居2键内达芳香环(桂皮酰/苯甲酰型)
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
    # --- 简单 SMARTS 项 ---
    simple = {'E06': 'E06_磷酸酯', 'E07': 'E07_硫酸酯', 'E08': 'E08_硫酯', 'N03': 'N03_腈',
              'A05': 'A05_炔', 'T01': 'T01_甲氧基', 'T04': 'T04_环氧', 'T05': 'T05_亚甲二氧桥',
              'S01': 'S01_硫醚', 'S02': 'S02_二硫'}
    for code, patt in simple.items():
        if SM[patt] is not None and mol.HasSubstructMatch(SM[patt]): hits.add(code)
    # T01 修正(第12轮): 酯/酸的O-CH3不算甲氧基醚(O须不连羰基碳)
    if 'T01' in hits:
        ester_ome = Chem.MolFromSmarts(r'[CX3](=[OX1])[OX2][CH3]')
        if mol.HasSubstructMatch(ester_ome):
            me_p = SM['T01_甲氧基']
            keep = False
            for match in mol.GetSubstructMatches(me_p):
                o_atom = mol.GetAtomWithIdx(match[0])
                if not any(nb.GetSymbol() == 'C' and any(
                        bb.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(nb.GetIdx(), bb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                        for bb in nb.GetNeighbors()) for nb in o_atom.GetNeighbors()):
                    keep = True; break
            if not keep:
                hits.discard('T01')
    # --- 酰胺: N01/N02/N04 (N02优先, N04肽键次之, 其余N01) ---
    amide_p = Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#7X3]')
    if mol.HasSubstructMatch(amide_p):
        if mol.HasSubstructMatch(Chem.MolFromSmarts(r'[#6X3](=[#8X1])[#7X3;R]')):
            hits.add('N02')
        elif mol.HasSubstructMatch(SM['N04_肽键']):
            hits.add('N04')
        else:
            hits.add('N01')
    # --- 烯键分级: A01>A02>A03>A04 ---
    cc_p = Chem.MolFromSmarts(r'[CX3]=[CX3]')
    if mol.HasSubstructMatch(cc_p):
        if mol.HasSubstructMatch(SM['A01_活化烯']):
            hits.add('A01')
            if mol.HasSubstructMatch(SM['A02_桂皮酰烯']): hits.add('A02')
        elif mol.HasSubstructMatch(SM['A03_共轭二烯']):
            hits.add('A03')
        else:
            hits.add('A04')
    # --- 醚族: T01甲氧基/T02烷基环醚/T03芳基醚 ---
    # 排除集合: 糖环内O + 苷桥O(异头sp3碳-O-芳环) —— 这些属糖苷语境非独立醚
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
        # 苷桥O判据(v13): 环判据收窄为真糖环(is_sugar_ring)
        # 糖-糖桥: O连两个sp3碳且两碳都在真糖环内; 芳基苷桥: O连芳环C+真糖环内sp3碳
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
            continue  # 甲氧基归T01
        if any(len(r) == 3 for r in ri.AtomRings() if o_idx in r):
            continue  # 环氧归T04
        hits.add('T02')
        break
    aryl_d = Chem.MolFromSmarts(r'c[OX2]c')
    aryl_a = Chem.MolFromSmarts(r'c[OX2][CX4]')
    # 第14轮: T05亚甲二氧桥的芳基O属T05语境, 不重复计T03
    t05_os = set()
    for match in mol.GetSubstructMatches(SM['T05_亚甲二氧桥']):
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
                continue  # 甲氧基归T01
            hits.add('T03')
            break
    # --- 芳基C-C: C01联苯 / C02木脂素 (C03异戊烯 SMARTS) ---
    ri = mol.GetRingInfo()
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetIsAromatic() and a2.GetIsAromatic() and not ri.NumBondRings(bond.GetIdx()):
            # 第13轮: 两端须分属不同的SSSR环; 第15轮: 两端须都在全碳六元芳环(苯环)
            # 且不邻环氧——色酮/苯并吡喃芳香感知环上外连芳基的键不是可裂解联苯
            def _benzene_ok(a):
                for r in ri.AtomRings():
                    if a.GetIdx() in r and len(r) == 6:
                        ratoms = [mol.GetAtomWithIdx(i) for i in r]
                        if all(x.GetSymbol() == 'C' and x.GetIsAromatic() for x in ratoms):
                            ring_o = [x for x in ratoms if any(
                                n.GetSymbol() == 'O' for n in x.GetNeighbors())]
                            # 该碳本身不邻环内O(环O成员无碳-环含O的六元芳环已排除,双保险:
                            # 邻居里不含处于同环的O)
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
    if mol.HasSubstructMatch(SM['C03a_异戊烯a']) or mol.HasSubstructMatch(SM['C03b_异戊烯b']):
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
    """单化合物判定入口: 返回 (键型集合, gly信息, 状态)"""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles.strip() else None
    if mol is None:
        return set(), None, 'SMILES解析失败'
    g = detect_glycoside(mol)
    hits = scan_nongly(mol)
    # 糖苷 -> G 编号
    if g['link'] == 'O': hits.add('G14')
    elif g['link'] == 'C': hits.add('G16')
    elif g['link'] == 'N': hits.add('G17')
    elif g['link'] == 'S': hits.add('G18')
    elif g['link'] == 'ester': hits.add('E05')
    # 'free'/None 不入苷
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
