# -*- coding: utf-8 -*-
"""02_accuracy_test.py: legacy交叉核对 + 三表汇总 + 抽样输出"""
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, r'D:\01\z_software\z_mine\02bond_confirm\01script')
import pandas as pd
from collections import Counter
import importlib.util
spec = importlib.util.spec_from_file_location('bj', r'D:\01\z_software\z_mine\02bond_confirm\01script\01_bond_judge.py')
bj = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bj)

BASE = r'D:\01\z_Meta\01biotransformation_result'
Z = r'D:\01\z_Meta\z'
OUT = r'D:\01\z_software\z_mine\02bond_confirm'

# ============ 1. legacy 交叉核对 (侧柏全量) ============
enz = pd.read_csv(f'{BASE}\\02_LableG_new_CBaver_v6_enzyme_integrated.csv', encoding='utf-8-sig', low_memory=False)
det = pd.read_csv(f'{OUT}\\01_键型判定_侧柏2680.tsv', encoding='utf-8-sig', sep='\t', low_memory=False)
det_map = dict(zip(det['编号'], det['检出键型']))
GSET = {'G14', 'G16', 'G17', 'G18', 'E05'}
rows = []
agree = disagree = 0
disagree_examples = []
for _, row in enz.iterrows():
    s = str(row.get('糖苷键型', ''))
    legacy_gly = ('糖苷' in s and '不含' not in s and '游离糖' not in s)
    legacy_legacy_kind = 'legacy苷' if legacy_gly else 'legacy非苷'
    hits = set(str(det_map.get(row['number'], '')).split(';')) - {''}
    det_gly = bool(hits & GSET)
    ok = (legacy_gly == det_gly)
    agree += ok
    disagree += (not ok)
    rows.append((row['number'], str(row['name'])[:40], legacy_legacy_kind,
                 ';'.join(sorted(hits)) or '无', '一致' if ok else '分歧'))
    if not ok and len(disagree_examples) < 40:
        disagree_examples.append((row['number'], str(row['name'])[:36], s[:36],
                                  ';'.join(sorted(hits)) or '无'))
print(f'=== legacy 交叉核对(侧柏2680 全量) ===')
print(f'一致 {agree} ({agree/2680:.1%}) | 分歧 {disagree}')
legacy_gly_total = sum(1 for r in rows if r[2] == 'legacy苷')
det_confirm = sum(1 for r in rows if r[2] == 'legacy苷' and r[4] == '一致')
print(f'legacy苷 {legacy_gly_total} 个中被检测器确认 {det_confirm} ({det_confirm/max(1,legacy_gly_total):.1%})')
pd.DataFrame(rows, columns=['编号', '名称', 'legacy判定', '检测器键型', '一致性']).to_csv(
    f'{OUT}\\06_交叉核对_侧柏.tsv', sep='\t', index=False, encoding='utf-8-sig')
print('\n--- 分歧样例(前12) ---')
for n, name, lk, h in disagree_examples[:12]:
    print(f'  {n} {name} | {lk} | 检测器:{h}')

# ============ 2. 三表键型汇总 ============
summary = []
for tag, path in [('侧柏2680', f'{OUT}\\01_键型判定_侧柏2680.tsv'),
                  ('粪便2024_4215', f'{OUT}\\02_键型判定_粪便2024_4215.tsv'),
                  ('粪便2025_4583', f'{OUT}\\03_键型判定_粪便2025_4583.tsv')]:
    df = pd.read_csv(path, encoding='utf-8-sig', sep='\t', low_memory=False)
    total = len(df)
    hit = (df['键型数'] > 0).sum()
    cnt = Counter()
    for s in df['检出键型'].fillna(''):
        for k in str(s).split(';'):
            if k: cnt[k] += 1
    summary.append((tag, total, hit, f'{hit/total:.1%}', cnt))
all_keys = sorted(set().union(*[set(c.keys()) for _, _, _, _, c in summary]),
                  key=lambda x: (x[0], int(x[1:]) if x[1:].isdigit() else 99))
out_rows = []
for k in all_keys:
    out_rows.append([k] + [c.get(k, 0) for _, _, _, _, c in summary])
out_rows.append(['至少1个键型'] + [hit for _, _, hit, _, _ in summary])
out_rows.append(['化合物总数'] + [t for _, t, _, _, _ in summary])
out_rows.append(['覆盖率'] + [r for _, _, _, r, _ in summary])
pd.DataFrame(out_rows, columns=['键型', '侧柏2680', '粪便2024_4215', '粪便2025_4583']).to_csv(
    f'{OUT}\\04_键型汇总_三表对比.tsv', sep='\t', index=False, encoding='utf-8-sig')
print('\n=== 三表键型汇总(前30行) ===')
print(pd.DataFrame(out_rows[:30], columns=['键型', '侧柏2680', '粪便2024', '粪便2025']).to_string(index=False))
