# -*- coding: utf-8 -*-
"""收敛抽样: 第2/3轮, 分层抽20例, 输出供人工核验"""
import sys, random
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import pandas as pd

FAMS = [('糖苷G', 'G'), ('酯E', 'E'), ('烯炔A', 'A'), ('醚T', 'T'),
        ('硫S', 'S'), ('碳碳C', 'C'), ('酰胺腈N', 'N')]

def sample_round(det_path, seed, tag, out):
    det = pd.read_csv(det_path, encoding='utf-8-sig', sep='\t', low_memory=False)
    det2 = det[det['键型数'] > 0]
    picks = []
    per = max(3, 20 // len(FAMS))
    for fam, pref in FAMS:
        sub = det2[det2['检出键型'].astype(str).str.startswith(tuple(f'{pref}{i:02d}' for i in range(1, 10)))]
        if len(sub):
            picks.append(sub.sample(min(per, len(sub)), random_state=seed))
    samp = pd.concat(picks).drop_duplicates(subset=['名称'])
    fec = pd.read_csv(r'D:\01\z_Meta\z\Ori_2024F_4215.csv', encoding='utf-8-sig', low_memory=False)
    smi_map = dict(zip(fec['MS2name'].astype(str).str[:60], fec['smiles']))
    lines = [f'# 抽样轮 {tag} (seed={seed}) n={len(samp)}', '键型\t名称\tSMILES(前85)']
    for _, r in samp.iterrows():
        smi = str(smi_map.get(str(r['名称'])[:60], ''))[:85]
        lines.append(f"{r['检出键型']}\t{str(r['名称'])[:40]}\t{smi}")
    open(out, 'w', encoding='utf-8-sig').write('\n'.join(lines))
    print(f'{tag}: 抽 {len(samp)} 例 -> {out}')

sample_round(r'D:\01\z_software\z_mine\02bond_confirm\02_键型判定_粪便2024_4215.tsv', 202, '第2轮',
             r'D:\01\z_software\z_mine\02bond_confirm\07_抽样核验_第2轮.tsv')
sample_round(r'D:\01\z_software\z_mine\02bond_confirm\02_键型判定_粪便2024_4215.tsv', 303, '第3轮',
             r'D:\01\z_software\z_mine\02bond_confirm\07_抽样核验_第3轮.tsv')
