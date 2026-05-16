# -*- coding: utf-8 -*-
"""
hs-CRP缺失敏感性分析三合一：
1. MAR检验：基线特征对比
2. IPW：逆概率加权校正
3. Tipping point：极端偏移情景测试
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from scipy import stats
from scipy.stats import chi2_contingency
import statsmodels.api as sm
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Load data
# ============================================================
df = pd.read_csv(r"C:\Users\Ws\Desktop\ai_drug_monitor\crosstalk_results\panel_data.csv")
n_total = df['id'].nunique()
print(f"Loaded: {df.shape[0]} rows, {n_total} persons")

# Mark if person has any CRP
crp_status = df.groupby('id')['Inflammation'].apply(lambda x: x.notna().any()).reset_index()
crp_status.columns = ['id', 'has_crp']
n_has = crp_status['has_crp'].sum()
n_no = n_total - n_has
pct = n_has / n_total * 100
print(f"Has CRP: {n_has} ({pct:.1f}%)")
print(f"No CRP:  {n_no} ({100-pct:.1f}%)\n")

# Encode sex as binary (1=male)
sex_map = {v: 1 if 'nan' not in str(v).lower() and v == list(df['sex'].unique())[0] else 0 for v in df['sex'].unique()}
# Actually explicit: first unique is male based on sample
male_char = df['sex'].iloc[0]  # should be the male char
df['sex_male'] = (df['sex'] == male_char).astype(int)

# Baseline = 2022
baseline = df[df['year'] == 2022].copy()
baseline = baseline.merge(crp_status, on='id')

# ============================================================
# 1. MAR test: baseline comparison
# ============================================================
print("=" * 70)
print("1. MAR Test: Baseline characteristics by hs-CRP availability")
print("=" * 70)

indicators = ['Liver', 'Kidney', 'Metabolic', 'Cardiovascular', 'Lipid']
names = {'Liver': 'ALT', 'Kidney': 'Creatinine', 'Metabolic': 'HbA1c',
         'Cardiovascular': 'SBP', 'Lipid': 'TC'}

def calc_smd(g1, g2):
    m1, m2 = np.nanmean(g1), np.nanmean(g2)
    v1, v2 = np.nanvar(g1), np.nanvar(g2)
    n1, n2 = (~np.isnan(g1)).sum(), (~np.isnan(g2)).sum()
    pv = ((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2) if (n1+n2) > 2 else 1
    return (m1 - m2) / np.sqrt(pv) if pv > 0 else np.nan

has = baseline[baseline['has_crp'] == True]
no = baseline[baseline['has_crp'] == False]

print(f"\n{'Variable':<25} {'CRP+ (n=' + str(len(has)) + ')':>22} {'CRP- (n=' + str(len(no)) + ')':>22} {'SMD':>8} {'p':>10}")
print("-" * 95)

# Age
a_h, a_n = has['age'], no['age']
t, p = stats.ttest_ind(a_h.dropna(), a_n.dropna())
s = calc_smd(a_h.dropna(), a_n.dropna())
sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
print(f"{'Age (years)':<25} {a_h.mean():>8.1f} ({a_h.std():.1f})    {a_n.mean():>8.1f} ({a_n.std():.1f})    {s:>8.3f} {p:>10.4f} {sig:>6}")

# Sex
male_h_pct = has['sex_male'].mean() * 100
male_n_pct = no['sex_male'].mean() * 100
ct = pd.crosstab(baseline['has_crp'], baseline['sex_male'])
chi2, p_sex, _, _ = chi2_contingency(ct)
s_sex = calc_smd(has['sex_male'].values, no['sex_male'].values)
sig_sex = '***' if p_sex<0.001 else '**' if p_sex<0.01 else '*' if p_sex<0.05 else 'ns'
print(f"{'Male (%)':<25} {male_h_pct:>8.1f}%         {male_n_pct:>8.1f}%         {s_sex:>8.3f} {p_sex:>10.4f} {sig_sex:>6}")

# Organ indicators
all_p_vals = [p, p_sex]
all_smds = [abs(s), abs(s_sex)]

for col in indicators:
    name = names[col]
    h_vals = has[col].dropna()
    n_vals = no[col].dropna()
    if len(h_vals) < 10 or len(n_vals) < 10:
        continue
    t, pv = stats.ttest_ind(h_vals, n_vals)
    smd = calc_smd(h_vals, n_vals)
    sig_m = '***' if pv<0.001 else '**' if pv<0.01 else '*' if pv<0.05 else 'ns'
    all_p_vals.append(pv)
    all_smds.append(abs(smd))
    print(f"{name + ' (z-score)':<25} {h_vals.mean():>8.3f} ({h_vals.std():.3f})  {n_vals.mean():>8.3f} ({n_vals.std():.3f})  {smd:>8.3f} {pv:>10.4f} {sig_m:>6}")

# Conclusion
max_smd = max(all_smds) if all_smds else 0
any_sig = any(pv < 0.05 for pv in all_p_vals)

print(f"\n--- MAR Assessment ---")
print(f"Max |SMD| = {max_smd:.3f}")
if max_smd < 0.1 and not any_sig:
    print("CONCLUSION: MCAR supported (missingness unrelated to observables)")
elif any_sig:
    sig_vars = []
    # identify which
    print("CONCLUSION: MAR indicated (missingness associated with baseline)")
    print("  IPW or MI sensitivity analysis warranted")
else:
    print("CONCLUSION: No significant differences, but some SMD > 0.1")

# ============================================================
# 2. IPW
# ============================================================
print("\n" + "=" * 70)
print("2. IPW (Inverse Probability Weighting)")
print("=" * 70)

ipw_data = baseline[['id', 'age', 'sex_male', 'Liver', 'Kidney', 'Metabolic', 'Cardiovascular', 'Lipid', 'has_crp']].copy()

for col in ['Liver', 'Kidney', 'Metabolic', 'Cardiovascular', 'Lipid']:
    ipw_data[col] = ipw_data[col].fillna(ipw_data[col].median())

features = ['age', 'sex_male', 'Liver', 'Kidney', 'Metabolic', 'Cardiovascular', 'Lipid']
X = ipw_data[features].values
y = ipw_data['has_crp'].astype(int).values

lr = LogisticRegression(max_iter=1000, random_state=42)
lr.fit(X, y)
propensity = lr.predict_proba(X)[:, 1]

has_mask = ipw_data['has_crp'].values
ipw_data['ipw_raw'] = np.where(has_mask, 1/propensity, 1/(1-propensity))
q99 = ipw_data['ipw_raw'].quantile(0.99)
ipw_data['ipw_trim'] = ipw_data['ipw_raw'].clip(upper=q99)
ipw_data['ipw_norm'] = ipw_data['ipw_trim'] / ipw_data['ipw_trim'].mean()

print(f"Propensity range: {propensity.min():.4f} - {propensity.max():.4f}")
print(f"IPW raw range:   {ipw_data['ipw_raw'].min():.2f} - {ipw_data['ipw_raw'].max():.2f}")
print(f"IPW trimmed:     {ipw_data['ipw_trim'].min():.2f} - {ipw_data['ipw_trim'].max():.2f}")
print(f"IPW normalized:  {ipw_data['ipw_norm'].min():.2f} - {ipw_data['ipw_norm'].max():.2f}")

# Merge weights
all_data = df.merge(crp_status, on='id')
all_data = all_data.merge(ipw_data[['id', 'ipw_norm']], on='id', how='left')
all_data['ipw_norm'] = all_data['ipw_norm'].fillna(1.0)

# CLPM: Inf -> organs, with and without IPW
targets_info = [
    ('Inflammation', 'Kidney'),
    ('Inflammation', 'Lipid'),
    ('Inflammation', 'Metabolic'),
    ('Inflammation', 'Liver'),
    ('Inflammation', 'Cardiovascular'),
]

print(f"\n{'Edge':<28} {'Orig beta':>10} {'Orig p':>10} {'IPW beta':>10} {'IPW p':>10} {'Delta':>10} {'N pairs':>10}")
print("-" * 82)

for source, target in targets_info:
    pairs = []
    for year in [2022, 2023, 2024]:
        t1 = all_data[all_data['year'] == year]
        t2 = all_data[all_data['year'] == year + 1]
        merged = t1[['id', source, target, 'ipw_norm']].merge(
            t2[['id', target]], on='id', suffixes=('_t', '_t1'))
        merged = merged.dropna(subset=[source, target + '_t', target + '_t1'])
        pairs.append(merged)
    
    panel = pd.concat(pairs, ignore_index=True)
    if len(panel) < 100:
        print(f"{source}->{target:<20} {'N/A':>10}")
        continue
    
    Xm = sm.add_constant(panel[[source, target + '_t']])
    ym = panel[target + '_t1']
    
    m1 = sm.OLS(ym, Xm).fit()
    m2 = sm.WLS(ym, Xm, weights=panel['ipw_norm']).fit()
    
    b1, p1 = m1.params[source], m1.pvalues[source]
    b2, p2 = m2.params[source], m2.pvalues[source]
    d = b2 - b1
    
    edge = f"{source}->{target}"
    print(f"{edge:<28} {b1:>10.5f} {p1:>10.5f} {b2:>10.5f} {p2:>10.5f} {d:>10.5f} {len(panel):>10}")

# Also reverse edges: organs -> Inflammation
print(f"\nReverse: {'Edge':<28} {'Orig beta':>10} {'Orig p':>10} {'IPW beta':>10} {'IPW p':>10} {'Delta':>10} {'N pairs':>10}")
print("-" * 82)

for source in ['Kidney', 'Lipid', 'Metabolic', 'Liver', 'Cardiovascular']:
    pairs = []
    for year in [2022, 2023, 2024]:
        t1 = all_data[all_data['year'] == year]
        t2 = all_data[all_data['year'] == year + 1]
        merged = t1[['id', source, 'Inflammation', 'ipw_norm']].merge(
            t2[['id', 'Inflammation']], on='id', suffixes=('_t', '_t1'))
        merged = merged.dropna(subset=[source, 'Inflammation_t', 'Inflammation_t1'])
        pairs.append(merged)
    
    panel = pd.concat(pairs, ignore_index=True)
    if len(panel) < 100:
        print(f"{source}->Inflammation {'N/A':>10}")
        continue
    
    Xm = sm.add_constant(panel[[source, 'Inflammation_t']])
    ym = panel['Inflammation_t1']
    
    m1 = sm.OLS(ym, Xm).fit()
    m2 = sm.WLS(ym, Xm, weights=panel['ipw_norm']).fit()
    
    b1, p1 = m1.params[source], m1.pvalues[source]
    b2, p2 = m2.params[source], m2.pvalues[source]
    d = b2 - b1
    
    edge = f"{source}->Inflammation"
    print(f"{edge:<28} {b1:>10.5f} {p1:>10.5f} {b2:>10.5f} {p2:>10.5f} {d:>10.5f} {len(panel):>10}")


# ============================================================
# 3. Tipping Point
# ============================================================
print("\n" + "=" * 70)
print("3. Tipping Point Analysis")
print("=" * 70)

offsets = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
core_edges = [('Inflammation', 'Kidney'), ('Inflammation', 'Lipid'), ('Inflammation', 'Metabolic')]

print(f"\n--- Betas across offset scenarios ---")
header = f"{'Edge':<22}"
for o in offsets:
    header += f"  d={o:<4}"
print(header)
print("-" * (22 + 10 * len(offsets)))

for source, target in core_edges:
    row = f"{source}->{target:<12}"
    for offset in offsets:
        df_off = df.copy()
        if offset > 0:
            for pid in df_off['id'].unique():
                pmask = df_off['id'] == pid
                pdata = df_off.loc[pmask, 'Inflammation']
                nonmiss = pdata.notna()
                if nonmiss.sum() == 0 or nonmiss.sum() == 4:
                    continue
                obs_mean = pdata[nonmiss].mean()
                miss_mask = pmask & df_off['Inflammation'].isna()
                df_off.loc[miss_mask, 'Inflammation'] = obs_mean + offset
        
        pairs = []
        for year in [2022, 2023, 2024]:
            t1 = df_off[df_off['year'] == year]
            t2 = df_off[df_off['year'] == year + 1]
            merged = t1[['id', source, target]].merge(
                t2[['id', target]], on='id', suffixes=('_t', '_t1'))
            merged = merged.dropna(subset=[source, target + '_t', target + '_t1'])
            pairs.append(merged)
        
        panel = pd.concat(pairs, ignore_index=True)
        if len(panel) < 100:
            row += f"  {'N/A':>8}"
            continue
        
        Xm = sm.add_constant(panel[[source, target + '_t']])
        ym = panel[target + '_t1']
        m = sm.OLS(ym, Xm).fit()
        b = m.params[source]
        pval = m.pvalues[source]
        sig_mark = '*' if pval < 0.05 else ''
        row += f"  {b:+.4f}{sig_mark:<4}"
    print(row)

# Detailed tipping point
print(f"\n--- Detailed Tipping Point per edge ---")
for source, target in core_edges:
    edge = f"{source}->{target}"
    print(f"\n{edge}:")
    found_tip = False
    for offset in offsets:
        df_off = df.copy()
        if offset > 0:
            for pid in df_off['id'].unique():
                pmask = df_off['id'] == pid
                pdata = df_off.loc[pmask, 'Inflammation']
                nonmiss = pdata.notna()
                if nonmiss.sum() == 0 or nonmiss.sum() == 4:
                    continue
                obs_mean = pdata[nonmiss].mean()
                miss_mask = pmask & df_off['Inflammation'].isna()
                df_off.loc[miss_mask, 'Inflammation'] = obs_mean + offset
        
        pairs = []
        for year in [2022, 2023, 2024]:
            t1 = df_off[df_off['year'] == year]
            t2 = df_off[df_off['year'] == year + 1]
            merged = t1[['id', source, target]].merge(
                t2[['id', target]], on='id', suffixes=('_t', '_t1'))
            merged = merged.dropna(subset=[source, target + '_t', target + '_t1'])
            pairs.append(merged)
        
        panel = pd.concat(pairs, ignore_index=True)
        Xm = sm.add_constant(panel[[source, target + '_t']])
        ym = panel[target + '_t1']
        m = sm.OLS(ym, Xm).fit()
        b = m.params[source]
        pval = m.pvalues[source]
        sig = '***' if pval<0.001 else '**' if pval<0.01 else '*' if pval<0.05 else 'ns'
        print(f"  delta={offset:.1f}SD: beta={b:+.5f}  p={pval:.5f}  [{sig}]  n={len(panel)}")
        
        if not found_tip and pval > 0.05 and offset > 0:
            print(f"  >>> TIPPING POINT: need delta > {offset:.1f} SD to lose significance")
            found_tip = True
    
    if not found_tip:
        print(f"  >>> ROBUST: remains significant even at delta={offsets[-1]}SD")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
