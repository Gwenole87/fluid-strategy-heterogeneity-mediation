# Subgroup interactions of the heterogeneity analysis:
#   1. the separate logistic models of Table E2, with 95% confidence intervals,
#   2. one joint model per outcome with the three interactions, Holm correction and a global
#      likelihood-ratio test,
#   3. the overlap between the three subgroups,
#   4. the power of the interaction tests: minimal detectable interaction from the observed standard
#      errors, and 1,000 simulated trials per scenario in which the strategy changes the probability
#      of the outcome within one subgroup and in the opposite direction outside it, with no overall
#      effect (Methods E1).
# Parts 2 and 3 and the minimal detectable interaction of part 4 answer questions of the reviewers:
# their numbers are in the response letter, not in the paper.
# The subgroups are defined as in the Figure 2 cell of CLOVERS_HTE.ipynb. Reads data/DATASET.csv and
# writes aggregate tables to output/subgroups/.
import os, json, warnings
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats
warnings.filterwarnings('ignore')
SEED = 20261027
OUT = os.path.join("output", "subgroups")
os.makedirs(OUT, exist_ok=True)
dataset_path = os.path.join("data", "DATASET.csv")
cols = ['id', 'intervalname', 'term_wd', 'rand_trt', 'term_d90stat', 'term_disch', 'term_dischdt',
        'term_deathdt', 'blvaso_yn', 'labs_lactate', 'charl_kidney', 'medhx_cdialysisyn', 'ardshx_randyn']
raw = pd.read_csv(dataset_path, usecols=cols, low_memory=False)
raw = raw[~raw['id'].isin(raw.loc[raw['term_wd'] == 'Yes', 'id'].unique())]
b = raw[raw['intervalname'] == 'Baseline'].groupby('id').first()
trt = b['rand_trt']
d = pd.DataFrame(index=b.index)
d['restrictive_fluid'] = ((trt == 2) if pd.api.types.is_numeric_dtype(trt) else (trt == 'Restrictive Fluid Group')).astype(int)
d['mortality_90d'] = (b['term_d90stat'] == 'Dead').astype(int)
disch = pd.to_numeric(b['term_dischdt'], errors='coerce'); death = pd.to_numeric(b['term_deathdt'], errors='coerce')
ad = ((b['term_disch'] == 'Yes') & (disch <= 28) & (death.isna() | (death > 28))).astype(float)
ad[b['term_disch'].isna()] = np.nan
d['alive_discharged_28d'] = ad
d['sg_shock'] = ((b['blvaso_yn'] == 'Yes') & (pd.to_numeric(b['labs_lactate'], errors='coerce') > 2)).astype(int)
d['sg_renal'] = (b['charl_kidney'].isin(['Present', 'Present with dialysis']) | (b['medhx_cdialysisyn'] == 'Yes')).astype(int)
d['sg_ards'] = (b['ardshx_randyn'] == 'Yes').astype(int)
SG = [('Septic shock', 'sg_shock'), ('Pre-existing renal impairment', 'sg_renal'), ('ARDS at randomisation', 'sg_ards')]
OUTC = [('primary', 'alive_discharged_28d'), ('mortality90', 'mortality_90d')]
res = {'n': int(len(d)), 'n_restrictive': int(d.restrictive_fluid.sum()), 'outcome_missing': {y: int(d[y].isna().sum()) for _, y in OUTC}}

# 1. Separate models (Table E2) with 95% CI
sep = []
for ck, y in OUTC:
    block = []
    for lab, sg in SG:
        dd = d[['restrictive_fluid', y, sg]].dropna()
        f = smf.logit(f"{y} ~ restrictive_fluid * {sg}", data=dd).fit(disp=0)
        t = f'restrictive_fluid:{sg}'
        lo, hi = f.conf_int().loc[t]
        block.append(dict(outcome=ck, subgroup=lab, model='separate', coef=f.params[t], se=f.bse[t], lo=lo, hi=hi, p=f.pvalues[t],
                          n_sub=int(dd[sg].sum()), n_ref=int((dd[sg] == 0).sum())))
    hp = multipletests([r['p'] for r in block], method='holm')[1]
    for r, h in zip(block, hp): r['p_holm'] = h
    sep += block
sep = pd.DataFrame(sep)

# 2. Joint model per outcome, Holm across its three interactions, global likelihood-ratio test
joint, glob = [], []
for ck, y in OUTC:
    dd = d[['restrictive_fluid', y, 'sg_shock', 'sg_renal', 'sg_ards']].dropna()
    f = smf.logit(f"{y} ~ restrictive_fluid * (sg_shock + sg_renal + sg_ards)", data=dd).fit(disp=0)
    f0 = smf.logit(f"{y} ~ restrictive_fluid + sg_shock + sg_renal + sg_ards", data=dd).fit(disp=0)
    lr = 2 * (f.llf - f0.llf)
    glob.append(dict(outcome=ck, lr_chi2=lr, df=3, p=stats.chi2.sf(lr, 3), n=int(len(dd))))
    block = []
    for lab, sg in SG:
        t = f'restrictive_fluid:{sg}'
        lo, hi = f.conf_int().loc[t]
        block.append(dict(outcome=ck, subgroup=lab, model='joint', coef=f.params[t], se=f.bse[t], lo=lo, hi=hi, p=f.pvalues[t]))
    hp = multipletests([r['p'] for r in block], method='holm')[1]
    for r, h in zip(block, hp): r['p_holm'] = h
    joint += block
joint = pd.DataFrame(joint); glob = pd.DataFrame(glob)

# 3. Overlap of the three subgroups
S = d[['sg_shock', 'sg_renal', 'sg_ards']]
combo = S.astype(str).agg(''.join, axis=1).value_counts().sort_index()
pairs = []
for (la, a), (lb, bb) in [((SG[0][0], 'sg_shock'), (SG[1][0], 'sg_renal')), ((SG[0][0], 'sg_shock'), (SG[2][0], 'sg_ards')), ((SG[1][0], 'sg_renal'), (SG[2][0], 'sg_ards'))]:
    both = int(((S[a] == 1) & (S[bb] == 1)).sum())
    phi = float(np.corrcoef(S[a], S[bb])[0, 1])
    _, pf = stats.fisher_exact(pd.crosstab(S[a], S[bb]).values)
    pairs.append(dict(pair=f"{la} / {lb}", n_both=both, phi=phi, fisher_p=pf))
pairs = pd.DataFrame(pairs)
n_any = int((S.sum(axis=1) >= 1).sum()); n_two_plus = int((S.sum(axis=1) >= 2).sum()); n_three = int((S.sum(axis=1) == 3).sum())

# 4a. Minimal detectable interaction at 80% power (two-sided), from the observed standard errors
z80 = stats.norm.ppf(0.80)
mde = []
for _, r in sep.iterrows():
    for lab_a, alpha in [('unadjusted, alpha 0.05', 0.05), ('Holm first step, alpha 0.05/3', 0.05 / 3)]:
        k = (stats.norm.ppf(1 - alpha / 2) + z80) * r['se']
        mde.append(dict(outcome=r['outcome'], subgroup=r['subgroup'], alpha=lab_a, se=r['se'], mde_logOR=k, mde_ratio_of_OR=float(np.exp(k)), mde_ratio_of_OR_inverse=float(np.exp(-k))))
mde = pd.DataFrame(mde)

# 4b. Simulation: effect in the subgroup offset by the opposite effect outside it, average effect null
rng = np.random.default_rng(SEED)
NSIM = 1000
Smat = S.values
base = {}
for ck, y in OUTC:
    dd = d.dropna(subset=[y])
    for lab, sg in SG:
        lib = dd[dd.restrictive_fluid == 0]
        base[(ck, sg)] = (float(lib.loc[lib[sg] == 1, y].mean()), float(lib.loc[lib[sg] == 0, y].mean()), float(dd[sg].mean()))
def fit_p(y, a, s):
    X = np.column_stack([np.ones_like(a), a, s, a * s]).astype(float)
    try:
        m = sm.Logit(y, X).fit(disp=0, maxiter=100)
        return float(m.pvalues[3]) if np.isfinite(m.pvalues[3]) else 1.0
    except Exception:
        return 1.0
sim = []
for ck, y in OUTC:
    for k, (lab, sg) in enumerate(SG):
        p1, p0, pi = base[(ck, sg)]
        for delta in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
            # direction chosen so that the subgroup moves away from the nearer bound
            sgn = 1.0 if p1 < 0.5 else -1.0
            q1 = p1 + sgn * delta; q0 = p0 - sgn * delta * pi / (1 - pi)
            if not (0.01 < q1 < 0.99 and 0.01 < q0 < 0.99):
                continue
            true_int = (np.log(q1 / (1 - q1)) - np.log(p1 / (1 - p1))) - (np.log(q0 / (1 - q0)) - np.log(p0 / (1 - p0)))
            rej, rej_h = 0, 0
            for _ in range(NSIM):
                idx = rng.integers(0, len(Smat), len(Smat)); Ss = Smat[idx]
                a = rng.permutation(np.r_[np.ones(len(Ss) // 2), np.zeros(len(Ss) - len(Ss) // 2)])
                s_t = Ss[:, k]
                prob = np.where(s_t == 1, np.where(a == 1, q1, p1), np.where(a == 1, q0, p0))
                # the other two subgroups keep the same risk as their stratum of the targeted one (no modification)
                yy = (rng.random(len(Ss)) < prob).astype(float)
                ps = [fit_p(yy, a, Ss[:, j]) for j in range(3)]
                rej += ps[k] < 0.05
                rej_h += multipletests(ps, method='holm')[1][k] < 0.05
            sim.append(dict(outcome=ck, subgroup=lab, prevalence=pi, risk_sub_liberal=p1, risk_rest_liberal=p0,
                            delta_subgroup_pp=100 * sgn * delta, delta_rest_pp=-100 * sgn * delta * pi / (1 - pi),
                            true_interaction_ratio_of_OR=float(np.exp(true_int)), power_unadjusted=rej / NSIM, power_holm=rej_h / NSIM, nsim=NSIM))
sim = pd.DataFrame(sim)
for nm, df_ in [('separate_models', sep), ('joint_models', joint), ('joint_global_lr', glob), ('overlap_pairs', pairs), ('mde', mde), ('simulation_balanced_modifier', sim)]:
    df_.to_csv(os.path.join(OUT, f'{nm}.csv'), index=False)
combo.rename_axis('shock_renal_ards').reset_index(name='n').to_csv(os.path.join(OUT, 'overlap_combinations.csv'), index=False)
res.update(dict(n_any=n_any, n_two_plus=n_two_plus, n_three=n_three, seed=SEED, nsim=NSIM))
json.dump(res, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)
pd.set_option('display.width', 250)
print(json.dumps(res))
print(sep.round(3).to_string(index=False)); print(joint.round(3).to_string(index=False)); print(glob.round(3).to_string(index=False))
print(combo.to_string()); print(pairs.round(3).to_string(index=False))
print(mde.round(3).to_string(index=False))
print(sim.round(3).to_string(index=False))
