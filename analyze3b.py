"""Extra analyses: alpha sweep, calibration-budget sweep, event ROC."""
import json, numpy as np
import analyze3 as A

out = {}
out['alpha'] = {}
for a in (0.005, 0.01, 0.02, 0.05, 0.1):
    res, _ = A.evaluate_alarms('fused', a)
    out['alpha'][str(a)] = {m: dict(mean=v['far_mean'], mx=v['far_max'], over=v['frac_wards_over'], sens=v['sens']) for m, v in res.items()}
out['budget'] = {}
full = list(A.CAL)
for frac in (0.25, 0.5, 1.0):
    vals = {m: [] for m in ('local', 'federated', 'pooled')}
    for rep in range(5):
        rng = np.random.default_rng(rep)
        sub = [r for r in full if rng.random() < frac] if frac < 1 else full
        A.CAL = sub
        res, _ = A.evaluate_alarms('fused', 0.01)
        for m in vals:
            vals[m].append((res[m]['far_max'], res[m]['frac_wards_over'], res[m]['far_mean']))
        if frac == 1.0:
            break
    out['budget'][str(frac)] = {m: [float(np.mean([x[i] for x in v])) for i in range(3)] for m, v in vals.items()}
A.CAL = full
roc = []
for fa in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
    ev = A.events(fa)
    roc.append((fa, ev['fa_per_hour'], ev['sens']))
out['roc'] = roc
json.dump(out, open('results/extra3.json', 'w'), indent=1)
print(json.dumps(out, indent=0)[:2500])
