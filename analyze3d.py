"""Wearable-availability sweep and per-type event detection."""
import json, pickle, numpy as np
import analyze3 as A

out = {}
base_R = A.R
for tag, lab in (('w0', 'no absent patches'), (None, '15% absent (main)'), ('w30', '30% absent'), ('w50', '50% absent'), ('wfreq', '15% absent, outages every 60 s')):
    R = base_R if tag is None else pickle.load(open(f'results/wards_{tag}.pkl', 'rb'))
    for r in R:
        r['idx'] = int(r['f'][-7:-4])
    A.CAL = [r for r in R if r['idx'] % 2 == 0]; A.TEST = [r for r in R if r['idx'] % 2 == 1]
    tc, tw = A.tune(A.CAL)
    f = A.fusion(A.TEST, tc, tw)
    al, _ = A.evaluate_alarms('fused', 0.01)
    out[lab] = dict(fusion=f, sens=al['pooled']['sens'], far_max=al['pooled']['far_max'])
    print(lab, {k: (round(v['coverage'], 3), round(v['within2'], 3)) for k, v in f.items()}, round(al['pooled']['sens'], 3))
json.dump(out, open('results/wear3.json', 'w'), indent=1)
