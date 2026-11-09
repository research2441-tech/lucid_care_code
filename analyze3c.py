"""External checks: (i) stress-test wards (posture dwell 90 s, up to 8 staff visits) with thresholds and detectors
calibrated on the main corpus; (ii) record-disjoint split of the main corpus."""
import json, pickle, numpy as np
import analyze3 as A

out = {}
OOD = pickle.load(open('results/wards_ood.pkl', 'rb'))
for r in OOD:
    r['idx'] = int(r['f'][-7:-4])
tc, tw = A.tune(A.CAL)
out['ood_fusion'] = A.fusion(OOD, tc, tw)
out['main_fusion'] = A.fusion(A.TEST, tc, tw)
TEST0 = A.TEST
A.TEST = OOD
res, _ = A.evaluate_alarms('fused', 0.01)
out['ood_alarms'] = {m: dict(mean=v['far_mean'], mx=v['far_max'], over=v['frac_wards_over'], sens=v['sens']) for m, v in res.items()}
out['ood_events'] = A.events(1.0)
A.TEST = TEST0
# record-disjoint split: records 0-26 calibrate, 27-52 test
for r in A.R:
    r['rec'] = r['idx'] % 53
A.CAL = [r for r in A.R if r['rec'] < 27]; A.TEST = [r for r in A.R if r['rec'] >= 27]
tc2, tw2 = A.tune(A.CAL)
out['disjoint_fusion'] = A.fusion(A.TEST, tc2, tw2)
res, _ = A.evaluate_alarms('fused', 0.01)
out['disjoint_alarms'] = {m: dict(mean=v['far_mean'], mx=v['far_max'], over=v['frac_wards_over'], sens=v['sens']) for m, v in res.items()}
out['disjoint_events'] = A.events(1.0)
json.dump(out, open('results/external3.json', 'w'), indent=1, default=float)
for k, v in out.items():
    print(k, v)
