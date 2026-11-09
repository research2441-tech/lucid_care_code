"""Error analysis by reference respiratory-rate band (bradypnoea <12, normal 12-20, tachypnoea >20 bpm)
and agreement statistics (Bland-Altman bias and limits) of accepted windows on the test half."""
import json, numpy as np
import analyze3 as A

tc, tw = json.load(open('results/analysis3.json'))['quality_thresholds']
a = A.rr_rows(A.TEST)
ref, rc, qc, rw, qw = a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
okc = np.isfinite(rc) & (qc >= tc); okw = np.isfinite(rw) & (qw >= tw)
pick_c = okc & (~okw | (qc / tc >= qw / tw))
fz = np.where(pick_c, rc, rw); okf = okc | okw
bands = {'low (<12)': ref < 12, 'normal (12-20)': (ref >= 12) & (ref <= 20), 'high (>20)': ref > 20}
out = {'bands': {}, 'ba': {}}
for b, m in bands.items():
    row = {'n': int(m.sum())}
    for name, ok, est in (('contactless', okc, rc), ('wearable', okw, rw), ('fused', okf, fz)):
        s = ok & m; e = est[s] - ref[s]
        row[name] = dict(coverage=float(s.sum() / max(m.sum(), 1)), within2=float((np.abs(e) < 2).mean()) if s.any() else float('nan'),
                         bias=float(e.mean()) if s.any() else float('nan'))
    out['bands'][b] = row
for name, ok, est in (('contactless', okc, rc), ('wearable', okw, rw), ('fused', okf, fz)):
    e = est[ok] - ref[ok]
    out['ba'][name] = dict(n=int(ok.sum()), bias=float(e.mean()), lo=float(e.mean() - 1.96 * e.std()), hi=float(e.mean() + 1.96 * e.std()),
                           med_abs=float(np.median(np.abs(e))), p95=float(np.quantile(np.abs(e), 0.95)))
# tachypnoea flag (>20 bpm) agreement for fused accepted windows
s = okf
tp = ((fz > 20) & (ref > 20) & s).sum(); fn = ((fz <= 20) & (ref > 20) & s).sum(); fp = ((fz > 20) & (ref <= 20) & s).sum()
out['tachy'] = dict(sens=float(tp / max(tp + fn, 1)), ppv=float(tp / max(tp + fp, 1)), n_pos=int(tp + fn))
json.dump(out, open('results/errband3.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
