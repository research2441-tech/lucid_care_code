"""Appendix B check of Proposition 2: per-ward total-variation distance between calibration score
distributions and the pooled mixture, against the DKW sampling term of the local quantile."""
import json, numpy as np
import analyze3 as A

def tv(a, b, edges):
    pa = np.histogram(a, edges)[0] / max(len(a), 1); pb = np.histogram(b, edges)[0] / max(len(b), 1)
    return 0.5 * np.abs(pa - pb).sum()

out = {}
for mode in ('contactless', 'wearable'):
    s, y, w, _ = A.apnea_windows(A.CAL, mode)
    neg = y == 0; sn, wn = s[neg], w[neg]
    edges = np.unique(np.quantile(sn, np.linspace(0, 1, 33)))
    edges[0] -= 1e-9; edges[-1] += 1e-9
    rows = {}
    for wd in A.WARDS:
        a = sn[wn == wd]; n = len(a); neff = n / 5.0   # 10 s windows with 2 s hop overlap fivefold
        rows[str(wd)] = dict(n=int(n), tv=float(tv(a, sn, edges)),
                             dkw=float(np.sqrt(np.log(2 / 0.1) / (2 * neff))))
    out[mode] = rows
json.dump(out, open('results/appendix3.json', 'w'), indent=1)
for m, r in out.items():
    tvs = [v['tv'] for v in r.values()]; dk = [v['dkw'] for v in r.values()]
    print(m, 'tv med %.3f max %.3f | dkw med %.3f max %.3f | wards tv<dkw %d/%d' % (
        np.median(tvs), max(tvs), np.median(dk), max(dk), sum(v['tv'] < v['dkw'] for v in r.values()), len(r)))
    for k, v in r.items(): print(' ', k, v)

# Sharper version: only the tail event at the pooled threshold matters.
alpha = 0.01
tail = {}
for mode in ('contactless', 'wearable'):
    s, y, w, _ = A.apnea_windows(A.CAL, mode)
    neg = y == 0; sn, wn = s[neg], w[neg]
    q = A.conformal_q(sn, alpha)
    r = {}
    for wd in A.WARDS:
        a = sn[wn == wd]; neff = len(a) / 5.0
        r[str(wd)] = dict(delta=float(abs((a > q).mean() - alpha)), samp=float(np.sqrt(alpha * (1 - alpha) / neff)))
    tail[mode] = r
    d = [v['delta'] for v in r.values()]; sm = [v['samp'] for v in r.values()]
    print(mode, 'tail delta med %.4f max %.4f | local sd med %.4f max %.4f | wards delta<sd %d/16' % (
        np.median(d), max(d), np.median(sm), max(sm), sum(v['delta'] < v['samp'] for v in r.values())))
out['tail'] = tail
json.dump(out, open('results/appendix3.json', 'w'), indent=1)
