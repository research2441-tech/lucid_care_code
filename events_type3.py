"""Per-type bed-exit / fall detection with latency, and onset latency of fused apnea alarms (pooled, alpha = 1%)."""
import json, numpy as np
import analyze3 as A
from sklearn.ensemble import HistGradientBoostingClassifier


def collect(rs):
    X, L = [], []
    for r in rs:
        keep = r['ev_on'] | (r['lab_w'] != 1)
        X.append(r['feats'][keep]); L.append(r['ev_on'][keep])
    return np.concatenate(X), np.concatenate(L)


tr = [r for r in A.CAL if r['idx'] % 4 == 0]; ca = [r for r in A.CAL if r['idx'] % 4 == 2]
X, L = collect(tr)
clf = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0).fit(X, L)
Xc, Lc = collect(ca)
sc = clf.predict_proba(Xc)[:, 1]
alpha = 1.0 * len(ca) * 8 / 60 / max((~Lc).sum(), 1) * 5
t = A.conformal_q(sc[~Lc], alpha)
out = {k: dict(n=0, det=0, lat=[]) for k in ('exit', 'fall')}
for r in A.TEST:
    p = clf.predict_proba(r['feats'])[:, 1]
    st = np.arange(len(p)) * 2.0
    for k in ('exit', 'fall'):
        if r['ev'][k] is None:
            continue
        on = r['ev'][k] / 50.0
        win = (st <= on + 2.0) & (st + 10.0 >= on)
        out[k]['n'] += 1
        hit = (p > t) & win
        if hit.any():
            out[k]['det'] += 1; out[k]['lat'].append(float(st[hit][0] + 10.0 - on))
res = {k: dict(n=v['n'], det=v['det'], lat_med=float(np.median(v['lat'])) if v['lat'] else float('nan')) for k, v in out.items()}

# apnea onsets: fused pooled thresholds at alpha = 1%
thc, _ = A.calibrate('contactless', 0.01); thw, _ = A.calibrate('wearable', 0.01)
lat = []; n_on = 0
for r in A.TEST:
    s = np.where(r['vw'], r['sw'], r['sc']); ok = (r['vc'] | r['vw']) & r['inb']
    th = np.where(r['vw'], thw['pooled'][r['ward']], thc['pooled'][r['ward']])
    alarm = ok & (s > th)
    a = r['ap'] > 0
    starts = np.flatnonzero(a & ~np.r_[False, a[:-1]])
    for s0 in starts:
        n_on += 1
        e0 = s0
        while e0 < len(a) and a[e0]:
            e0 += 1
        hit = np.flatnonzero(alarm[s0:e0])
        if len(hit):
            lat.append(2.0 * hit[0])
res['apnea_latency_med'] = float(np.median(lat)); res['apnea_latency_p90'] = float(np.quantile(lat, 0.9))
res['apnea_onsets_detected'] = len(lat); res['apnea_onsets'] = n_on
json.dump(res, open('results/events_type3.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
