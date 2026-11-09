"""Fusion, conformal alarm calibration (local / pooled / federated-normalised) and event detection summaries."""
import pickle, json, numpy as np
from scipy.stats import wilcoxon

R = pickle.load(open('results/wards.pkl', 'rb'))
for r in R:
    r['idx'] = int(r['f'][-7:-4])
CAL = [r for r in R if r['idx'] % 2 == 0]
TEST = [r for r in R if r['idx'] % 2 == 1]
WARDS = sorted({r['ward'] for r in R})


def rr_rows(rs):
    a = np.concatenate([np.array(r['rows']) for r in rs if len(r['rows'])])
    return a                                        # t0, ref, rc, qc, rw, qw


def fusion(rs, tc, tw):
    a = rr_rows(rs)
    ref, rc, qc, rw, qw = a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
    okc = np.isfinite(rc) & (qc >= tc); okw = np.isfinite(rw) & (qw >= tw)
    out = {}
    for name, ok, est in (('contactless', okc, rc), ('wearable', okw, rw)):
        e = np.abs(est[ok] - ref[ok])
        out[name] = dict(coverage=float(ok.mean()), within2=float((e < 2).mean()), mae=float(e.mean()))
    pick_c = okc & (~okw | (qc / tc >= qw / tw))
    est = np.where(pick_c, rc, rw); ok = okc | okw
    e = np.abs(est[ok] - ref[ok])
    out['fused'] = dict(coverage=float(ok.mean()), within2=float((e < 2).mean()), mae=float(e.mean()))
    return out


def tune(rs, target=0.9):
    a = rr_rows(rs)
    best = {}
    for k, (est, q) in (('c', (a[:, 2], a[:, 3])), ('w', (a[:, 4], a[:, 5]))):
        cand = []
        for t in np.quantile(q[np.isfinite(est)], np.linspace(0, 0.95, 40)):
            ok = np.isfinite(est) & (q >= t)
            if ok.sum() < 30: continue
            w2 = (np.abs(est[ok] - a[ok, 1]) < 2).mean()
            cand.append((t, w2, ok.mean()))
        good = [c for c in cand if c[1] >= target]
        best[k] = min(good, key=lambda c: c[0])[0] if good else max(cand, key=lambda c: c[1])[0]
    return best['c'], best['w']


def coverage_curve(rs):
    a = rr_rows(rs); out = {}
    for k, (est, q) in (('contactless', (a[:, 2], a[:, 3])), ('wearable', (a[:, 4], a[:, 5]))):
        pts = []
        for t in np.quantile(q[np.isfinite(est)], np.linspace(0, 0.95, 20)):
            ok = np.isfinite(est) & (q >= t)
            pts.append((float(ok.mean()), float((np.abs(est[ok] - a[ok, 1]) < 2).mean())))
        out[k] = pts
    return out


def apnea_windows(rs, mode):
    S, Y, W, E = [], [], [], []
    for r in rs:
        sc, sw, vc, vw = r['sc'], r['sw'], r['vc'], r['vw']
        if mode == 'contactless':
            s, v = sc, vc
        elif mode == 'wearable':
            s, v = sw, vw
        else:
            s = np.where(vc & vw, np.maximum(sc, sw), np.where(vc, sc, sw)); v = vc | vw
        m = v & r['inb']
        lab = np.where(r['ap'] >= 0.8, 1, np.where(r['ap'] == 0, 0, -1))
        S.append(s[m]); Y.append(lab[m]); W.append(np.full(m.sum(), r['ward'])); E.append(np.full(m.sum(), r['idx']))
    return np.concatenate(S), np.concatenate(Y), np.concatenate(W), np.concatenate(E)


def conformal_q(s, alpha):
    n = len(s)
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return np.sort(s)[min(k, n) - 1] if n else np.inf


def calibrate(mode, alpha):
    s, y, w, _ = apnea_windows(CAL, mode)
    neg = y == 0
    th = {}
    th['fixed'] = {wd: 1.0 for wd in WARDS}
    g = conformal_q(s[neg], alpha); th['pooled'] = {wd: g for wd in WARDS}
    th['local'] = {wd: conformal_q(s[neg & (w == wd)], alpha) for wd in WARDS}
    mu = {wd: np.median(s[neg & (w == wd)]) for wd in WARDS}
    sd = {wd: 1.4826 * np.median(np.abs(s[neg & (w == wd)] - mu[wd])) + 1e-9 for wd in WARDS}
    z = np.concatenate([(s[neg & (w == wd)] - mu[wd]) / sd[wd] for wd in WARDS])
    qz = conformal_q(z, alpha)
    th['federated'] = {wd: mu[wd] + sd[wd] * qz for wd in WARDS}
    return th, {wd: int((neg & (w == wd)).sum()) for wd in WARDS}


def evaluate_alarms(mode, alpha):
    if mode == 'fused':
        thc, ncal = calibrate('contactless', alpha); thw, _ = calibrate('wearable', alpha)
        S, Y, W, A = [], [], [], []
        for r in TEST:
            m = (r['vc'] | r['vw']) & r['inb']
            S.append(np.where(r['vw'], r['sw'], r['sc'])[m]); A.append(r['vw'][m])
            Y.append(np.where(r['ap'] >= 0.8, 1, np.where(r['ap'] == 0, 0, -1))[m]); W.append(np.full(m.sum(), r['ward']))
        s, y, w, aw = map(np.concatenate, (S, Y, W, A))
        out = {}
        for meth in thc:
            tt = np.array([thw[meth][x] if a else thc[meth][x] for x, a in zip(w, aw)])
            far = {wd: float(((s > tt) & (y == 0) & (w == wd)).sum() / max(((y == 0) & (w == wd)).sum(), 1)) for wd in WARDS}
            out[meth] = dict(far_mean=float(np.mean(list(far.values()))), far_max=float(np.max(list(far.values()))),
                             frac_wards_over=float(np.mean([v > 2 * alpha for v in far.values()])),
                             sens=float(((s > tt) & (y == 1)).sum() / max((y == 1).sum(), 1)), far=far)
        return out, ncal
    th, ncal = calibrate(mode, alpha)
    s, y, w, ep = apnea_windows(TEST, mode)
    out = {}
    for meth, t in th.items():
        tt = np.array([t[x] for x in w])
        far = {wd: float(((s > tt) & (y == 0) & (w == wd)).sum() / max(((y == 0) & (w == wd)).sum(), 1)) for wd in WARDS}
        sens = float(((s > tt) & (y == 1)).sum() / max((y == 1).sum(), 1))
        out[meth] = dict(far_mean=float(np.mean(list(far.values()))), far_max=float(np.max(list(far.values()))),
                         frac_wards_over=float(np.mean([v > 2 * alpha for v in far.values()])), sens=sens, far=far)
    return out, ncal


def events(fa_per_hour=1.0):
    """Bed-exit / fall detection: gradient boosting on change features (trained on half of the calibration episodes),
    threshold from conformal calibration on the other half for a target false-alarm rate per hour."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    def collect(rs):
        X, L, M, O = [], [], [], []
        for r in rs:
            keep = r['ev_on'] | (r['lab_w'] != 1)
            X.append(r['feats'][keep]); L.append(r['ev_on'][keep]); M.append(r['mov'][keep]); O.append(r['occl'][keep])
        return map(np.concatenate, (X, L, M, O))
    tr = [r for r in CAL if r['idx'] % 4 == 0]; ca = [r for r in CAL if r['idx'] % 4 == 2]
    X, L, _, _ = collect(tr)
    clf = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0).fit(X, L)
    Xc, Lc, _, _ = collect(ca)
    sc = clf.predict_proba(Xc)[:, 1]
    hours = len(ca) * 8 / 60
    alpha = fa_per_hour * hours / max((~Lc).sum(), 1) * 5
    t = conformal_q(sc[~Lc], alpha)
    det, n_ev, fa, fa_mov, fa_occ, lat = 0, 0, 0, 0, 0, []
    for r in TEST:
        p = clf.predict_proba(r['feats'])[:, 1]
        st = np.arange(len(p)) * 2.0
        for k in ('exit', 'fall'):
            if r['ev'][k] is not None:
                n_ev += 1
                on = r['ev'][k] / 50.0
                win = (st <= on + 2.0) & (st + 10.0 >= on)          # windows overlapping this event's onset only
                det += int(((p > t) & win).any())
        al = (p > t) & ~r['ev_on'] & (r['lab_w'] != 1)
        fa += al.sum(); fa_mov += (al & r['mov']).sum(); fa_occ += (al & r['occl']).sum()
    h = len(TEST) * 8 / 60
    return dict(threshold=float(t), events=n_ev, detected=det, sens=det / max(n_ev, 1), fa_per_hour=float(fa / 5 / h),
                fa_share_posture=float(fa_mov / max(fa, 1)), fa_share_occlusion=float(fa_occ / max(fa, 1)))


if __name__ == '__main__':
    tc, tw = tune(CAL)
    out = dict(quality_thresholds=[float(tc), float(tw)], fusion_all=fusion(TEST, -np.inf, -np.inf), fusion_q=fusion(TEST, tc, tw),
               curve=coverage_curve(TEST))
    out['fusion_by_cond'] = {}
    for cond in sorted({(r['d'], r['amb']) for r in TEST}):
        rs = [r for r in TEST if (r['d'], r['amb']) == cond]
        out['fusion_by_cond'][f'{cond[0]}|{cond[1]}'] = fusion(rs, tc, tw)
    out['alarms'] = {}
    for mode in ('contactless', 'wearable', 'fused'):
        for a in (0.01, 0.05):
            res, ncal = evaluate_alarms(mode, a)
            out['alarms'][f'{mode}|{a}'] = res
    out['ncal'] = ncal
    out['events'] = events()
    out['n'] = dict(episodes=len(R), cal=len(CAL), test=len(TEST), wards=len(WARDS), wear_absent=int(sum(r['wear_absent'] for r in R)))
    json.dump(out, open('results/analysis3.json', 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ('curve',)}, indent=1, default=float)[:3500])
