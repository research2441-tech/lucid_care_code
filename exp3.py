"""LUCID-Care experiments: rate-CRB-accuracy trade-off, events, wearable fusion, federated conformal alarms."""
import sys, glob, pickle, json, numpy as np
from multiprocessing import Pool
import isac, sense

import os
FILES = sorted(glob.glob(os.environ.get('LC_OUT', 'eps') + '/ep_*.pkl'))
FS = sense.FS


def rr_eval(e, Y, nsd, valid):
    W = isac.rr_eigen(Y, valid, nsd)
    rows = []
    for t0, t1, rr, q in W:
        fr = slice(int(t0 * FS), int(t1 * FS))
        if e['apnea'][fr].any() or (e['lab'][fr] != 0).any():
            continue
        ref = isac.rr_reference(e['breaths'], t0, t1)
        if np.isfinite(ref):
            rows.append((t0, ref, rr, q))
    return rows


def pareto_job(a):
    f, d, rp, amb = a
    e = pickle.load(open(f, 'rb')); rng = np.random.default_rng(abs(hash((f, d, rp, amb))) % 2 ** 32)
    Y = isac.estimate(e, d, rp, rng, amb)
    nsd = np.tile(np.sqrt(isac.crb(e['Hdc'], d, rp, amb)), e['H'].shape[1])
    valid = ~isac.artifact_mask(Y) & (e['lab'] == 0)
    return (d, rp, amb), rr_eval(e, Y, nsd, valid)


def pareto():
    grid = [(d, rp, amb) for amb in (5e-5, 5e-6) for d in (0.2, 0.35, 0.5) for rp in (0.01, 0.03, 0.1, 0.2, 0.4)]
    jobs = [(f, *g) for g in grid for f in FILES[::4]]
    out = {g: [] for g in grid}
    with Pool(2) as p:
        for g, rows in p.imap_unordered(pareto_job, jobs, chunksize=4):
            out[g] += rows
    res = {}
    for (d, rp, amb), rows in out.items():
        a = np.array(rows) if rows else np.zeros((0, 4))
        ok = np.isfinite(a[:, 2]) if len(a) else np.zeros(0, bool)
        err = np.abs(a[ok, 2] - a[ok, 1]) if len(a) else np.zeros(0)
        e0 = pickle.load(open(FILES[0], 'rb'))
        res[f'{d}|{rp}|{amb}'] = dict(d=d, rp=rp, amb=amb, rate=float(isac.rate(d, rp)), crb=float(np.mean(isac.crb(e0['Hdc'], d, rp, amb))),
                                     mae=float(err.mean()) if len(err) else np.nan, within2=float((err < 2).mean()) if len(err) else np.nan,
                                     coverage=float(ok.mean()) if len(a) else 0.0, n=int(len(a)))
        print(res[f'{d}|{rp}|{amb}'], flush=True)
    json.dump(res, open('results/pareto.json', 'w'), indent=1)



CONDS = [(0.2, 5e-6), (0.35, 5e-6), (0.5, 5e-6), (0.35, 5e-5), (0.5, 5e-5)]       # (dimming, ambient) per ward
RHO = 0.1


def ward_cond(layout):
    return CONDS[layout % len(CONDS)]


def band_power(sig_mat, nsd, win, hop):
    """Coherent respiratory-band power per window: mean over band of the largest cross-spectral eigenvalue."""
    fs = FS; n, h = int(win * fs), int(hop * fs)
    Z = isac.sosfiltfilt(isac.SOS, sig_mat - sig_mat.mean(0), axis=0) / (nsd + 1e-30)
    out = []
    for s in range(0, len(Z) - n + 1, h):
        X = Z[s:s + n]
        C = X.T @ X / n
        out.append(np.linalg.eigvalsh(C)[-1] if X.shape[1] > 1 else C[0, 0])
    return np.array(out)


def ward_job(f):
    e = pickle.load(open(f, 'rb'))
    rng = np.random.default_rng(int(f[-7:-4]) + 99)
    d, amb = ward_cond(e['layout'])
    Y = isac.estimate(e, d, RHO, rng, amb)
    nsd = np.tile(np.sqrt(isac.crb(e['Hdc'], d, RHO, amb)), e['H'].shape[1])
    art = isac.artifact_mask(Y)
    valid_c = ~art & (e['lab'] == 0)
    # wearable chest patch: impedance respiration with sensor noise; absent episodes and uplink outages
    T = len(e['resp'])
    avail = np.ones(T, bool)
    if rng.random() < float(os.environ.get('LC_ABSENT', 0.15)):
        avail[:] = False
    else:
        t = 0; on = True
        while t < T:
            L = int(rng.exponential(float(os.environ.get('LC_UP', 300.0)) if on else float(os.environ.get('LC_DOWN', 20.0))) * FS) + 1
            avail[t:t + L] = on; t += L; on = not on
    wsig = e['resp'] + rng.standard_normal(T) * 0.1 * e['resp'].std() * np.sqrt(10)
    wsig = np.where(avail, wsig, 0.0)
    Wc = isac.rr_eigen(Y, valid_c, nsd)
    Ww = isac.rr_eigen(wsig[:, None], avail & (e['lab'] == 0), np.ones(1))
    rows = []
    for (t0, t1, rc, qc), (_, _, rw, qw) in zip(Wc, Ww):
        fr = slice(int(t0 * FS), int(t1 * FS))
        if e['apnea'][fr].any() or (e['lab'][fr] != 0).any():
            continue
        ref = isac.rr_reference(e['breaths'], t0, t1)
        if np.isfinite(ref):
            rows.append((t0, ref, rc, qc, rw, qw))
    # apnea scores on 10 s windows, hop 2 s, relative to a causal 120 s baseline
    win, hop = 10.0, 2.0
    pc = band_power(np.where(valid_c[:, None], Y, Y.mean(0)), nsd, win, hop)
    pw = band_power(np.where(avail[:, None], wsig[:, None], 0.0), np.ones(1), win, hop)
    nb = int(120 / hop)
    def rel(p):
        base = np.array([np.median(p[max(0, i - nb):i]) if i >= 5 else np.median(p[:5]) for i in range(len(p))])
        return np.log(base + 1e-12) - np.log(p + 1e-12)
    sc, sw = rel(pc), rel(pw)
    n = int(win * FS); h = int(hop * FS)
    st = np.arange(0, T - n + 1, h)
    ap = np.array([e['apnea'][s:s + n].mean() for s in st]); inb = np.array([(e['lab'][s:s + n] == 0).all() for s in st])
    vc = np.array([valid_c[s:s + n].mean() > 0.95 for s in st]); vw = np.array([avail[s:s + n].mean() > 0.95 for s in st])
    # gross events (bed exit / fall): change of the slow total body-reflection relative to its noise
    tot = isac.sosfiltfilt(isac.SOS_LP, Y.mean(1)); tn = np.sqrt(np.mean(nsd ** 2) / Y.shape[1]) * np.sqrt(0.6 / 25)
    g = np.abs(tot[int(2 * FS):] - tot[:-int(2 * FS)]) / tn
    gscore = np.array([g[s:s + n].max() if s + n <= len(g) else g[s:].max() for s in st])
    Ls = isac.sosfiltfilt(isac.SOS_LP, Y, axis=0)
    w5 = int(5 * FS)
    feats = []
    for s0 in st:
        c = s0 + n // 2
        a = Ls[max(0, c - 2 * w5):c - w5 // 2].mean(0) if c - w5 // 2 > 0 else Ls[:1].mean(0)
        b = Ls[c + w5 // 2:c + 2 * w5].mean(0) if c + w5 // 2 < T else Ls[-1:].mean(0)
        dv = (b - a) / (np.abs(e['Hs']).reshape(-1).mean() * 0 + np.abs(Ls).mean(0) * 0 + nsd * np.sqrt(0.6 / 25) + 1e-30)
        sl = np.abs(np.diff(Ls[max(0, c - w5):c + w5:int(FS / 2)].mean(1))).max() / (np.sqrt(np.mean(nsd ** 2) / Y.shape[1]) + 1e-30)
        feats.append(np.concatenate([[np.log1p(np.linalg.norm(dv)), np.log1p(sl)], np.sign(dv) * np.log1p(np.abs(dv))]))
    feats = np.array(feats, np.float32)
    ev_on = np.zeros(len(st), bool)
    for k in ('exit', 'fall'):
        if e['ev'][k] is not None:
            ev_on |= (st <= e['ev'][k] + 2 * FS) & (st + n >= e['ev'][k])
    return dict(f=f, ward=e['layout'], d=d, amb=amb, rows=rows, sc=sc, sw=sw, ap=ap, inb=inb, vc=vc, vw=vw, gscore=gscore,
                ev_on=ev_on, feats=feats, mov=np.array([e['mov'][s:s + n].any() for s in st]), lab_w=np.array([np.bincount(e['lab'][s:s + n], minlength=4).argmax() for s in st]), occl=np.array([e['occl'][s:s + n].any() for s in st]), ev=e['ev'], wear_absent=not avail.any())


def wards():
    with Pool(2) as p:
        R = list(p.imap_unordered(ward_job, FILES, chunksize=4))
    pickle.dump(R, open(os.environ.get('LC_WARDS', 'results/wards.pkl'), 'wb'))
    print('done', len(R))



def abl_job(a):
    f, var = a
    e = pickle.load(open(f, 'rb')); rng = np.random.default_rng(int(f[-7:-4]) + 7)
    d, rp, amb = 0.5, 0.1, 5e-6
    Y = isac.estimate(e, d, rp, rng, amb)
    nsd = np.tile(np.sqrt(isac.crb(e['Hdc'], d, rp, amb)), e['H'].shape[1])
    valid = (e['lab'] == 0) & (np.ones(len(Y), bool) if var == 'nogate' else ~isac.artifact_mask(Y))
    if var == 'best':
        k = np.argmax(e['H'].reshape(len(Y), -1).mean(0)); Y, nsd = Y[:, [k]], nsd[[k]]
    if var == 'pca':
        W, _ = isac.rr_track(Y, valid, noise_sd=nsd)
    elif var == 'nowhite':
        W = isac.rr_eigen(Y, valid, np.full_like(nsd, nsd.mean()) * (Y.std(0) / Y.std(0).mean()))
    else:
        W = isac.rr_eigen(Y, valid, nsd)
    rows = []
    for t0, t1, rr, q in W:
        fr = slice(int(t0 * FS), int(t1 * FS))
        if e['apnea'][fr].any() or (e['lab'][fr] != 0).any():
            continue
        ref = isac.rr_reference(e['breaths'], t0, t1)
        if np.isfinite(ref):
            rows.append((ref, rr))
    return var, rows


def ablation():
    vars_ = ['full', 'nogate', 'best', 'pca', 'nowhite']
    out = {v: [] for v in vars_}
    with Pool(2) as p:
        for v, rows in p.imap_unordered(abl_job, [(f, v) for v in vars_ for f in FILES[1::4]], chunksize=4):
            out[v] += rows
    res = {}
    for v, rows in out.items():
        a = np.array(rows); ok = np.isfinite(a[:, 1]); err = np.abs(a[ok, 1] - a[ok, 0])
        res[v] = dict(coverage=float(ok.mean()), within2=float((err < 2).mean()), mae=float(err.mean()), n=int(len(a)))
    json.dump(res, open('results/ablation3.json', 'w'), indent=1); print(res)



def design_job(a):
    f, var, val = a
    e = pickle.load(open(f, 'rb')); rng = np.random.default_rng(int(f[-7:-4]) + 11)
    d, rp, amb = 0.5, 0.1, 5e-6
    A0 = isac.A_S
    if var == 'area':
        isac.A_S = val
    Y = isac.estimate(e, d, rp, rng, amb)
    nsd = np.tile(np.sqrt(isac.crb(e['Hdc'], d, rp, amb)), e['H'].shape[1])
    isac.A_S = A0
    win = 32.0
    if var == 'pairs':
        S = e['H'].shape[1]; keep = [s * S + k for s in range(val) for k in range(val)]
        Y, nsd = Y[:, keep], nsd[keep]
    if var == 'win':
        win = val
    k = val if var == 'gate' else 6.0
    valid = ~isac.artifact_mask(Y, k=k) & (e['lab'] == 0)
    W = isac.rr_eigen(Y, valid, nsd, win=win, seg=min(16.0, win / 2))
    rows = []
    Hm = e['H'].reshape(len(e['H']), -1).max(1)
    for t0, t1, rr, q in W:
        fr = slice(int(t0 * FS), int(t1 * FS))
        if e['apnea'][fr].any() or (e['lab'][fr] != 0).any():
            continue
        ref = isac.rr_reference(e['breaths'], t0, t1)
        if np.isfinite(ref):
            rows.append((ref, rr, float(Hm[fr].mean()), e['amp']))
    return (var, val), rows


def design():
    grid = [('area', v) for v in (1e-4, 4e-4, 9e-4)] + [('pairs', v) for v in (1, 2, 4, 6)] + [('win', v) for v in (16.0, 32.0, 64.0)] + [('gate', v) for v in (4.0, 6.0, 8.0)]
    out = {g: [] for g in grid}
    with Pool(2) as p:
        for g, rows in p.imap_unordered(design_job, [(f, *g) for g in grid for f in FILES[1::4]], chunksize=4):
            out[g] += rows
    res = {}
    for (var, val), rows in out.items():
        a = np.array(rows); ok = np.isfinite(a[:, 1]); err = np.abs(a[ok, 1] - a[ok, 0])
        r = dict(var=var, val=val, coverage=float(ok.mean()), within2=float((err < 2).mean()), mae=float(err.mean()))
        if var == 'win' and val == 32.0:
            r['by_ref'] = [float((err[(a[ok, 0] >= lo) & (a[ok, 0] < hi)] < 2).mean()) for lo, hi in ((0, 12), (12, 20), (20, 40))]
            r['n_ref'] = [int(((a[ok, 0] >= lo) & (a[ok, 0] < hi)).sum()) for lo, hi in ((0, 12), (12, 20), (20, 40))]
            g = a[ok, 2]; qs = np.quantile(g, [0, 0.25, 0.5, 0.75, 1.0])
            r['by_gain'] = [float((err[(g >= qs[i]) & (g <= qs[i + 1])] < 2).mean()) for i in range(4)]
            r['gain_q'] = [float(x) for x in qs]
            am = a[ok, 3]
            r['by_amp'] = [float((err[(am >= lo) & (am < hi)] < 2).mean()) for lo, hi in ((0.004, 0.0067), (0.0067, 0.0093), (0.0093, 0.0121))]
        res[f'{var}|{val}'] = r
        print(r, flush=True)
    json.dump(res, open('results/design3.json', 'w'), indent=1)


if __name__ == '__main__':
    import os
    os.makedirs('results', exist_ok=True)
    globals()[sys.argv[1]]()
