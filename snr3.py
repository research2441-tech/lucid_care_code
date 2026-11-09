"""Link-budget check: predicted window SNR of the respiratory component vs. measured contactless accuracy per episode."""
import glob, pickle, json, numpy as np
from multiprocessing import Pool
import isac, sense

def job(f):
    e = pickle.load(open(f, 'rb')); rng = np.random.default_rng(int(f[-7:-4]) + 5)
    d, rp, amb = 0.5, 0.1, 5e-6
    crb = np.tile(isac.crb(e['Hdc'], d, rp, amb), e['H'].shape[1])
    H = e['H'].reshape(len(e['H']), -1)
    m = (e['lab'] == 0) & ~e['occl'] & ~e['mov'] & ~e['apnea']
    Z = isac.sosfiltfilt(isac.SOS, H - H.mean(0), axis=0)[m]
    frac = (0.8 - 0.09) / (sense.FS / 2)
    snr = (Z.var(0) / (crb * frac)).sum()                       # coherent-combining SNR per frame in band
    n = int(32 * sense.FS)
    snr_win = 10 * np.log10(snr * n * frac + 1e-12)              # after a 32 s window (band-limited samples)
    Y = isac.estimate(e, d, rp, rng, amb); nsd = np.sqrt(crb)
    W = isac.rr_eigen(Y, ~isac.artifact_mask(Y) & (e['lab'] == 0), nsd)
    ok = []
    for t0, t1, rr, q in W:
        fr = slice(int(t0 * 50), int(t1 * 50))
        if e['apnea'][fr].any() or (e['lab'][fr] != 0).any() or not np.isfinite(rr): continue
        ref = isac.rr_reference(e['breaths'], t0, t1)
        if np.isfinite(ref): ok.append(abs(rr - ref) < 2)
    return float(snr_win), float(np.mean(ok)) if ok else np.nan, len(ok)

if __name__ == '__main__':
    with Pool(2) as p:
        R = p.map(job, sorted(glob.glob('eps/ep_*.pkl'))[1::2])
    a = np.array(R)
    bins = [-30, 0, 5, 10, 15, 40]
    out = dict(points=a.tolist(), bins=bins, acc=[float(np.nanmean(a[(a[:, 0] >= lo) & (a[:, 0] < hi), 1])) for lo, hi in zip(bins[:-1], bins[1:])],
               n=[int(((a[:, 0] >= lo) & (a[:, 0] < hi)).sum()) for lo, hi in zip(bins[:-1], bins[1:])],
               spearman=float(__import__('scipy.stats', fromlist=['spearmanr']).spearmanr(a[:, 0], a[:, 1], nan_policy='omit')[0]))
    json.dump(out, open('results/snr3.json', 'w'), indent=1); print({k: v for k, v in out.items() if k != 'points'})
