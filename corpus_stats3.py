"""Corpus statistics (Table III) and per-episode processing time of the contactless estimator."""
import glob, pickle, json, time, numpy as np
import isac

F = sorted(glob.glob('eps/ep_*.pkl'))
ap, ex, fa, oc, mv, amp, rr = 0, 0, 0, [], [], [], []
for f in F:
    e = pickle.load(open(f, 'rb'))
    a = e['apnea'] > 0
    ap += int((a & ~np.r_[False, a[:-1]]).sum())
    ex += e['ev']['exit'] is not None; fa += e['ev']['fall'] is not None
    oc.append(e['occl'].mean()); mv.append(e['mov'].mean()); amp.append(e['amp'])
    b = e['breaths']; rr.append(60 * (len(b) - 1) / (b[-1] - b[0]))
out = dict(episodes=len(F), apnea_segments=ap, exits=int(ex), falls=int(fa), occl=float(np.mean(oc)), mov=float(np.mean(mv)),
           amp=[float(min(amp)), float(max(amp))], rr=[float(min(rr)), float(np.median(rr)), float(max(rr))])
e = pickle.load(open(F[1], 'rb')); rng = np.random.default_rng(0)
Y = isac.estimate(e, 0.5, 0.1, rng, 5e-6)
nsd = np.tile(np.sqrt(isac.crb(e['Hdc'], 0.5, 0.1, 5e-6)), e['H'].shape[1])
t0 = time.perf_counter(); v = ~isac.artifact_mask(Y) & (e['lab'] == 0); t1 = time.perf_counter()
W = isac.rr_eigen(Y, v, nsd); t2 = time.perf_counter()
out.update(t_gate_s=t1 - t0, t_rr_s=t2 - t1, windows=len(W), t_episode_s=t2 - t0)
json.dump(out, open('results/stats3.json', 'w'), indent=1); print(out)
