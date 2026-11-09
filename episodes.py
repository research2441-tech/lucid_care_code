"""LUCID-Care episodes: an in-bed patient under the ceiling of a ward, BIDMC respiration, inserted apnea,
bed-exit and fall events, staff occlusion; output per-frame body-reflection and static sensing gains."""
import os, sys, pickle, numpy as np
from multiprocessing import Pool
import raytrace as rt, rooms, ward, sense

_OU = ward._ou


def _quiet_ou(r, T, tau, sd):
    """Resting patients: trunk wobble slower and smaller than in the communication study (sd x0.02, tau >= 120 s)."""
    return _OU(r, T, max(tau, 120.0), sd * 0.02)


ward._ou = _quiet_ou

T_S = 480.0
FS = sense.FS


def static_gain(room):
    """Luminaire -> room surface -> photodiode first-order gain (S x K), same receiver model as the body path."""
    c, n, a, rho = rooms.patches(room)
    lum = room['lum']
    d1 = c[:, None] - lum[None]; r1 = np.linalg.norm(d1, axis=-1)
    E = 2 / (2 * np.pi * r1 ** 2) * np.clip(-d1[..., 2] / r1, 0, None) * np.clip(-(d1 * n[:, None]).sum(-1) / r1, 0, None)
    g = rt.conc_gain(1.5, sense.FOV_PD)
    d2 = lum[None] - c[:, None]; r2 = np.linalg.norm(d2, axis=-1)
    cth = np.clip((d2 * n[:, None]).sum(-1) / r2, 0, None); cps = np.clip(d2[..., 2] / r2, 0, None)
    R = rho[:, None] / np.pi * cth * cps * sense.A_PD * g / r2 ** 2 * (cps >= np.cos(sense.FOV_PD))
    return np.einsum('ps,pk,p->sk', E, R, a)


DWELL = float(os.environ.get('LC_DWELL', 300.0))
STAFF_MAX = int(os.environ.get('LC_STAFF', 4))


def kinematics(room, b, T, r, resp, ev):
    """Torso hip/neck/front/lateral over T frames with respiration displacement and events."""
    bd = room['beds'][b]
    nodes, occ, post = ward.in_bed(room, b, int(T * 100 / FS) + 10, r, dwell=DWELL)
    idx = (np.arange(T) * 100 / FS).astype(int)
    hip, neck = occ[0][0][idx].copy(), occ[0][1][idx].copy()
    front = nodes['chest'][1][idx].copy()
    s, xc = bd['sgn'], bd['xc']
    lab = np.zeros(T, np.int8)                        # 0 in bed, 1 transition, 2 out of bed (standing), 3 fall / on floor
    if ev['exit'] is not None:
        t0 = ev['exit']; dur = int(r.uniform(4, 8) * FS)
        side = r.choice([-1, 1])
        stand_hip = np.array([xc + side * 0.8, bd['y_head'] + s * 1.0, 0.95])
        for t in range(t0, T):
            a = min(1.0, (t - t0) / dur)
            h = hip[t0] * (1 - a) + stand_hip * a
            hip[t] = h; neck[t] = h + np.array([0, 0, 0.55]) * a + (neck[t0] - hip[t0]) * (1 - a)
            f = np.array([side, 0, 0.0]) * a + front[t0] * (1 - a); front[t] = f / np.linalg.norm(f)
            lab[t] = 1 if a < 1 else 2
        if ev['fall'] is not None:
            tf = ev['fall']; fd = int(r.uniform(0.6, 1.2) * FS); dirv = np.array([0, r.choice([-1, 1]) * 1.0, 0])
            h0, n0 = hip[tf].copy(), neck[tf].copy()
            for t in range(tf, T):
                a = min(1.0, (t - tf) / fd)
                hip[t] = h0 * (1 - a) + np.array([h0[0], h0[1], 0.12]) * a
                top = np.array([h0[0], h0[1], 0.15]) + dirv * 0.55
                neck[t] = n0 * (1 - a) + top * a
                f = np.array([0, 0, 1.0]) * a + np.array([side, 0, 0.0]) * (1 - a); front[t] = f / np.linalg.norm(f)
                lab[t] = 3
    axis = neck - hip
    lat = np.cross(axis, front); lat /= np.linalg.norm(lat, axis=1, keepdims=True)
    pst = post[idx]
    mov = np.zeros(T, bool); ch = np.nonzero(np.diff(pst) != 0)[0]
    for c in ch:
        mov[max(0, c - int(3 * FS)):c + int(3 * FS)] = True
    return hip, neck, front, lat, lab, mov


def episode(seed, rec, pool_layout=(3000, 16), p_exit=0.3, p_fall=0.5):
    r = np.random.default_rng(seed)
    room = ward.layout(pool_layout[0] + r.integers(pool_layout[1]))
    b = int(r.integers(len(room['beds'])))
    T = int(T_S * FS)
    x = rec['resp'][:int(T_S * rec['fs'])]
    t_src = np.arange(len(x)) / rec['fs']; t = np.arange(T) / FS
    x = np.interp(t, t_src, x)
    from scipy.signal import butter, sosfiltfilt
    x = sosfiltfilt(butter(2, [0.05, 1.0], btype='band', fs=FS, output='sos'), x)      # remove electrode baseline wander
    x = x / (2 * np.sqrt(2) * x.std() + 1e-12)                                          # unit peak-to-peak breath scale
    apnea = np.zeros(T, bool)
    for _ in range(r.integers(1, 3)):
        a0 = int(r.uniform(60, T_S - 60) * FS); al = int(r.uniform(10, 25) * FS)
        apnea[a0:a0 + al] = True
    x = np.where(apnea, 0.0, x)
    amp = r.uniform(0.004, 0.012)                      # peak-to-peak tidal excursion of the abdominal wall (m)
    ev = dict(exit=None, fall=None)
    if r.random() < p_exit:
        ev['exit'] = int(r.uniform(200, 380) * FS)
        if r.random() < p_fall:
            ev['fall'] = ev['exit'] + int(r.uniform(15, 60) * FS)
    hip, neck, front, lat, lab, mov = kinematics(room, b, T, r, x, ev)
    pos, area, w = sense.torso_patches(hip, neck, front, lat)
    disp = amp * x[:, None] * (0.6 * w[:, 0] + 1.0 * w[:, 1])[None]          # chest moves ~60% of the abdomen
    disp = np.where((lab == 0)[:, None], disp, 0.0)
    pos = pos + disp[..., None] * front[:, None]
    nrm = np.repeat(front[:, None], pos.shape[1], 1)
    near = np.argsort(np.linalg.norm(room['lum'][:, :2] - hip[0, :2], axis=1))[:6]
    lum = room['lum'][near]
    H = np.concatenate([sense.reflect_gain(lum, lum, pos[i:i + 2000], nrm[i:i + 2000], area) for i in range(0, T, 2000)])
    # staff occlusion of the body paths (one staff member visiting this bed now and then)
    blk = np.zeros((T, len(near)), bool)
    for _ in range(r.integers(0, STAFF_MAX)):
        s0 = int(r.uniform(0, T_S - 40) * FS); sl = int(r.uniform(8, 40) * FS)
        sp = hip[s0, :2] + np.array([r.choice([-1, 1]) * 0.8, 0])
        ca = np.zeros((1, 1, 3)); ca[..., :2] = sp; cb = ca.copy(); cb[..., 2] = 1.75
        bb = rt.blocked(hip[s0:s0 + 1] + front[s0:s0 + 1] * 0.16, lum, ca, cb, np.array([0.22]), [])
        blk[s0:s0 + sl] |= bb[0]
    H = H * (1 - 0.95 * (blk[:, :, None] | blk[:, None, :]))
    Hs = static_gain(room)[np.ix_(near, near)]
    tot_static = static_gain(room).sum(0)[near]          # DC illumination of every luminaire seen by each PD
    return dict(H=H.astype(np.float32), Hs=Hs, Hdc=tot_static, resp=x.astype(np.float32), apnea=apnea, lab=lab,
                breaths=rec['breaths'][rec['breaths'] < T_S], amp=amp, occl=blk.any(1), mov=mov, ev=ev, layout=int(room['L'] * 100), n_staff_blk=int(blk.any(1).sum()))


def job(a):
    i, rid, seed = a
    fn = f"{os.environ.get('LC_OUT', 'eps')}/ep_{i:03d}.pkl"
    if os.path.exists(fn):
        return fn
    recs = sense.bidmc()
    e = episode(seed, recs[rid]); e.update(rec=rid, seed=seed)
    pickle.dump(e, open(fn, 'wb'))
    return fn


if __name__ == '__main__':
    out = os.environ.get('LC_OUT', 'eps'); os.makedirs(out, exist_ok=True)
    n = int(os.environ.get('LC_N', 212)); s0 = int(os.environ.get('LC_SEED', 7000))
    jobs = [(i, i % 53, s0 + i) for i in range(n)]          # 4 episodes per BIDMC recording
    with Pool(2) as p:
        for f in p.imap_unordered(job, jobs):
            print(f, flush=True)
