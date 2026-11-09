"""Visible-light sensing of in-bed patients: luminaire light reflected by the torso surface into ceiling
photodiodes. Respiration displaces chest and abdomen patches along the torso normal (BIDMC impedance
waveform); bed-exit and falls move the whole torso. Gains per (luminaire, photodiode) pair."""
import numpy as np, scipy.io as sio
import raytrace as rt, rooms, ward

BIDMC = __import__('os').environ.get('BIDMC_MAT', './bidmc_data.mat')
FS = 50.0                      # sensing frame rate (Hz)
RHO_BODY = 0.45                # diffuse reflectance of gown / blanket
A_PD = 1e-4                    # sensing photodiode area (m^2)
FOV_PD = np.deg2rad(70)


def bidmc():
    m = sio.loadmat(BIDMC, squeeze_me=True, struct_as_record=False)['data']
    out = []
    for d in m:
        r = np.asarray(d.ref.resp_sig.imp.v, float); fs = float(d.ref.resp_sig.imp.fs)
        b = np.asarray(d.ref.breaths.ann1, float) / fs
        out.append(dict(resp=r, fs=fs, breaths=b, id=str(d.fix.id) if hasattr(d.fix, 'id') else ''))
    return out


def torso_patches(hip, neck, front, lat, nx=3, ny=4):
    """Front-surface patches of a torso (T x P x 3 positions, normals, areas, and a chest/abdomen weight)."""
    axis = neck - hip
    L = np.linalg.norm(axis, axis=-1, keepdims=True)
    u = axis / L
    pos, w = [], []
    for i in range(ny):
        a = (i + 0.5) / ny
        for j in range(nx):
            b = (j + 0.5) / nx - 0.5
            pos.append(hip + a * axis + b * 0.32 * lat + 0.16 * front)
            w.append((a, 1.0 - a))           # chest weight grows towards the neck
    pos = np.stack(pos, 1)
    area = np.full(pos.shape[1], (L.mean() / ny) * (0.32 / nx))
    return pos, area, np.array(w)


def reflect_gain(lum, pd, pos, nrm, area, m=1.0):
    """First-order luminaire -> body patch -> photodiode gain. lum: S x 3, pd: K x 3, pos: T x P x 3, nrm: T x P x 3."""
    g = rt.conc_gain(1.5, FOV_PD)
    d1 = pos[:, :, None, :] - lum[None, None]                 # T P S 3
    r1 = np.linalg.norm(d1, axis=-1)
    cphi = np.clip(-d1[..., 2] / r1, 0, None)
    cth1 = np.clip(-(d1 * nrm[:, :, None]).sum(-1) / r1, 0, None)
    E = (m + 1) / (2 * np.pi * r1 ** 2) * cphi ** m * cth1     # irradiance per W at each patch from each luminaire
    d2 = pd[None, None] - pos[:, :, None, :]                    # T P K 3
    r2 = np.linalg.norm(d2, axis=-1)
    cth2 = np.clip((d2 * nrm[:, :, None]).sum(-1) / r2, 0, None)
    cpsi = np.clip(d2[..., 2] / r2, 0, None)
    vis = cpsi >= np.cos(FOV_PD)
    R = RHO_BODY / np.pi * cth2 * cpsi * A_PD * g / r2 ** 2 * vis   # radiance-to-PD factor
    return np.einsum('tps,tpk,p->tsk', E, R, area)              # T x S x K
