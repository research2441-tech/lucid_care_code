"""Dual-function DCO-OFDM luminaires: pilot-based estimation of the luminaire->body->photodiode channel,
illumination/dimming constraints, rate-CRB trade-off, and the contactless respiration/apnea/event pipeline."""
import numpy as np, glob, pickle
from scipy.signal import butter, sosfiltfilt, welch
import phy, sense

Q = 1.602e-19
P_MAX, KAPPA, RESP, I_BG, B = 2.0, 3.2, 0.5, 50e-6, 50e6
TF = 1 / sense.FS
LK = phy.Link('dl')


def ac_std(d):
    """AC (OFDM) optical standard deviation allowed by the dimming level d (fraction of P_MAX) and LED range."""
    return np.minimum(d, 1 - d) * P_MAX / KAPPA


A_S = 4e-4                     # sensing photodiode area (m^2); traced gains used 1 cm^2 and scale linearly
K_B, TEMP = 1.381e-23, 295.0


def noise_density(Hdc, d, rho_p, amb=50e-6):
    """Noise current density (A^2/Hz) seen by a pilot band occupying the lowest rho_p*B of the spectrum:
    shot noise of the DC illumination and of the ambient light (amb per cm^2) plus Komine-Nakagawa thermal noise,
    whose FET term grows with frequency, averaged over the pilot band."""
    Bp = rho_p * B
    eta, I2, I3, G, Gam, gm = 112e-12 / 1e-4, 0.562, 0.0868, 10.0, 1.5, 30e-3
    ca = eta * A_S
    shot = 2 * Q * (RESP * d * P_MAX * Hdc * (A_S / sense.A_PD) + amb * (A_S / 1e-4))
    th = 8 * np.pi * K_B * TEMP * ca * I2 * Bp / G + 16 * np.pi ** 2 * K_B * TEMP * Gam * ca ** 2 * I3 * Bp ** 2 / gm
    return shot + th


def crb(Hdc, d, rho_p, amb=50e-6):
    """CRB of the per-frame estimate of the 1 cm^2-normalised channel gain (variance)."""
    v = noise_density(Hdc, d, rho_p, amb) / (RESP ** 2 * ac_std(d) ** 2 * rho_p * TF)
    return v / (A_S / sense.A_PD) ** 2


def rate(d, rho_p, Hc=2e-5):
    """Downlink data rate (bit/s) with the remaining (1 - rho_p) subcarriers, Paper-1 bit loading at SER 1e-3."""
    s = (LK.R * LK.K * ac_std(d)) ** 2 * Hc ** 2 * LK.led2
    c = (LK.R * ac_std(d)) ** 2 * LK.clip_var * Hc ** 2 * LK.led2
    n = LK.noise_var(np.array(d * P_MAX * Hc * 3))
    g = s / (n + c)
    bits = LK.load(g, phy.snr_threshold(1e-3))
    return (1 - rho_p) * LK.se(bits, 2) * B


def estimate(e, d, rho_p, rng, amb=50e-6):
    """Noisy per-frame channel estimates of all luminaire-photodiode pairs (T x S*K), static part removed."""
    H = e['H'] + e['Hs'][None]
    v = crb(e['Hdc'], d, rho_p, amb)                 # per photodiode
    noise = rng.standard_normal(H.shape) * np.sqrt(v)[None, None, :]
    return (H + noise).reshape(H.shape[0], -1)


SOS = butter(3, [0.09, 0.8], btype='band', fs=sense.FS, output='sos')


def resp_component(Y):
    """Band-pass every pair and project on the first principal component (sign-free)."""
    Z = sosfiltfilt(SOS, Y - Y.mean(0), axis=0)
    U, S, Vt = np.linalg.svd(Z[::5] / (Z[::5].std(0) + 1e-30), full_matrices=False)
    w = Vt[0] / (Z.std(0) + 1e-30)
    return Z @ w, S[0] ** 2 / (S ** 2).sum()


def rr_windows(sig, win=32.0, hop=8.0):
    fs = sense.FS; n = int(win * fs); h = int(hop * fs)
    out = []
    for s in range(0, len(sig) - n + 1, h):
        f, P = welch(sig[s:s + n], fs=fs, nperseg=int(16 * fs), noverlap=int(12 * fs), nfft=4096)
        m = (f > 0.08) & (f < 0.8)
        k = np.argmax(P[m])
        out.append((s / fs, s / fs + win, 60 * f[m][k], P[m][k] / P[m].sum()))
    return np.array(out)


def rr_reference(breaths, t0, t1):
    b = breaths[(breaths >= t0) & (breaths < t1)]
    return 60 * (len(b) - 1) / (b[-1] - b[0]) if len(b) >= 3 else np.nan


SOS_LP = butter(2, 0.3, btype='low', fs=sense.FS, output='sos')


def artifact_mask(Y, k=6.0):
    """Gross-motion / occlusion frames. Each pair is detrended with a 20 s running median, low-passed at 0.3 Hz
    and scaled by its median absolute deviation; frames whose pooled robust score exceeds k are flagged and the
    flag is dilated by 2 s."""
    from scipy.ndimage import median_filter, binary_dilation
    Ys = Y[::5]
    base = median_filter(Ys, size=(int(20 * sense.FS / 5), 1), mode='nearest')
    base = np.repeat(base, 5, 0)[:len(Y)]
    L = sosfiltfilt(SOS_LP, Y - base, axis=0)
    mad = np.median(np.abs(L - np.median(L, 0)), 0) + 1e-30
    sc = np.sqrt(np.mean((L / mad) ** 2, 1))
    return binary_dilation(sc > k, iterations=int(2 * sense.FS))


def rr_track(Y, valid, win=32.0, hop=8.0, max_art=0.05, noise_sd=None):
    """Gated respiration tracking: PCA over artifact-free frames, Welch peak per window; windows with more than
    max_art artifact frames are withheld (no estimate)."""
    fs = sense.FS
    Z = sosfiltfilt(SOS, Y - Y.mean(0), axis=0)
    if valid.sum() < 10 * fs:
        return np.array([(s / fs, s / fs + win, np.nan, 0.0) for s in range(0, len(Y) - int(win * fs) + 1, int(hop * fs))]), np.zeros(len(Y))
    Zc = Z[valid]
    sd = (noise_sd if noise_sd is not None else Zc.std(0)) + 1e-30          # noise-whitened PCA when the CRB is known
    U, S, Vt = np.linalg.svd(Zc[::5] / sd, full_matrices=False)
    sig = Z @ (Vt[0] / sd)
    n = int(win * fs); h = int(hop * fs); out = []
    for s in range(0, len(sig) - n + 1, h):
        art = 1 - valid[s:s + n].mean()
        if art > max_art:
            out.append((s / fs, s / fs + win, np.nan, 0.0)); continue
        f, P = welch(np.where(valid[s:s + n], sig[s:s + n], 0.0), fs=fs, nperseg=int(16 * fs), noverlap=int(12 * fs), nfft=4096)
        m = (f > 0.08) & (f < 0.8); k = np.argmax(P[m])
        out.append((s / fs, s / fs + win, 60 * f[m][k], P[m][k] / P[m].sum()))
    return np.array(out), sig


def rr_eigen(Y, valid, noise_sd, win=32.0, hop=8.0, max_art=0.05, seg=16.0):
    """Window-wise spectral eigen-combining: the cross-spectral matrix of the noise-whitened, band-passed pair
    signals is estimated from overlapping segments; its largest eigenvalue per frequency is the coherent
    respiration spectrum. Returns (t0, t1, rr, quality) per window; quality = peak / median eigenvalue."""
    fs = sense.FS
    Z = sosfiltfilt(SOS, Y - Y.mean(0), axis=0) / (noise_sd + 1e-30)
    Z = np.where(valid[:, None], Z, 0.0)
    n, h, ns = int(win * fs), int(hop * fs), int(seg * fs)
    nfft = 2048
    f = np.fft.rfftfreq(nfft, 1 / fs); m = (f > 0.08) & (f < 0.8)
    wnd = np.hanning(ns)[:, None]
    out = []
    for s in range(0, len(Z) - n + 1, h):
        if 1 - valid[s:s + n].mean() > max_art:
            out.append((s / fs, s / fs + win, np.nan, 0.0)); continue
        X = [np.fft.rfft(Z[a:a + ns] * wnd, nfft, axis=0)[m] for a in range(s, s + n - ns + 1, int(4 * fs))]
        X = np.stack(X)                                  # segs x F x P
        C = np.einsum('sfp,sfq->fpq', X, X.conj()) / len(X)
        ev = np.linalg.eigvalsh(C)[:, -1]
        k = np.argmax(ev)
        out.append((s / fs, s / fs + win, 60 * f[m][k], ev[k] / np.median(ev)))
    return np.array(out)
