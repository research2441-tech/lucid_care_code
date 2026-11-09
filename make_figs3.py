"""Figures for LUCID-Care (black and white)."""
import os, sys, json, pickle, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import isac

plt.rcParams.update({'font.family': 'serif', 'font.size': 7.5, 'axes.linewidth': 0.6, 'legend.frameon': False, 'mathtext.fontset': 'dejavuserif'})
LS = ['-', '--', ':', '-.', (0, (5, 1, 1, 1))]; GR = ['black', '0.35', '0.55', '0.7', '0.2']


def fig_arch():
    fig, ax = plt.subplots(figsize=(7.16, 3.3)); ax.set_xlim(0, 100); ax.set_ylim(0, 50); ax.axis('off')
    def box(x, y, w, h, t):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.2,rounding_size=0.8', fc='white', ec='black', lw=0.85))
        ax.text(x + w / 2, y + h / 2, t, ha='center', va='center', fontsize=6.4, linespacing=1.3)
        return (x, y, w, h)
    def frame(x, y, w, h, t):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.2,rounding_size=1.0', fc='#f3f3f3', ec='black', lw=0.6, ls='--'))
        ax.text(x + 1.2, y + h - 2.0, t, fontsize=6.8, style='italic', va='center')
    def ar(p, q, ls='-'):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=7, lw=0.8, color='black', ls=ls, shrinkA=0, shrinkB=0))
    frame(1, 26, 98, 23, 'Each luminaire: dual-function DCO-OFDM (illumination, data, low-band sensing pilots)')
    frame(1, 1, 98, 22, 'Ward controller: per-patient fusion and federated alarm calibration')
    a = box(3, 30, 16, 12, 'Dimming level $d$\nsets AC range\n$\\sigma(d)$, eq. (2)')
    b = box(23, 30, 16, 12, 'Pilot fraction $\\rho$\nin low band\n(rate vs. CRB)')
    c = box(43, 30, 16, 12, 'Ceiling\nphotodiodes:\nestimates $\\hat{h}_{\\ell k,t}$')
    dd = box(63, 30, 16, 12, 'Artifact gating\n+ noise-whitened\neigen-combining')
    e = box(83, 30, 14, 12, 'Contactless RR,\napnea, exit/fall\nscores')
    for p, q in ((a, b), (b, c), (c, dd), (dd, e)):
        ar((p[0] + p[2], p[1] + p[3] / 2), (q[0], q[1] + q[3] / 2))
    f = box(3, 5, 18, 12, 'Wearable chest patch\n(when its uplink\nis available)')
    g = box(26, 5, 18, 12, 'Quality-gated fusion\nof respiratory rate')
    h = box(49, 5, 22, 12, 'Mondrian conformal\nthresholds per ward\nand modality')
    i = box(76, 5, 21, 12, 'Federated\nquantile sketch\nacross wards')
    ar((f[0] + f[2], 11), (g[0], 11)); ar((g[0] + g[2], 11), (h[0], 11)); ar((i[0], 11), (h[0] + h[2], 11))
    ar((90, 30), (90, 17.2)); ax.text(91, 23.5, 'scores', fontsize=5.8)
    ar((70, 30), (35, 17.2)); ax.text(30, 23.9, 'contactless RR + quality', fontsize=5.8, va='center')
    fig.savefig('figs/fig1_arch.png', dpi=300, bbox_inches='tight')


def fig_signal():
    e = pickle.load(open('eps/ep_017.pkl', 'rb')); rng = np.random.default_rng(1)
    Y = isac.estimate(e, 0.5, 0.1, rng, 5e-6); nsd = np.tile(np.sqrt(isac.crb(e['Hdc'], 0.5, 0.1, 5e-6)), e['H'].shape[1])
    Z = isac.sosfiltfilt(isac.SOS, Y - Y.mean(0), axis=0) / nsd
    C = Z.T @ Z / len(Z); w = np.linalg.eigh(C)[1][:, -1]; s = Z @ w
    t = np.arange(len(s)) / 50
    R = isac.sosfiltfilt(isac.SOS, e['resp'])
    fig, ax = plt.subplots(3, 1, figsize=(7.16, 3.9), sharex=True)
    k = np.argmax(e['H'].reshape(len(Y), -1).mean(0))
    ax[0].plot(t, 1e9 * e['H'].reshape(len(Y), -1)[:, k], color='black', lw=0.6); ax[0].set_ylabel('body path\ngain ($\\times10^{-9}$)')
    ax[1].plot(t, R / R.std(), color='0.5', lw=0.6, label='BIDMC impedance (reference)')
    vm = ~isac.artifact_mask(Y) & (e['lab'] == 0) & ~e['apnea']
    sg = np.sign(np.corrcoef(s[vm], R[vm])[0, 1])
    ax[1].plot(t, np.where(vm | e['apnea'], sg * s / s[vm].std(), np.nan), color='black', lw=0.5, label='contactless estimate (artifact frames hidden)')
    ax[1].set_ylabel('normalised'); ax[1].legend(fontsize=6, loc='upper right', ncol=2); ax[1].set_ylim(-4, 6)
    for a in ax:
        for name, m, hatch in (('apnea', e['apnea'], 'black'), ('staff', e['occl'], '0.6')):
            idx = np.nonzero(np.diff(np.r_[0, m.astype(int), 0]))[0].reshape(-1, 2)
            for i0, i1 in idx:
                a.axvspan(i0 / 50, i1 / 50, color=hatch, alpha=0.15, lw=0)
    W = isac.rr_eigen(Y, ~isac.artifact_mask(Y) & (e['lab'] == 0), nsd)
    ax[2].plot((W[:, 0] + W[:, 1]) / 2, W[:, 2], 'o', ms=2.5, color='black', label='contactless RR')
    ref = [isac.rr_reference(e['breaths'], a0, a1) for a0, a1 in W[:, :2]]
    ax[2].plot((W[:, 0] + W[:, 1]) / 2, ref, '-', color='0.5', lw=0.8, label='annotated RR')
    ax[2].set_ylabel('RR (bpm)'); ax[2].set_xlabel('time (s)'); ax[2].legend(fontsize=6, loc='upper right', ncol=2); ax[2].set_xlim(0, 480)
    fig.tight_layout(h_pad=0.4); fig.savefig('figs/fig2_signal.png', dpi=300, bbox_inches='tight')


def fig_pareto():
    P = json.load(open('results/pareto.json'))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.16, 2.8))
    i = 0
    for amb, al in ((5e-6, 'night'), (5e-5, 'day')):
        for d in (0.2, 0.35, 0.5):
            pts = sorted([v for v in P.values() if v['amb'] == amb and v['d'] == d], key=lambda v: v['rp'])
            a.plot([v['rate'] / 1e6 for v in pts], [100 * v['within2'] for v in pts], color=GR[i % 5] if amb == 5e-6 else '0.6',
                   ls=LS[[0.2, 0.35, 0.5].index(d)], marker='o' if amb == 5e-6 else 's', ms=3, label=f'{al}, $d$={d}')
            b.loglog([v['rp'] for v in pts], [v['crb'] for v in pts], color='black' if amb == 5e-6 else '0.6', ls=LS[[0.2, 0.35, 0.5].index(d)], marker='o' if amb == 5e-6 else 's', ms=3)
            i += 1
    a.set_xlabel('downlink rate (Mbit/s)'); a.set_ylabel('RR windows within 2 bpm (%)'); a.legend(fontsize=5.6, ncol=2, loc='upper left'); a.set_title('(a)', fontsize=8)
    a.set_ylim(10, 85)
    b.set_xlabel('pilot fraction $\\rho$'); b.set_ylabel('CRB of $\\hat{h}$ per frame'); b.set_title('(b)', fontsize=8)
    fig.tight_layout(w_pad=1.5); fig.savefig('figs/fig3_pareto.png', dpi=300, bbox_inches='tight')


def fig_fusion():
    A = json.load(open('results/analysis3.json'))
    fig, (a, b, c3) = plt.subplots(1, 3, figsize=(7.16, 2.7))
    W = pickle.load(open('results/wards.pkl', 'rb'))
    for j, (key, vk, lab) in enumerate((('sw', 'vw', 'wearable'), ('sc', 'vc', 'contactless'))):
        ap = np.concatenate([r[key][(r['ap'] >= 0.8) & r[vk] & r['inb']] for r in W]); no = np.concatenate([r[key][(r['ap'] == 0) & r[vk] & r['inb']] for r in W])
        bins = np.linspace(-2, 6, 50)
        c3.hist(no, bins, density=True, histtype='step', color='0.6' if j else 'black', ls='-', label=f'{lab}, normal')
        c3.hist(ap, bins, density=True, histtype='step', color='0.6' if j else 'black', ls='--', label=f'{lab}, apnea')
    c3.set_xlabel('apnea score $s_{\\mathrm{ap}}$'); c3.set_ylabel('density'); c3.legend(fontsize=5.4, loc='upper right'); c3.set_title('(c)', fontsize=8); c3.set_ylim(0, 3.2)
    for i, k in enumerate(('contactless', 'wearable')):
        c = np.array(A['curve'][k]); a.plot(100 * c[:, 0], 100 * c[:, 1], color=GR[i], ls=LS[i], marker='os'[i], ms=3, label=k)
    fq = A['fusion_q']['fused']; a.plot([100 * fq['coverage']], [100 * fq['within2']], marker='*', ms=9, color='black', ls='none', label='fused (quality-gated)')
    a.set_xlabel('coverage (% of windows with an estimate)'); a.set_ylabel('within 2 bpm (%)'); a.legend(fontsize=6, loc='center', bbox_to_anchor=(0.42, 0.72)); a.set_title('(a)', fontsize=8)
    al = A['alarms']['fused|0.01']
    meths = [('fixed', 'fixed'), ('local', 'local'), ('federated', 'fed. normalised'), ('pooled', 'fed. pooled')]
    data = [100 * np.array(list(al[m]['far'].values())) for m, _ in meths]
    bp = b.boxplot(data, widths=0.5, patch_artist=True, medianprops=dict(color='black'))
    for p in bp['boxes']: p.set_facecolor('white')
    for j, dd in enumerate(data): b.plot(np.full(len(dd), j + 1) + np.random.default_rng(j).uniform(-0.12, 0.12, len(dd)), dd, 'o', ms=2, color='0.4')
    b.axhline(1.0, color='black', ls=':', lw=0.7); b.set_xticks(range(1, 5)); b.set_xticklabels(['fixed', 'local', 'fed.\nnorm.', 'fed.\npooled'], fontsize=6)
    b.set_ylabel('per-ward FAR (%)'); b.set_title('(b)', fontsize=8); b.tick_params(axis='x', labelsize=5.6)
    fig.tight_layout(w_pad=1.5); fig.savefig('figs/fig4_fusion_alarms.png', dpi=300, bbox_inches='tight')



def fig_design():
    Dz = json.load(open('results/design3.json'))
    fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.7))
    b = Dz['win|32.0']
    x = np.arange(4); ax[0].bar(x, [100 * v for v in b['by_gain']], 0.6, color='white', edgecolor='black', hatch='////', lw=0.6)
    ax[0].set_xticks(x); ax[0].set_xticklabels(['Q1', 'Q2', 'Q3', 'Q4'])
    ax[0].set_xlabel('body-path gain quartile'); ax[0].set_ylabel('within 2 bpm (%)'); ax[0].set_ylim(0, 100); ax[0].set_title('(a)', fontsize=8)
    x = np.arange(3); ax[1].bar(x, [100 * v for v in b['by_ref']], 0.6, color='0.3', edgecolor='black', lw=0.6)
    ax[1].set_xticks(x); ax[1].set_xticklabels(['<12', '12–20', '>20']); ax[1].set_xlabel('reference RR (bpm)'); ax[1].set_ylim(0, 100); ax[1].set_title('(b)', fontsize=8)
    for i, (var, vals, lab, sc) in enumerate((('area', (1e-4, 4e-4, 9e-4), 'PD area (cm$^2$)', 1e4), ('pairs', (1, 2, 4, 6), 'luminaires × PDs per side', 1), ('win', (16.0, 32.0, 64.0), 'window (s)', 1))):
        v = [Dz[f'{var}|{x}']['within2'] * 100 for x in vals]
        ax[2].plot(range(len(vals)), v, color=GR[i], ls=LS[i], marker='os^'[i], ms=3, label=lab + ': ' + ', '.join(f'{x * sc:g}' for x in vals))
    ax[2].set_xticks([0, 1, 2, 3]); ax[2].set_xticklabels(['1st', '2nd', '3rd', '4th']); ax[2].set_xlabel('setting (in the order listed)')
    ax[2].set_ylim(40, 70); ax[2].legend(fontsize=5.4, loc='lower left'); ax[2].set_title('(c)', fontsize=8)
    fig.tight_layout(w_pad=1.0); fig.savefig('figs/fig5_design.png', dpi=300, bbox_inches='tight')


def fig_calib():
    X = json.load(open('results/extra3.json'))
    fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.7))
    al = sorted(float(a) for a in X['alpha'])
    for i, (m, lab) in enumerate((('local', 'local'), ('federated', 'fed. normalised'), ('pooled', 'fed. pooled'), ('fixed', 'fixed'))):
        ax[0].loglog(al, [100 * X['alpha'][str(a)][m]['mx'] for a in al], color=GR[i], ls=LS[i], marker='os^d'[i], ms=3, label=lab)
    ax[0].loglog(al, [100 * a for a in al], color='0.6', lw=0.6, ls=':'); ax[0].set_xlabel('target $\\alpha$'); ax[0].set_ylabel('worst-ward FAR (%)')
    ax[0].set_xticks(al); ax[0].set_xticklabels([f'{a:g}' for a in al], fontsize=6); ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[0].set_yticks([1, 3, 10, 20]); ax[0].set_yticklabels(['1', '3', '10', '20']); ax[0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[0].legend(fontsize=5.6); ax[0].set_title('(a)', fontsize=8)
    fr = ['0.25', '0.5', '1.0']
    for i, m in enumerate(('local', 'federated', 'pooled')):
        ax[1].plot([25, 50, 100], [100 * X['budget'][f][m][0] for f in fr], color=GR[i], ls=LS[i], marker='os^'[i], ms=3, label=m)
    ax[1].set_xlabel('calibration data used (%)'); ax[1].set_ylabel('worst-ward FAR (%)'); ax[1].set_title('(b)', fontsize=8); ax[1].set_ylim(0, 18); ax[1].legend(fontsize=5.6, loc='upper right', ncol=3)
    r = np.array(X['roc']); ax[2].plot(r[:, 1], 100 * r[:, 2], 'o-', color='black', ms=3)
    for a, b_, c in r:
        ax[2].annotate(f'{a:g}', (b_, 100 * c), textcoords='offset points', xytext=(3, -8), fontsize=5.5)
    ax[2].set_xscale('log'); ax[2].set_xticks([0.3, 1, 3, 7]); ax[2].set_xticklabels(['0.3', '1', '3', '7']); ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[2].set_xlabel('false alarms per hour'); ax[2].set_ylabel('events detected (%)'); ax[2].set_title('(c)', fontsize=8); ax[2].set_ylim(40, 102)
    fig.tight_layout(w_pad=1.0); fig.savefig('figs/fig6_calib.png', dpi=300, bbox_inches='tight')



def fig_snr():
    o = json.load(open('results/snr3.json')); a = np.array(o['points']); a = a[np.isfinite(a[:, 1])]
    fig, ax = plt.subplots(figsize=(3.5, 2.7))
    ax.scatter(a[:, 0], 100 * a[:, 1], s=10, facecolors='white', edgecolors='black', lw=0.6)
    qs = np.quantile(a[:, 0], [0, 0.25, 0.5, 0.75, 1])
    cx = [(qs[i] + qs[i + 1]) / 2 for i in range(4)]; cy = [100 * a[(a[:, 0] >= qs[i]) & (a[:, 0] <= qs[i + 1]), 1].mean() for i in range(4)]
    ax.plot(cx, cy, 's-', color='0.4', ms=4, label='quartile means')
    ax.set_xlabel('combining SNR bound per 32 s window (dB)'); ax.set_ylabel('windows within 2 bpm (%)'); ax.set_ylim(-5, 105); ax.legend(fontsize=6, loc='lower right')
    fig.tight_layout(); fig.savefig('figs/fig7_snr.png', dpi=300, bbox_inches='tight')



def fig_analytic():
    e = pickle.load(open('eps/ep_000.pkl', 'rb')); hdc = e['Hdc'].mean()
    eta, I2, I3, G, Gam, gm = 112e-12 / 1e-4, 0.562, 0.0868, 10.0, 1.5, 30e-3
    ca = eta * isac.A_S; sc = isac.A_S / 1e-4
    c2 = 16 * np.pi ** 2 * isac.K_B * isac.TEMP * Gam * ca ** 2 * I3 * isac.B ** 2 / gm
    c1 = 8 * np.pi * isac.K_B * isac.TEMP * ca * I2 * isac.B / G
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.16, 2.8))
    amb = np.logspace(-6, -4, 40)
    for i, d in enumerate((0.2, 0.5, 0.8)):
        c0 = 2 * isac.Q * (isac.RESP * d * isac.P_MAX * hdc * sc + amb * sc)
        a.semilogx(amb * 1e6, np.sqrt(c0 / c2), color=GR[i], ls=LS[i], label=f'$d$={d}')
    a.set_xlabel('ambient photocurrent ($\\mu$A per cm$^2$)'); a.set_ylabel('optimal pilot fraction $\\rho^{\\star}$'); a.legend(fontsize=6); a.set_title('(a)', fontsize=8)
    dd = np.linspace(0.05, 0.95, 60)
    for i, (am, lab) in enumerate(((5e-6, 'night'), (5e-5, 'day'))):
        c0 = 2 * isac.Q * (isac.RESP * dd * isac.P_MAX * hdc * sc + am * sc)
        rs = np.sqrt(c0 / c2)
        crb = (2 * np.sqrt(c0 * c2) + c1) / (isac.RESP ** 2 * isac.ac_std(dd) ** 2 * isac.TF) / sc ** 2
        b.semilogy(dd, crb, color=GR[i * 2], ls=LS[i], label=f'min CRB, {lab}')
    b2 = b.twinx(); b2.plot(dd, [isac.rate(x, 0.1) / 1e6 for x in dd], color='0.6', ls=':', label='rate at $\\rho$=0.1')
    b.set_xlabel('dimming level $d$'); b.set_ylabel('minimum CRB'); b2.set_ylabel('downlink rate (Mbit/s)')
    h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels(); b.legend(h1 + h2, l1 + l2, fontsize=5.8, loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=3)
    b.set_title('(b)', fontsize=8)
    fig.tight_layout(w_pad=1.5); fig.savefig('figs/fig8_analytic.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    os.makedirs('figs', exist_ok=True)
    for fn in sys.argv[1:]:
        globals()[fn]()
