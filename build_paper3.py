"""Assemble the LUCID-Care manuscript (IEEE template)."""
import json, pickle, numpy as np
from scipy.stats import wilcoxon
import docbuild as D
import cite as C
import paper3_text as P

TEMPLATE = './IEEE_template.docx'
OUT = './out/LUCID_Care_Q1_Paper.docx'
PA = json.load(open('results/pareto.json'))
AN = json.load(open('results/analysis3.json'))
AB = json.load(open('results/ablation3.json'))
ST = json.load(open('results/stats3.json'))
WR = pickle.load(open('results/wards.pkl', 'rb'))
F = './figs/'


def T(x): return C.resolve(x)
def pc(x): return f'{100 * x:.1f}'
def pa(d, rp, amb): return PA[f'{d}|{rp}|{amb}']


def add(d, items):
    for k, c in items:
        d.para(T(c)) if k == 'p' else d.eq(c)


# ------------------------------------------------------------------ statistics: fused vs contactless per episode (test half)
def per_episode():
    tc, tw = AN['quality_thresholds']
    rows = []
    for r in WR:
        if int(r['f'][-7:-4]) % 2 == 0 or not r['rows']:
            continue
        a = np.array(r['rows']); ref = a[:, 1]
        okc = np.isfinite(a[:, 2]) & (a[:, 3] >= tc); okw = np.isfinite(a[:, 4]) & (a[:, 5] >= tw)
        pick = okc & (~okw | (a[:, 3] / tc >= a[:, 5] / tw)); est = np.where(pick, a[:, 2], a[:, 4]); ok = okc | okw
        rows.append(dict(cov_c=okc.mean(), cov_f=ok.mean(), cov_w=okw.mean(),
                         acc_c=(np.abs(a[okc, 2] - ref[okc]) < 2).mean() if okc.any() else np.nan,
                         acc_f=(np.abs(est[ok] - ref[ok]) < 2).mean() if ok.any() else np.nan))
    return rows


EP = per_episode()


def wil(a, b):
    a, b = np.asarray(a), np.asarray(b); m = np.isfinite(a) & np.isfinite(b)
    res = wilcoxon(a[m], b[m], method='approx'); return float(res.pvalue), float(abs(res.zstatistic) / np.sqrt(m.sum())), int(m.sum())


S_COV = wil([e['cov_f'] for e in EP], [e['cov_w'] for e in EP])
S_COVC = wil([e['cov_f'] for e in EP], [e['cov_c'] for e in EP])
S_ACC = wil([e['acc_f'] for e in EP], [e['acc_c'] for e in EP])


def holm(ps):
    o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for i, j in enumerate(o):
        run = max(run, (m - i) * ps[j]); adj[j] = min(1, run)
    return adj


HP = holm([S_COV[0], S_COVC[0], S_ACC[0]])


def fp(p): return '$<10^{-4}$' if p < 1e-4 else f'{p:.4f}'


# ------------------------------------------------------------------ tables
PARAMS = [['Parameter', 'Value'],
          ['Luminaire $P_{\\max}$; clipping $\\kappa$; bandwidth $B$; $N_{\\mathrm{FFT}}$', '2 W; 3.2; 50 MHz; 128'],
          ['Dimming $d$; ambient $I_{\\mathrm{amb}}$ (per cm$^2$)', '0.2–0.5; 5 µA (night), 50 µA (day)'],
          ['Sensing photodiode: area, FOV, responsivity', '4 cm$^2$, 70°, 0.5 A/W'],
          ['Pilot fraction $\\rho$; frame $T_f$', '0.01–0.4; 20 ms (50 Hz)'],
          ['Torso: patches, reflectance $\\rho_b$; tidal excursion $a$', '12, 0.45; 4–12 mm'],
          ['Pairs per bed (luminaires × photodiodes)', '6 × 6'],
          ['Windows: RR / apnea / events', '32 s (hop 8 s) / 10 s (hop 2 s) / 10 s'],
          ['Wearable: absent episodes; outages', '15%; exp. 300 s on, 20 s off'],
          ['Alarm level $\\alpha$; target event false alarms', '1% and 5% per window; 1 per hour']]


APX = json.load(open('results/appendix3.json'))
APX_TVMED = {m: float(np.median([v['tv'] for v in APX[m].values()])) for m in ('contactless', 'wearable')}
APX_TVOVER = {m: sum(v['tv'] > v['dkw'] for v in APX[m].values()) for m in ('contactless', 'wearable')}
APX_POOL = {m: sum(v['delta'] < v['samp'] for v in APX['tail'][m].values()) for m in ('contactless', 'wearable')}


def APX_TABLE():
    rows = [['Ward', '$n_w$', 'Contactless $\\delta_w$ / $s_w$', 'Wearable $\\delta_w$ / $s_w$', 'Preferred (C / W)']]
    tc, tw = APX['tail']['contactless'], APX['tail']['wearable']
    for k in sorted(tc, key=lambda x: APX['contactless'][x]['n']):
        c, w = tc[k], tw[k]
        pref = ('pool' if c['delta'] < c['samp'] else 'local') + ' / ' + ('pool' if w['delta'] < w['samp'] else 'local')
        rows.append([k, str(APX['contactless'][k]['n']), '%.1f / %.1f' % (1e3 * c['delta'], 1e3 * c['samp']),
                     '%.1f / %.1f' % (1e3 * w['delta'], 1e3 * w['samp']), pref])
    return rows



EB = json.load(open('results/errband3.json'))


def BAND_TABLE():
    rows = [['Reference band (bpm)', 'Windows', 'Contactless', 'Wearable', 'Fused']]
    for b, v in EB['bands'].items():
        cells = [b.replace('-', '–'), str(v['n'])]
        for m in ('contactless', 'wearable', 'fused'):
            x = v[m]; cells.append('%.0f / %.0f / %+.1f' % (100 * x['coverage'], 100 * x['within2'], x['bias']))
        rows.append(cells)
    return rows


def RES_BAND():
    b, ba, t = EB['bands'], EB['ba'], EB['tachy']
    lo, no, hi = b['low (<12)'], b['normal (12-20)'], b['high (>20)']
    return [
        f"Table XIV splits the test windows by the reference rate. Slow breathing is the weak spot of every source. Below 12 bpm, the fused estimate is within 2 bpm in only {pc(lo['fused']['within2'])}% of accepted windows, against {pc(no['fused']['within2'])}% in the normal band and {pc(hi['fused']['within2'])}% above 20 bpm. Even the wearable drops to {pc(lo['wearable']['within2'])}%, with a positive bias of {lo['wearable']['bias']:.1f} bpm. A 32 s window holds fewer than six cycles of breathing below 12 bpm, so the Welch spectrum has few averages and its second harmonic, which falls inside the band, can win the peak search. The positive bias in this band, against a negative bias of the contactless estimate in the other two bands ({no['contactless']['bias']:+.1f} and {hi['contactless']['bias']:+.1f} bpm), points to two distinct failure modes: harmonic locking at low rates and locking onto residual posture motion below the true rate at higher rates.",
        f"Over all accepted windows, the Bland–Altman bias and 95% limits of agreement are {ba['fused']['bias']:+.2f} bpm and [{ba['fused']['lo']:.1f}, {ba['fused']['hi']:.1f}] bpm for the fused estimate, [{ba['wearable']['lo']:.1f}, {ba['wearable']['hi']:.1f}] bpm for the wearable and [{ba['contactless']['lo']:.1f}, {ba['contactless']['hi']:.1f}] bpm for the contactless estimate. The median absolute error is close to 0.5 bpm for all three, so the wide contactless limits come from a minority of gross errors (the 95th percentile of its absolute error is {ba['contactless']['p95']:.1f} bpm), which is the error pattern that the quality gate is meant to remove and does not remove completely. For the clinically relevant flag of a rate above 20 bpm, the fused estimate has a sensitivity of {pc(t['sens'])}% and a positive predictive value of {pc(t['ppv'])}% over {t['n_pos']} accepted tachypnoeic windows. A longer window or a harmonic check in the slow band are the obvious remedies, and the released code makes both easy to test."]



def CORPUS():
    return [['Item', 'Value'],
            ['Episodes (8 min, 50 Hz) / BIDMC records', f"{ST['episodes']} / 53 (four placements each)"],
            ['Wards (layouts) / calibration–test split', '16 / alternate episodes (106 + 106)'],
            ['Reference respiratory rate (min / median / max)', f"{ST['rr'][0]:.1f} / {ST['rr'][1]:.1f} / {ST['rr'][2]:.1f} bpm"],
            ['Inserted apneas (10–25 s)', f"{ST['apnea_segments']}"],
            ['Bed exits / falls', f"{ST['exits']} / {ST['falls']}"],
            ['Frames with staff occlusion / posture change', f"{pc(ST['occl'])}% / {pc(ST['mov'])}%"],
            ['Episodes without a wearable', f"{AN['n']['wear_absent']}"]]


def PARETO_TABLE():
    r = [['Ambient', '$d$', '$\\rho$', 'Rate (Mb/s)', 'CRB ($10^{-19}$)', 'Within 2 bpm (%)', 'MAE (bpm)', 'Coverage (%)']]
    for amb, lab in ((5e-6, 'Night'), (5e-5, 'Day')):
        for d in (0.2, 0.35, 0.5):
            for rp in (0.01, 0.03, 0.1, 0.2, 0.4):
                v = pa(d, rp, amb)
                r.append([lab, f'{d}', f'{rp}', f"{v['rate'] / 1e6:.0f}", f"{v['crb'] / 1e-19:.2f}", pc(v['within2']), f"{v['mae']:.2f}", pc(v['coverage'])])
    return r


def FUSION_TABLE():
    r = [['Condition ($d$ / ambient)', 'Contactless', 'Wearable', 'Fused']]
    def cell(x): return f"{pc(x['coverage'])} / {pc(x['within2'])} / {x['mae']:.2f}"
    r.append(['All wards, no quality gate'] + [cell(AN['fusion_all'][k]) for k in ('contactless', 'wearable', 'fused')])
    r.append(['All wards, quality-gated'] + [cell(AN['fusion_q'][k]) for k in ('contactless', 'wearable', 'fused')])
    for k, v in AN['fusion_by_cond'].items():
        d, amb = k.split('|')
        r.append([f"{d} / {'night' if float(amb) < 1e-5 else 'day'}"] + [cell(v[m]) for m in ('contactless', 'wearable', 'fused')])
    return r


def ABL_TABLE():
    lab = dict(full='Full estimator (gating + whitening + eigen-combining)', nogate='Without artifact gating', pca='PCA over the whole episode instead of per-window eigen-combining',
               nowhite='Without noise whitening', best='Single best pair selected with ground truth (oracle)')
    r = [['Variant', 'Coverage (%)', 'Within 2 bpm (%)', 'MAE (bpm)']]
    for k in ('full', 'nogate', 'pca', 'nowhite', 'best'):
        v = AB[k]; r.append([lab[k], pc(v['coverage']), pc(v['within2']), f"{v['mae']:.2f}"])
    return r


def ALARM_TABLE():
    r = [['Modality / $\\alpha$', 'Method', 'Mean FAR (%)', 'Max ward FAR (%)', 'Wards > 2$\\alpha$ (%)', 'Apnea sensitivity (%)']]
    names = dict(fixed='Fixed threshold', local='Local per ward', federated='Federated, normalised', pooled='Federated, pooled')
    for mode in ('wearable', 'contactless', 'fused'):
        for a in ('0.01',):
            v = AN['alarms'][f'{mode}|{a}']
            for m in ('fixed', 'local', 'federated', 'pooled'):
                x = v[m]; r.append([f'{mode} / {a}', names[m], f"{100 * x['far_mean']:.2f}", f"{100 * x['far_max']:.2f}", f"{100 * x['frac_wards_over']:.0f}", pc(x['sens'])])
    v = AN['alarms']['fused|0.05']
    for m in ('local', 'pooled'):
        x = v[m]; r.append(['fused / 0.05', names[m], f"{100 * x['far_mean']:.2f}", f"{100 * x['far_max']:.2f}", f"{100 * x['frac_wards_over']:.0f}", pc(x['sens'])])
    return r


DZ = json.load(open('results/design3.json'))
XA = json.load(open('results/extra3.json'))


def DESIGN_TABLE():
    lab = {'area': 'Photodiode area {} cm$^2$', 'pairs': 'Nearest {0} luminaires × {0} photodiodes', 'win': 'Window length {} s', 'gate': 'Gate threshold {}'}
    r = [['Setting', 'Coverage (%)', 'Within 2 bpm (%)', 'MAE (bpm)']]
    for k, v in DZ.items():
        var, val = k.split('|'); val = float(val)
        txt = lab[var].format(f'{val * 1e4:g}' if var == 'area' else f'{val:g}')
        r.append([txt, pc(v['coverage']), pc(v['within2']), f"{v['mae']:.2f}"])
    return r


def RES_DESIGN():
    b = DZ['win|32.0']
    return [f"Table VIII and Fig. 5 vary the estimator and the receiver one factor at a time. Accuracy depends mostly on how much light the torso returns. In the lowest quartile of the strongest body-path gain only {pc(b['by_gain'][0])}% of windows fall within 2 bpm, against {pc(b['by_gain'][2])}% in the third quartile. Low gains correspond to lateral postures and to beds whose luminaire stands to one side, so posture explains much of the spread in Fig. 2. Slow breathing is harder: only {pc(b['by_ref'][0])}% of windows below 12 bpm are within 2 bpm ({b['n_ref'][0]} windows), against {pc(b['by_ref'][1])}% between 12 and 20 bpm, because fewer cycles fit in a 32 s window and the lowest rates approach the band edge. A 64 s window raises accuracy to {pc(DZ['win|64.0']['within2'])}% but lowers coverage to {pc(DZ['win|64.0']['coverage'])}%, because long windows are more often hit by an artifact. Enlarging the photodiode from 1 to 9 cm$^2$ raises accuracy from {pc(DZ['area|0.0001']['within2'])}% to {pc(DZ['area|0.0009']['within2'])}%. Using fewer pairs helps: the nearest luminaire and photodiode alone give {pc(DZ['pairs|1']['within2'])}% and the full 6 × 6 set gives {pc(DZ['pairs|6']['within2'])}%, which confirms the ablation finding that distant pairs add more noise than signal to the combiner. The gate threshold trades coverage against accuracy smoothly."]


def BUDGET_TABLE():
    r = [['Calibration data used', 'Local', 'Federated, normalised', 'Federated, pooled']]
    for f in ('0.25', '0.5', '1.0'):
        r.append([f'{float(f) * 100:.0f}%'] + [f"{100 * XA['budget'][f][m][0]:.1f} / {100 * XA['budget'][f][m][1]:.0f}" for m in ('local', 'federated', 'pooled')])
    return r


def RES_BUDGET():
    a = XA['alpha']; b = XA['budget']
    return [f"Fig. 6(a) repeats the comparison for targets from 0.5% to 10%. Federated pooling keeps the worst ward closest to the target at every level, for example {100 * a['0.005']['pooled']['mx']:.1f}% for a 0.5% target, where local calibration gives {100 * a['0.005']['local']['mx']:.1f}%. The fixed threshold ignores the target. Table X and Fig. 6(b) shrink the calibration budget. With a quarter of the calibration episodes, the worst local ward reaches {100 * b['0.25']['local'][0]:.1f}% and {100 * b['0.25']['local'][1]:.0f}% of wards exceed $2\\alpha$, while the pooled threshold stays at {100 * b['0.25']['pooled'][0]:.1f}%. This is the regime described by Proposition 2: each ward alone has too few windows, and the differences between wards are smaller than the sampling error of a local quantile. Fig. 6(c) shows the operating curve of the event detector. Allowing two false alarms per hour raises detection to {100 * XA['roc'][3][2]:.0f}%, and a quarter of an alarm per hour still finds {100 * XA['roc'][0][2]:.0f}% of events."]


SN = json.load(open('results/snr3.json'))


def WARD_TABLE():
    tc, tw = AN['quality_thresholds']
    al = AN['alarms']['fused|0.01']
    r = [['Ward', 'Condition', 'Episodes', 'Coverage (%)', 'Within 2 bpm (%)', 'FAR local (%)', 'FAR pooled (%)']]
    wards = sorted({x['ward'] for x in WR})
    for i, w in enumerate(wards):
        rs = [x for x in WR if x['ward'] == w and int(x['f'][-7:-4]) % 2 == 1 and len(x['rows'])]
        if not rs: continue
        a = np.concatenate([np.array(x['rows']) for x in rs]); ref = a[:, 1]
        okc = np.isfinite(a[:, 2]) & (a[:, 3] >= tc); okw = np.isfinite(a[:, 4]) & (a[:, 5] >= tw)
        pick = okc & (~okw | (a[:, 3] / tc >= a[:, 5] / tw)); est = np.where(pick, a[:, 2], a[:, 4]); ok = okc | okw
        cond = f"{rs[0]['d']} / {'night' if rs[0]['amb'] < 1e-5 else 'day'}"
        r.append([str(i + 1), cond, str(len(rs)), pc(ok.mean()), pc((np.abs(est[ok] - ref[ok]) < 2).mean()) if ok.any() else '—',
                  f"{100 * al['local']['far'][str(w)] if str(w) in al['local']['far'] else 100 * al['local']['far'][w]:.2f}",
                  f"{100 * al['pooled']['far'][str(w)] if str(w) in al['pooled']['far'] else 100 * al['pooled']['far'][w]:.2f}"])
    return r


def RES_WARD():
    pts = np.array(SN['points']); pts = pts[np.isfinite(pts[:, 1])]
    qs = np.quantile(pts[:, 0], [0, 0.25, 0.5, 0.75, 1])
    acc = [pts[(pts[:, 0] >= qs[i]) & (pts[:, 0] <= qs[i + 1]), 1].mean() for i in range(4)]
    return [f"Table XI lists the outcome ward by ward. Coverage and accuracy follow the lighting condition of each ward, and the pooled thresholds keep every ward's false-alarm rate within a few percent, while local thresholds leave several wards far above the target. Fig. 7 compares contactless accuracy with a link-budget bound. For each test episode the bound is the in-band variance of the noise-free body signal summed over pairs, divided by the in-band CRB, as if every pair were combined perfectly. The bound exceeds 6 dB per window in every episode and has a median of {np.median(pts[:, 0]):.0f} dB. Yet accuracy changes little across its quartiles ({pc(acc[0])}%, {pc(acc[1])}%, {pc(acc[2])}%, {pc(acc[3])}%; Spearman $\\rho_s={SN['spearman']:.2f}$). Most errors therefore do not come from noise in the sense of (6). They come from the estimator, which cannot combine pairs perfectly when many of them carry little respiratory signal, and from residual slow motion inside the respiratory band. This agrees with the finding of Table VIII that fewer, stronger pairs work better, and it shows where the remaining accuracy lies: in better pair selection and motion modelling rather than in more pilot power."]


EX = json.load(open('results/external3.json'))


def EXT_TABLE():
    r = [['Metric', 'Main test half', 'Stress-test wards', 'Record-disjoint test']]
    f = lambda x: f"{pc(x['coverage'])} / {pc(x['within2'])}"
    r.append(['Contactless RR: coverage / within 2 bpm (%)', f(EX['main_fusion']['contactless']), f(EX['ood_fusion']['contactless']), f(EX['disjoint_fusion']['contactless'])])
    r.append(['Fused RR: coverage / within 2 bpm (%)', f(EX['main_fusion']['fused']), f(EX['ood_fusion']['fused']), f(EX['disjoint_fusion']['fused'])])
    m = AN['alarms']['fused|0.01']
    for meth, lab in (('pooled', 'Apnea FAR, pooled: mean / max (%)'), ('local', 'Apnea FAR, local: mean / max (%)')):
        r.append([lab, f"{100 * m[meth]['far_mean']:.2f} / {100 * m[meth]['far_max']:.1f}", f"{100 * EX['ood_alarms'][meth]['mean']:.2f} / {100 * EX['ood_alarms'][meth]['mx']:.1f}",
                  f"{100 * EX['disjoint_alarms'][meth]['mean']:.2f} / {100 * EX['disjoint_alarms'][meth]['mx']:.1f}"])
    ev = AN['events']
    r.append(['Events: detected (%) / false alarms per hour', f"{pc(ev['sens'])} / {ev['fa_per_hour']:.2f}", f"{pc(EX['ood_events']['sens'])} / {EX['ood_events']['fa_per_hour']:.2f}",
              f"{pc(EX['disjoint_events']['sens'])} / {EX['disjoint_events']['fa_per_hour']:.2f}"])
    return r


def RES_EXT():
    o, dj = EX['ood_fusion'], EX['disjoint_fusion']
    return [f"Two checks probe generalisation (Table XII). The stress-test wards contain 106 new episodes with posture changes every 90 s on average instead of 300 s and up to eight staff visits instead of three. They are scored with the quality gates, alarm thresholds and event detector fitted on the main corpus. Contactless accuracy drops from {pc(EX['main_fusion']['contactless']['within2'])}% to {pc(o['contactless']['within2'])}% within 2 bpm, because more windows contain residual motion, but fused accuracy holds at {pc(o['fused']['within2'])}% with {pc(o['fused']['coverage'])}% coverage. The gate removes most of the extra motion, and the wearable carries the rest. Alarm thresholds transfer less well. The pooled false-alarm rate rises to {100 * EX['ood_alarms']['pooled']['mean']:.1f}% on average, worst {100 * EX['ood_alarms']['pooled']['mx']:.1f}%, against {100 * EX['ood_alarms']['local']['mx']:.1f}% for local thresholds. Event detection degrades most, to {pc(EX['ood_events']['sens'])}% at {EX['ood_events']['fa_per_hour']:.1f} false alarms per hour, with {pc(EX['ood_events']['fa_share_posture'])}% of false alarms caused by the more frequent posture changes. Thresholds and detectors must therefore be recalibrated when a ward's activity pattern changes. The online calibration of the companion link controller would be a natural mechanism for this. In the record-disjoint split, no patient's respiration appears in both calibration and test, and the results stay close to the main split: fused coverage {pc(dj['fused']['coverage'])}% with {pc(dj['fused']['within2'])}% within 2 bpm, and a pooled worst-ward false-alarm rate of {100 * EX['disjoint_alarms']['pooled']['mx']:.1f}%. Sharing records between the halves therefore did not inflate the main results."]


WE = json.load(open('results/wear3.json'))
ET = json.load(open('results/events_type3.json'))


def WEAR_TABLE():
    r = [['Wearable availability', 'Wearable', 'Contactless', 'Fused', 'Apnea sensitivity (%)']]
    for k, v in WE.items():
        f = v['fusion']; c = lambda x: f"{pc(x['coverage'])} / {pc(x['within2'])}"
        r.append([k, c(f['wearable']), c(f['contactless']), c(f['fused']), pc(v['sens'])])
    return r


def RES_WEAR():
    a, b = WE['no absent patches'], WE['50% absent']
    fq = WE['15% absent, outages every 60 s']
    return [f"Table XIII varies the availability of the wearable while the contactless channel is unchanged. The gain from fusion grows as the patch becomes less reliable. With every patch present, fusion adds {pc(a['fusion']['fused']['coverage'] - a['fusion']['wearable']['coverage'])} points of coverage. With half the patches absent it adds {pc(b['fusion']['fused']['coverage'] - b['fusion']['wearable']['coverage'])} points, and with outages every minute it adds {pc(fq['fusion']['fused']['coverage'] - fq['fusion']['wearable']['coverage'])} points. Accuracy of the fused stream degrades gracefully, from {pc(a['fusion']['fused']['within2'])}% to {pc(b['fusion']['fused']['within2'])}% within 2 bpm. Apnea sensitivity, by contrast, falls in proportion to wearable availability, from {pc(a['sens'])}% to {pc(b['sens'])}%. The contactless channel cannot replace the patch for apnea."]


def STAT_TABLE():
    return [['Comparison (per test episode)', 'n', '$p_{\\mathrm{Holm}}$', '$r$'],
            ['Coverage: fused vs. wearable only', str(S_COV[2]), fp(HP[0]), f'{S_COV[1]:.2f}'],
            ['Coverage: fused vs. contactless only', str(S_COVC[2]), fp(HP[1]), f'{S_COVC[1]:.2f}'],
            ['Within 2 bpm: fused vs. contactless only', str(S_ACC[2]), fp(HP[2]), f'{S_ACC[1]:.2f}']]


# ------------------------------------------------------------------ text
def ABSTRACT():
    b = pa(0.5, 0.1, 5e-6); w = pa(0.2, 0.1, 5e-5); q = AN['fusion_q']; al = AN['alarms']['fused|0.01']; ev = AN['events']
    return (f"Respiratory rate is the earliest sign of deterioration on hospital wards, yet wearable patches go missing and manual counts are rare. "
            f"This paper turns ceiling luminaires into dual-function visible-light sensing and communication nodes. Body-reflected DCO-OFDM pilots "
            f"are estimated by ceiling photodiodes and fused with an intermittent wearable patch. Thresholds are calibrated per ward so that every "
            f"ward keeps its promised false-alarm rate. Placing the sensing pilots in the lowest part of the band, where receiver thermal noise is smallest, gives the channel-estimation "
            f"Cramér–Rao bound an interior minimum in the pilot fraction, and this sets a rate–accuracy frontier constrained by dimming and ambient light. "
            f"Artifact gating followed by noise-whitened spectral eigen-combining estimates respiratory rate from 36 luminaire–photodiode pairs. "
            f"Evaluation uses 212 ray-traced eight-minute episodes in which the impedance respiration of the 53 BIDMC recordings moves the torso of a patient lying under the ceiling, with inserted "
            f"apneas, bed exits, falls, posture changes and staff occlusions. At night with half dimming, 10% pilots keep a {b['rate'] / 1e6:.0f} Mb/s downlink and place "
            f"{pc(b['within2'])}% of contactless estimates within 2 bpm; in daylight at 20% dimming this falls to {pc(w['within2'])}%. Quality-gated fusion "
            f"reports a rate for {pc(q['fused']['coverage'])}% of windows, against {pc(q['wearable']['coverage'])}% for the patch alone, with {pc(q['fused']['within2'])}% within 2 bpm. "
            f"Contactless sensing could not detect apneas at these signal levels, so alarms rely on the patch when it is present. Federated pooled calibration "
            f"kept the worst ward's false-alarm rate at {100 * al['pooled']['far_max']:.1f}% for a 1% target, against {100 * al['local']['far_max']:.1f}% for local calibration. "
            f"A boosted detector found {pc(ev['sens'])}% of bed exits and falls at {ev['fa_per_hour']:.2f} false alarms per hour.")


def RES():
    b, w = pa(0.5, 0.1, 5e-6), pa(0.2, 0.1, 5e-5)
    lo, hi = pa(0.5, 0.01, 5e-6), pa(0.5, 0.4, 5e-6)
    out = {}
    out['A'] = [
f"Fig. 3 and Table V give the frontier. For every condition, the share of windows within 2 bpm first rises with the pilot fraction and then falls. At night with $d=0.5$, it rises from {pc(lo['within2'])}% at $\\rho=0.01$ to {pc(b['within2'])}% at $\\rho=0.1$ and returns to {pc(hi['within2'])}% at $\\rho=0.4$. The rate falls from {lo['rate'] / 1e6:.0f} to {b['rate'] / 1e6:.0f} and {hi['rate'] / 1e6:.0f} Mb/s over the same range. The CRB in Fig. 3(b) has its minimum at $\\rho\\approx0.1$ for every dimming level, as predicted by (5)–(6). Beyond that point the pilot band reaches frequencies where the transistor noise of the receiver grows with $\\rho^2$. Dimming acts through the AC swing (1): $d=0.2$ allows 40% of the swing of $d=0.5$, raises the CRB by a factor of about six and costs about 20 points of accuracy. Daylight adds shot noise and costs a further 12 points or so. Moving from $d=0.5$ to $d=0.2$ in daylight reduces the share within 2 bpm from {pc(pa(0.5, 0.1, 5e-5)['within2'])}% to {pc(w['within2'])}%. Night-time monitoring with dimmed but not dark lighting is therefore the favourable regime, and sensing cannot be added to deep-dimmed luminaires without shortening the dark period.",
f"Fig. 2 shows one episode. For the first three minutes the patient lies in a posture whose torso reflects little light towards the photodiodes ($4\\times10^{{-9}}$), and the contactless estimates scatter. After a posture change the body-path gain rises tenfold, and the estimates follow the annotated rate within 1 bpm until the end of the episode. Staff occlusion and the posture change are withheld by the gate. The contactless channel is therefore reliable for some postures and not others, which is why it is combined with the wearable rather than used alone.",
    ]
    q, a0 = AN['fusion_q'], AN['fusion_all']
    out['B'] = [
f"Table VI lists coverage, accuracy and mean absolute error for each source. Without a quality gate, the contactless estimator reports for {pc(a0['contactless']['coverage'])}% of windows but only {pc(a0['contactless']['within2'])}% of them fall within 2 bpm. The wearable reports for {pc(a0['wearable']['coverage'])}% of windows, because 15% of episodes have no patch and the rest have uplink outages, and {pc(a0['wearable']['within2'])}% of its estimates are within 2 bpm. The quality gates chosen on the calibration half keep {pc(q['contactless']['coverage'])}% of contactless windows at {pc(q['contactless']['within2'])}% accuracy. Fusion reports for {pc(q['fused']['coverage'])}% of windows, {pc(q['fused']['coverage'] - q['wearable']['coverage'])} points more than the patch alone, with {pc(q['fused']['within2'])}% within 2 bpm and a mean absolute error of {q['fused']['mae']:.2f} bpm. Fig. 4(a) shows that contactless quality is not monotone: the windows with the very highest peak ratio are often those where residual motion produced a strong low-frequency peak. The gate is therefore set at a moderate quality level instead of the highest.",
f"The benefit depends on the lighting condition of the ward. Under $d=0.5$ at night, fusion covers {pc(AN['fusion_by_cond']['0.5|5e-06']['fused']['coverage'])}% of windows with {pc(AN['fusion_by_cond']['0.5|5e-06']['fused']['within2'])}% within 2 bpm. Under $d=0.2$ at night, the contactless share is lower, and fusion covers {pc(AN['fusion_by_cond']['0.2|5e-06']['fused']['coverage'])}% at {pc(AN['fusion_by_cond']['0.2|5e-06']['fused']['within2'])}%.",
    ]
    out['C'] = [
f"Table VII separates the parts of the estimator at night with $d=0.5$ and $\\rho=0.1$. Artifact gating raises the share within 2 bpm from {pc(AB['nogate']['within2'])}% to {pc(AB['full']['within2'])}%, at the cost of withholding {pc(1 - AB['full']['coverage'])}% of windows. Noise whitening matters most: without it the combiner weights pairs by their total variance, which is dominated by noise on weak pairs, and accuracy collapses to {pc(AB['nowhite']['within2'])}%. Per-window eigen-combining is slightly better than a single principal component over the whole episode ({pc(AB['full']['within2'])}% against {pc(AB['pca']['within2'])}%), because it adapts to posture changes. Selecting the single strongest pair with knowledge of the true body gain reaches {pc(AB['best']['within2'])}%, so the combiner still loses some accuracy by weighting weak pairs. An adaptive pair selection that does not need ground truth is a natural next step.",
    ]
    al1 = AN['alarms']['fused|0.01']; alw = AN['alarms']['wearable|0.01']; alc = AN['alarms']['contactless|0.01']
    out['D'] = [
f"Table IX reports alarm calibration on the test half. The wearable apnea score separates apneas cleanly: at $\\alpha=1\\%$ every calibration method detects about {100 * alw['pooled']['sens']:.0f}% of apnea windows. The contactless score does not. Its apnea sensitivity stays at {pc(alc['pooled']['sens'])}% at 1% false alarms, because a 10 s window holds too little respiratory energy relative to estimation noise to tell an apnea from a quiet breath. Fused alarms therefore detect {pc(al1['pooled']['sens'])}% of apnea windows, the share in which the patch is present. This is a negative result: contactless monitoring adds rate coverage but not apnea detection at ceiling distance. Fig. 4(c) shows why. The wearable score of apnea windows lies far to the right of its normal windows, because an apnea removes nearly all in-band power from a clean chest signal. The contactless score of apnea windows overlaps its normal windows almost entirely. Its respiratory-band power in a 10 s window is dominated by estimation noise and by residual slow motion, so stopping breathing changes it by less than the window-to-window spread. Longer apnea windows would separate the classes better but would delay alarms beyond clinical usefulness.",
f"Calibration decides whether each ward keeps its promise. A fixed threshold gives a mean false-alarm rate of {100 * al1['fixed']['far_mean']:.1f}% and exceeds $2\\alpha$ in {100 * al1['fixed']['frac_wards_over']:.0f}% of wards. Local per-ward conformal thresholds meet the target on average ({100 * al1['local']['far_mean']:.2f}%) but with only a few hundred calibration windows per ward the worst ward reaches {100 * al1['local']['far_max']:.1f}% and {100 * al1['local']['frac_wards_over']:.0f}% of wards exceed $2\\alpha$. The federated pooled threshold lowers the worst ward to {100 * al1['pooled']['far_max']:.1f}%, and only {100 * al1['pooled']['frac_wards_over']:.0f}% of wards exceed $2\\alpha$. The normalised federated variant lies in between ({100 * al1['federated']['far_max']:.1f}%). Its per-ward medians and scales are themselves estimated from small samples, which adds variance instead of removing heterogeneity. Fig. 4(b) shows the per-ward spread. Sharing quantile sketches is therefore worthwhile even when wards differ, as long as each ward keeps its own modality partition.",
    ]
    ev = AN['events']
    out['E'] = [
f"The bed-exit and fall detector, trained on one quarter of the episodes and calibrated on another quarter, found {ev['detected']} of {ev['events']} events in the test half ({pc(ev['sens'])}%) at {ev['fa_per_hour']:.2f} false alarms per hour, close to the target of one per hour. Of the false alarms, {pc(ev['fa_share_posture'])}% coincided with posture changes in bed and {pc(ev['fa_share_occlusion'])}% with staff occlusion. Both resemble an exit in the light reflected towards the ceiling. Wearable motion data, when present, could separate these cases. Falls are found more reliably than exits ({ET['fall']['det']} of {ET['fall']['n']} against {ET['exit']['det']} of {ET['exit']['n']}), because a fall moves the whole torso by more than a metre within a second, whereas a slow exit can resemble a posture change. The median delay from onset to the end of the first alarming window is {ET['fall']['lat_med']:.0f} s for falls and {ET['exit']['lat_med']:.0f} s for exits. Fused apnea alarms at $\\alpha=1$% catch {ET['apnea_onsets_detected']} of the {ET['apnea_onsets']} apnea onsets in the test half, and the first alarming 10 s window starts a median of {ET['apnea_latency_med']:.0f} s (90th percentile {ET['apnea_latency_p90']:.0f} s) after the first window that overlaps the apnea.",
f"Table XV gives paired Wilcoxon tests on the test episodes. Fusion increases coverage over the patch alone ($p_{{\\mathrm{{Holm}}}}$ {fp(HP[0])}, $r={S_COV[1]:.2f}$) and over the contactless channel alone ($p_{{\\mathrm{{Holm}}}}$ {fp(HP[1])}, $r={S_COVC[1]:.2f}$). It also raises the per-episode share within 2 bpm over the gated contactless estimate ($p_{{\\mathrm{{Holm}}}}$ {fp(HP[2])}, $r={S_ACC[1]:.2f}$).",
f"Generating one episode with the ray tracer takes {ST['t_episode_s']:.1f} s. Gating one 8 min episode of 36 pairs takes {ST['t_gate_s']:.2f} s and estimating its {ST['windows']} respiratory-rate windows takes {ST['t_rr_s']:.2f} s on one CPU core, more than 700 times faster than real time. The per-window eigen-decomposition involves a $36\\times36$ matrix at about 100 frequencies, which a luminaire controller can compute.",
    ]
    return out


DISC = {
"Principal Findings": [
"Ceiling luminaires can sense breathing while they continue to light the ward and carry data, but only in part of the operating space. With mid-level dimming at night and about 10% low-band pilots, contactless estimates fall within 2 bpm in roughly two of three windows and a 150 Mb/s downlink remains. With deep dimming or daylight the contactless channel weakens sharply. Its main value is coverage: it reports when the patch is absent or disconnected, and quality gating keeps the fused stream accurate.",
"Two results limit what contactless sensing can promise. Apnea could not be detected from ceiling reflections at these signal levels, and contactless accuracy depends strongly on posture. A patient lying on the side presents little reflecting area to the ceiling. Clinical claims for this channel should therefore be restricted to respiratory-rate coverage and gross events.",
],
"Implications for Design": [
"Fig. 8 turns Proposition 1 into design curves. The optimal pilot share rises from about 0.05 in a dark ward to 0.3 in bright daylight and barely depends on dimming. The minimum CRB is lowest at half dimming and grows steeply towards either end of the dimming range, while the downlink rate falls in steps as the bit loading drops. Table XVI condenses the results into a design guide. The pilot placement result is general. When receiver noise grows with frequency, sensing pilots belong at the bottom of the band, and their share has an optimum that can be computed from the receiver model before deployment. Lighting policy and monitoring quality are coupled through the AC swing. A ward that dims to 20% at night gives up most of its contactless sensing, while one that dims to 50% keeps it. The federated result is practical: sharing small histogram sketches across wards gave more uniform false-alarm control than calibrating each ward alone, without moving any patient data.",
],
"Comparison with Other Contactless Modalities": [
"Doppler radar resolves chest motion through phase and reaches sub-millimetre sensitivity [@li13]. Ceiling light offers only intensity, so its sensitivity to a few millimetres of excursion at 2 m is low, as Section III-A predicts. Its advantages lie elsewhere. It adds no emitter, because the luminaire is already there. It adds no radio energy near medical equipment, and it reveals no images of the patient, unlike camera-based monitoring. The results suggest a division of labour. Ceiling sensing provides coverage and gross-event detection, wearables provide precision and apnea detection, and radio sensing is reserved for beds where contact is impossible and radio emission is acceptable.",
],
"Deployment Considerations": [
"A deployment would place a 4 cm$^2$ photodiode with a low-noise sensing front end beside each luminaire above a bed. The low-band pilot share would follow the ambient light level, which a luminaire with a daylight sensor already measures, and the controller would use the nearest one or two luminaires instead of every visible one. Quantile sketches of alarm scores are small and can be exchanged between ward controllers or through a hospital server, without moving any physiological data. Alarm thresholds would be recalibrated whenever the lighting schedule of a ward changes, because the noise level and hence the score distribution change with it.",
],
"Clinical Integration and Privacy": [
"A ward deploying LUCID-Care would see three kinds of output on its central station. The first is a respiratory-rate stream with an explicit source label (patch, ceiling or both), so staff know when a value comes from the less precise channel. The second is apnea and desaturation alarms that come from the patch. The third is bed-exit and fall alerts from the ceiling at about one false alert per hour per bed at the chosen operating point. That rate is still high for a busy ward and would need the motion data of the patch or a bed sensor to fall further. Because the ceiling photodiodes measure only scalar gains, no image of the patient is formed, and the system raises fewer privacy concerns than camera-based fall detection. The federated calibration keeps raw scores in each ward, which fits hospital data-governance rules.",
],
"Threats to Validity": [
"The sensing channel is ray-traced with first-order body reflections, and the torso is a set of Lambertian patches. Real clothing, blankets and specular components will change the gains. The respiration waveforms are real, but the mapping from impedance to wall excursion uses a range taken from the literature rather than measured on the same patients. Falls and bed exits come from a kinematic model, and apneas were inserted, because the camera-only mirror of UP-Fall [@martinez19] lacked the sensor streams and the BIDMC recordings contain no annotated apneas. Calibration and test episodes share the 53 recordings across different wards and placements, so the respiration of one patient can appear in both halves.",
],
}

LIMITS = ["The ambient-light model uses two fixed photocurrent levels; real wards see daylight that varies by the minute and by window orientation, which would move the optimal pilot share over the day. The receiver model follows the Komine–Nakagawa formulation of the companion link study, and a dedicated low-noise sensing front end could change the thermal terms and hence the optimum. No bench-top or clinical measurement was possible in this study, so every sensing result is simulation-based and should be confirmed in a ward mock-up with real luminaires and photodiodes. The wearable channel is a noisy copy of the reference respiration with modelled outages, not a separate device. The event detector was trained on 53 episodes and calibrated on 53, which limits the precision of its false-alarm estimate. Federated calibration was evaluated with 16 wards and one split; the advantage of pooling may shrink when wards have more calibration data or larger systematic differences."]


def CONCL():
    b = pa(0.5, 0.1, 5e-6); q = AN['fusion_q']; al = AN['alarms']['fused|0.01']
    ro = [f"The three research objectives can now be answered. For RO1, sensing costs a ceiling luminaire about 10% of its subcarriers at night and 20% by day, which lowers the downlink by 9–19%; dimming below about 30% removes most of the sensing benefit. For RO2, the contactless channel adds about 11 points of rate coverage with a patch present 85% of the time, and up to 20 points when the patch is often missing, while the fused stream keeps about 85–90% of its estimates within 2 bpm. For RO3, pooled federated calibration kept the worst ward within a few percent of a 1% target, where local calibration left some wards above 10%."]
    return ro + [f"LUCID-Care adds respiratory sensing to the luminaires that already light and connect a hospital ward. Low-band sensing pilots give the channel-estimation bound an interior optimum, and at night with mid-level dimming the system keeps a {b['rate'] / 1e6:.0f} Mb/s downlink while placing {pc(b['within2'])}% of contactless respiratory-rate estimates within 2 bpm. Quality-gated fusion with an intermittent patch reported rates for {pc(q['fused']['coverage'])}% of windows instead of {pc(q['wearable']['coverage'])}%. Federated pooled calibration kept the worst ward at {100 * al['pooled']['far_max']:.1f}% false alarms for a 1% target. The contactless channel could not detect apneas, and its accuracy depends on posture. Future work should measure the body channel in a ward mock-up, select pairs adaptively, and add wearable motion to the bed-exit detector."]


def build():
    d = D.Doc(TEMPLATE)
    d.title(P.TITLE); d.authors(P.AUTHORS); d.affil(P.AFFIL); d.abstract(ABSTRACT())
    d.keywords('Visible light communication, integrated sensing and communication, contactless respiration monitoring, DCO-OFDM, '
               'Cramér–Rao bound, wearable fusion, conformal prediction, federated calibration')
    d.end_front()
    d.section('Introduction')
    for p in P.INTRO: d.para(T(p))
    d.subsection('Research Objectives')
    for p in ["**RO1.** Quantify how much downlink rate a ceiling luminaire must give up to estimate respiration from body-reflected pilots, under dimming and ambient light.",
              "**RO2.** Determine how much monitoring coverage and accuracy a contactless channel adds to an intermittent wearable.",
              "**RO3.** Set alarm thresholds that keep a target false-alarm rate in every ward when wards differ and calibration data are scarce."]:
        d.para(T(p))
    d.subsection('Contributions')
    for p in P.CONTRIB: d.para(T(p))
    d.para('Section II reviews related work, Section III presents the system, Section IV the data and protocol, Section V the results, Section VI the discussion, and Sections VII and VIII the limitations and conclusions.')
    d.section('Related Work')
    for sub, ps in P.RELATED.items():
        d.subsection(sub[3:])
        for p in ps: d.para(T(p))
    d.table([[T(c) for c in r] for r in P.RW_TABLE], 'Comparison with the closest sensing, ISAC and calibration studies', widths=[0.8, 1.4, 1.5, 1.0, 0.9, 1.4], wide=True, fs=7.2)
    d.section('LUCID-Care System')
    for k in 'ABCDE':
        name, items = P.METHOD[k]; d.subsection(name); add(d, items)
        if k == 'C':
            d.algorithm('Algorithm 1: Contactless respiration per 32 s window', [
                'Input: per-frame estimates $\\hat{h}_{\\ell k,t}$ of all pairs, CRB per photodiode',
                'Detrend each pair (20 s median), low-pass at 0.3 Hz, scale by MAD; compute $\\chi_t$ by (9); flag frames with $\\chi_t>6$ and dilate by 2 s',
                'if more than 5% of the window is flagged then withhold the window',
                'Band-pass 0.09–0.8 Hz; divide each pair by its CRB standard deviation',
                'Form $\\mathbf{S}(f)$ from 16 s Hann segments with 4 s hop by (10); take $\\lambda(f)$',
                'Output $\\widehat{\\mathrm{RR}}$ and quality $q_c$ by (11); apnea score by (12) on 10 s windows'])
        if k == 'D':
            d.algorithm('Algorithm 2: Federated Mondrian calibration of alarm thresholds', [
                'Each ward $w$ collects non-event calibration scores per modality $m$ (wearable, contactless)',
                'Each ward sends a histogram sketch of its scores (pooled variant) or its standardised scores with $(\\mu_{w,m},\\varsigma_{w,m})$ (normalised variant)',
                'Server merges sketches per modality and returns the $\\lceil (n+1)(1-\\alpha)\\rceil$-th order statistic',
                'Each ward sets $\\tau_{w,m}$ by (14) or (16)',
                'At run time: use the wearable score when the patch is available, otherwise the contactless score; alarm if the score exceeds $\\tau_{w,m}$'])
        if k == 'A':
            d.table(P.NOTATION, 'Principal notation', widths=[1.0, 2.6], fs=7.2)
            d.figure(F + 'fig1_arch.png', 'Architecture of LUCID-Care. Every luminaire emits illumination, data and low-band sensing pilots within the AC range allowed by its dimming level. Ceiling photodiodes estimate the luminaire–body–photodiode gains, which are gated and combined into contactless respiration, apnea and event scores. The ward controller fuses them with the wearable patch and calibrates per-ward, per-modality alarm thresholds from federated quantile sketches.', wide=True)
    d.section('Data and Experimental Protocol')
    d.subsection('Sensing Corpus')
    d.para(T("Episodes reuse the procedurally generated wards of the companion hybrid-access study. Each episode places one patient in a bed under a ward ceiling and traces the gains between the six nearest luminaires and the six photodiodes co-located with them. The patient lies still, with a trunk drift below 0.2°, and changes posture with an exponential dwell time of mean 300 s. The respiration of one BIDMC record [@pimentel17], band-passed to 0.05–1 Hz to remove electrode drift, drives the torso patches through (3). One or two apneas of 10–25 s are inserted. With probability 0.3 the patient sits up and stands beside the bed, and half of those patients then fall to the floor within 0.6–1.2 s. Staff stand at the bedside up to three times and block the paths between the torso and the luminaires on their side. Table III summarises the corpus."))
    d.table(CORPUS(), 'Sensing corpus', widths=[2.2, 1.6], fs=7.2)
    d.subsection('Reference Respiratory Rate and Metrics')
    d.para(T("The reference rate of a window comes from the manual breath annotations $b_1<\\dots<b_K$ that fall inside it:"))
    d.eq('\\mathrm{RR}_{\\mathrm{ref}}=\\frac{60\\,(K-1)}{b_K-b_1}')
    d.para(T("Windows containing an apnea or any frame out of bed are excluded from respiratory-rate scoring. Metrics are coverage (the share of windows with an estimate), the share within 2 bpm and the mean absolute error. Apnea windows are those with at least 80% apnea frames, and non-event windows have none. The false-alarm rate is computed per ward. Events count as detected if any alarm overlaps the 10 s after onset. To convert a target of $A_h$ false alarms per hour into a window level for the event detector, the overlap of windows (hop 2 s, so five windows per 10 s) is accounted for:"))
    d.eq('\\alpha_{\\mathrm{ev}}=\\frac{5\\,A_h\\,H_{\\mathrm{cal}}}{N_{\\mathrm{cal}}}')
    d.para(T("Here $H_{\\mathrm{cal}}$ is the number of calibration hours and $N_{\\mathrm{cal}}$ the number of calibration non-event windows. Episodes with even indices form the calibration half and odd indices the test half, so every ward appears in both. Wards differ in operating condition: five combinations of dimming (0.2, 0.35, 0.5) and ambient light (night, day) are assigned across the 16 layouts. Paired Wilcoxon signed-rank tests over test episodes use Holm correction and effect size $r=|Z|/\\sqrt{n}$ [@wilcoxon45, holm79]."))
    d.para(T("The wearable channel is modelled as the reference impedance respiration plus white noise at 10 dB SNR, sampled at the sensing rate. Its availability alternates between exponential up periods (mean 300 s) and down periods (mean 20 s), reflecting the uplink outages and brief detachments observed in the companion link study, and 15% of episodes have no patch at all. This availability model is deliberately simple, and Section V-G varies it. The contactless and wearable channels are processed with the same estimator, so differences between them reflect signal quality, not processing."))
    d.table(PARAMS, 'System and protocol parameters', widths=[2.2, 1.6], fs=7.2)
    d.section('Results')
    R = RES()
    d.subsection('Rate–CRB–Accuracy Frontier')
    for p in R['A']: d.para(T(p))
    d.figure(F + 'fig2_signal.png', 'Episode 17 at night with $d=0.5$ and $\\rho=0.1$. Top: noise-free body-path gain of the strongest pair, showing a posture change near 180 s and a later change near 310 s. Middle: BIDMC reference respiration and the contactless eigen-combined signal; shaded spans are apneas (dark) and staff occlusions (light). Bottom: contactless rate per 32 s window against the annotated rate.', wide=True)
    d.figure(F + 'fig3_pareto.png', 'Rate–CRB–accuracy frontier on 53 episodes. (a) Share of windows within 2 bpm against downlink rate as the pilot fraction varies from 0.4 (left end) to 0.01 (right end). (b) CRB of the per-frame gain estimate against the pilot fraction; black: night, grey: day; line style: dimming level as in (a).', wide=True)
    d.table(PARETO_TABLE(), 'Rate, CRB and contactless accuracy over the full operating grid (53 episodes per point)', widths=[0.6, 0.4, 0.4, 0.7, 0.8, 0.9, 0.7, 0.8], wide=True, fs=7.2)
    d.subsection('Wearable Fusion')
    for p in R['B']: d.para(T(p))
    d.table(FUSION_TABLE(), 'Respiratory-rate coverage (%) / within 2 bpm (%) / MAE (bpm) on the test half', widths=[1.5, 1.4, 1.4, 1.4], wide=True, fs=7.2)
    d.figure(F + 'fig4_fusion_alarms.png', '(a) Coverage against accuracy as the quality threshold varies, for the contactless and wearable estimators, and the fused operating point. (b) Per-ward false-alarm rates of fused apnea alarms at $\\alpha=1\\%$ (dotted line) for four calibration methods. (c) Distributions of the apnea score in apnea and normal windows for the wearable and contactless channels (test half).', wide=True)
    d.subsection('Ablation Study')
    for p in R['C']: d.para(T(p))
    d.table(ABL_TABLE(), 'Ablation of the contactless estimator (night, $d=0.5$, $\\rho=0.1$, 53 episodes)', widths=[2.6, 0.8, 0.9, 0.8], wide=True, fs=7.2)
    d.subsection('Design Sensitivity')
    for p in RES_DESIGN(): d.para(T(p))
    d.table(DESIGN_TABLE(), 'Design sensitivity of the contactless estimator (night, $d=0.5$, $\\rho=0.1$, 53 episodes)', widths=[2.0, 0.8, 0.9, 0.8], wide=True, fs=7.2)
    d.figure(F + 'fig5_design.png', 'Contactless accuracy (night, $d=0.5$, $\\rho=0.1$) by (a) quartile of the strongest body-path gain, a proxy for posture, and (b) reference respiratory rate; (c) effect of photodiode area, of the number of luminaires and photodiodes per side used, and of window length.', wide=True)
    d.subsection('Alarm Calibration across Wards')
    for p in R['D']: d.para(T(p))
    d.table(ALARM_TABLE(), 'Apnea alarm calibration on the test half: false-alarm rate (FAR) per ward and sensitivity', widths=[1.1, 1.3, 0.8, 0.9, 0.9, 0.9], wide=True, fs=7.2)
    for p in RES_BUDGET(): d.para(T(p))
    d.table(BUDGET_TABLE(), 'Worst-ward false-alarm rate (%) / share of wards above $2\\alpha$ (%) of fused apnea alarms at $\\alpha=1\\%$ for reduced calibration budgets (mean of five random subsamples)', widths=[1.4, 1.3, 1.3, 1.3], wide=True, fs=7.2)
    d.figure(F + 'fig6_calib.png', '(a) Worst-ward false-alarm rate against the target $\\alpha$ (dotted: equality). (b) Worst-ward false-alarm rate at $\\alpha=1\\%$ against the share of calibration episodes used. (c) Bed-exit and fall detection against false alarms per hour as the hourly target (labels) varies.', wide=True)
    d.subsection('Per-Ward Results and Link Budget')
    for p in RES_WARD(): d.para(T(p))
    d.table(WARD_TABLE(), 'Per-ward results on the test half (fused respiratory rate; fused apnea alarms at $\\alpha=1\\%$)', widths=[0.5, 0.9, 0.8, 0.9, 0.8, 0.8, 0.8], wide=True, fs=7.0)
    d.figure(F + 'fig7_snr.png', 'Contactless accuracy per test episode against the coherent-combining SNR bound of the respiratory component (night, $d=0.5$, $\\rho=0.1$).')
    d.subsection('Generalization/External Validation')
    for p in RES_EXT(): d.para(T(p))
    d.table(EXT_TABLE(), 'External checks: thresholds, quality gates and detectors calibrated on the main corpus and applied to stress-test wards, and a record-disjoint split of the main corpus', widths=[1.6, 1.0, 1.0, 1.0], wide=True, fs=7.2)
    for p in RES_WEAR(): d.para(T(p))
    d.table(WEAR_TABLE(), 'Respiratory-rate coverage / within 2 bpm (%) and fused apnea sensitivity under different wearable availability', widths=[1.6, 0.9, 0.9, 0.9, 0.9], wide=True, fs=7.2)
    d.subsection('Error Analysis by Respiratory-Rate Band')
    for p_ in RES_BAND(): d.para(T(p_))
    d.table(BAND_TABLE(), 'Accepted-window coverage (%) / within 2 bpm (%) / mean bias (bpm) by reference respiratory-rate band on the test half', widths=[1.1, 0.5, 1.2, 1.2, 1.2], wide=True, fs=7.2)
    d.subsection('Events, Statistics and Efficiency')
    for p in R['E']: d.para(T(p))
    d.table(STAT_TABLE(), 'Paired Wilcoxon signed-rank tests on test episodes', widths=[2.3, 0.4, 0.8, 0.5], fs=7.2)
    d.section('Discussion')
    for sub, ps in DISC.items():
        d.subsection(sub)
        for p in ps: d.para(T(p))
    d.table([['Ward condition', 'Pilot share', 'Pairs used', 'Expected contactless accuracy', 'Downlink rate'],
             ['Night, $d=0.5$', '$\\rho\\approx0.1$', 'nearest 1–2', '≈60–65% within 2 bpm', '≈150 Mb/s'],
             ['Night, $d=0.2$', '$\\rho\\approx0.1$', 'nearest 1–2', '≈45% within 2 bpm', '≈90 Mb/s'],
             ['Day, $d=0.5$', '$\\rho\\approx0.2$', 'nearest 1–2', '≈55–58% within 2 bpm', '≈135 Mb/s'],
             ['Day, $d=0.2$', 'sensing not advised', '—', '≈35% within 2 bpm', '≈80 Mb/s']],
            'Design guide derived from Table V, Proposition 1 and Table VIII', widths=[1.2, 1.0, 0.8, 1.3, 0.8], wide=True, fs=7.2)
    d.figure(F + 'fig8_analytic.png', 'Analytic design curves from Proposition 1. (a) CRB-optimal pilot fraction against ambient light for three dimming levels. (b) Minimum CRB at the optimal pilot fraction (left axis) and downlink rate at $\\rho=0.1$ (right axis) against dimming.', wide=True)
    d.section('Limitations')
    for p in LIMITS: d.para(T(p))
    d.section('Conclusion')
    for p in CONCL(): d.para(T(p))
    d.unnumbered('Appendix A: Derivation of the CRB (6)')
    d.para(T("Consider one luminaire–photodiode pair and one frame. The pilot band carries $M=\\rho BT_f$ real samples. In each sample the photocurrent is $i_m=R\\,h\\,u_m+n_m$, where $u_m$ is the known pilot waveform with mean power $\\sigma(d)^2$ and $n_m$ is Gaussian noise of variance $N_k(\\rho)\\rho B$. The log-likelihood in $h$ is quadratic, and its Fisher information is"))
    d.eq('\\mathcal{I}(h)=\\sum_{m=1}^{M}\\frac{R^{2}u_m^{2}}{N_k(\\rho)\\rho B}=\\frac{R^{2}\\sigma(d)^{2}\\rho BT_f}{N_k(\\rho)\\rho B}')
    d.para(T("Inverting it gives $\\mathrm{CRB}=N_k(\\rho)/(R^{2}\\sigma(d)^{2}T_f)$ per unit pilot share. Since only the fraction $\\rho$ of the AC power is spent on pilots, the pilot power is $\\rho\\sigma(d)^2$, and this yields (6). The least-squares estimator $\\hat{h}=\\sum_m u_m i_m/(R\\sum_m u_m^2)$ is unbiased and attains the bound. The bound in Proposition 2 follows from the coverage of the pooled quantile under the mixture $\\bar{P}$ and from the fact that the probability of any event under $P_w$ and $\\bar{P}$ differs by at most their total-variation distance."))
    d.unnumbered('Appendix B: Ward-Level Check of Proposition 2')
    d.para(T("For a threshold $\\tau$, the false-alarm rate of ward $w$ is $P_w(s>\\tau)$, so only one event of the score space enters the comparison between pooled and local calibration. With the pooled conformal threshold $\\bar{\\tau}$, the heterogeneity term of the bound is therefore the tail discrepancy"))
    d.eq('\\delta_w=\\left|P_w\\left(s>\\bar{\\tau}\\right)-\\alpha\\right|\\le d_{\\mathrm{TV}}\\left(P_w,\\bar{P}\\right)')
    d.para(T("which is never larger than the total-variation distance and is usually much smaller. The local threshold, set from $n_w$ calibration windows, has a sampling error whose standard deviation is approximately"))
    d.eq('s_w=\\sqrt{\\frac{\\alpha\\left(1-\\alpha\\right)}{n_w/\\kappa}}')
    d.para(T(("where $\\kappa=5$ accounts for the fivefold overlap of the 10 s apnea windows taken with a 2 s hop. Pooling is preferred in ward $w$ when $\\delta_w<s_w$. Table XVII lists both quantities on the calibration half at $\\alpha=1$%%. Measured with 32 quantile bins, the full total-variation distance (median %.3f for contactless and %.3f for wearable scores) exceeds the DKW sampling term $\\sqrt{\\ln(2/\\eta)\\kappa/(2n_w)}$ at $\\eta=0.1$ in %d of 16 and %d of 16 wards, so the condition of Proposition 2 in its total-variation form is not met. The tail form is met in %d of 16 wards for contactless scores and %d of 16 for wearable scores. For wearable scores, the three smallest wards (1171, 1188 and 1201, each with fewer than 300 calibration windows) all favour pooling, because there the local 99th percentile rests on two or three effective exceedances. For contactless scores, two of them do, while ward 1201 has a large tail discrepancy that neither method handles well. The small wards are also the ones that set the worst-ward false-alarm rate of local calibration in Table IX. In the larger wards a local threshold is the better choice, which suggests a hybrid rule that pools only below a ward-size threshold. The rule follows directly from the table, but it was not evaluated on the test half and is left for future work.") % (APX_TVMED['contactless'], APX_TVMED['wearable'], APX_TVOVER['contactless'], APX_TVOVER['wearable'], APX_POOL['contactless'], APX_POOL['wearable'])))
    d.table(APX_TABLE(), 'Tail discrepancy $\\delta_w$ of the pooled threshold and sampling deviation $s_w$ of the local threshold on calibration episodes ($\\alpha=1\\%$; values in units of $10^{-3}$)', widths=[0.6, 0.6, 1.0, 1.0, 0.9], fs=7.0)
    d.unnumbered('Data Availability')
    d.para('The BIDMC PPG and Respiration Dataset is openly available from PhysioNet (Open Data Commons Attribution License). The CMU motion-capture clips are public. All sensing episodes and results are regenerated deterministically by the released code from the listed seeds.')
    d.unnumbered('Code Availability')
    d.para('The source code (ward geometry and ray tracing, sensing episodes, dual-function waveform model, estimators, fusion, federated calibration, statistics and figures) is provided as a supplementary archive with a run script. [Repository URL and archived release DOI to be added upon acceptance.]')
    d.unnumbered('Ethics Statement')
    d.para('Only de-identified, openly licensed recordings were used, and no new human-subject data were collected.')
    d.unnumbered('Author Contributions')
    d.para('[To be completed by the authors using the CRediT taxonomy.]')
    d.unnumbered('Conflict of Interest')
    d.para('The authors declare no competing interests.')
    d.unnumbered('Acknowledgment')
    d.para('[Funding information to be added by the authors.]')
    d.unnumbered('References')
    d.references(C.bibliography())
    d.render(OUT)


if __name__ == '__main__':
    import os
    os.makedirs('./out', exist_ok=True)
    build()
    print('refs', len(C.ORDER), 'of', len(C.REFS))
