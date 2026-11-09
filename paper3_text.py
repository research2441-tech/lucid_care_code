"""Static manuscript content for Paper 3 (LUCID-Care)."""

TITLE = ('LUCID-Care: Dual-Function Visible-Light Sensing and Communication for Contactless-Plus-Wearable '
         'Patient Monitoring with Federated Conformal Alarm Calibration across Hospital Wards')
AUTHORS = '[First Author], [Second Author], [Third Author], and [Corresponding Author]'
AFFIL = '[Department, Institution, City, Country]  |  Corresponding author: [email]'

INTRO = [
"Respiratory rate is the vital sign that changes earliest before clinical deterioration, and it is also the one measured least reliably on general wards. Track-and-trigger systems such as NEWS2 assign it the largest score range [@rcp17, smith13]. Yet manual counts are infrequent, and continuous wearable monitoring is still mostly at the feasibility stage [@leenen20]. A chest patch measures breathing well while it is attached, charged and connected. When the patch is removed, its battery is flat or its uplink is shadowed, the bed becomes unmonitored.",
"The ceiling of the ward offers a second sensor that needs no attachment. Each luminaire already illuminates the patient below it. If the luminaire also carries a modulated downlink, as in optical medical body sensor networks [@torkestani12, zwaag21, antaki25], a photodiode next to it can estimate the channel through the patient's body. The light reflected from the chest and abdomen changes with every breath, with every turn in bed and when the patient leaves the bed. Visible-light sensing of human posture has been demonstrated with ceiling photodiodes [@li15], radio sensing of breathing is well established [@li13], and integrated sensing and communication (ISAC) is now a design paradigm for 6G systems [@liu22, xiong23].",
"Turning the luminaire into a clinical sensor raises three questions that earlier work does not answer. The first is how much of the lighting and data function must be given up for sensing. The waveform must respect illumination, dimming and flicker constraints [@rajagopal12], and every pilot used to estimate the body channel is a subcarrier taken from data. The second is how a contactless estimate that is reliable only some of the time should be combined with a wearable estimate that is missing some of the time. The third is how to set alarm thresholds that hold a promised false-alarm rate in every ward, when wards differ in layout, lighting and staff traffic, and when each ward alone has too little data to calibrate reliably. Monitor alarms are already a major burden [@drew14].",
"The gap can be stated precisely. No existing system (i) sizes the sensing pilots of a dual-function luminaire against its data rate under the illumination and dimming constraints of a ward, (ii) quantifies how much monitoring coverage an ambient optical channel adds to an intermittent wearable, and (iii) gives a per-ward false-alarm promise for the resulting alarms when each ward has little calibration data. Each of these is needed before a ward can rely on its lighting as part of its monitoring.",
"This paper addresses the three questions with a system called LUCID-Care. Its contributions are as follows.",
]

CONTRIB = [
"**A dual-function DCO-OFDM luminaire with a low-band pilot design.** The body channel is estimated from pilots placed in the lowest part of the band. In this region the receiver's frequency-dependent thermal noise is smallest, so the Cramér–Rao bound (CRB) of the estimate has an interior minimum in the pilot fraction. This yields a rate–CRB–accuracy frontier that depends on dimming and ambient light.",
"**A ray-traced sensing corpus driven by clinical respiration.** The corpus contains 212 eight-minute episodes. In each, the impedance respiration of the 53 BIDMC recordings [@pimentel17] moves the torso of a patient lying under a ray-traced ward ceiling. The episodes also contain inserted apneas, bed exits, falls, posture changes and staff occlusions.",
"**Gated, noise-whitened spectral eigen-combining for contactless respiratory rate, with quality-gated fusion.** The fusion combines the contactless estimate with a wearable patch whose availability is intermittent.",
"**Federated conformal alarm calibration.** Wards share only quantile sketches, and thresholds are kept per ward and per modality (Mondrian). Per-ward false-alarm rates are compared against local and fixed thresholds.",
]

RELATED = {
"A. Optical Wireless Links in Hospitals": [
"Indoor optical channels combine a Lambertian line of sight with diffuse reflections [@kahn97, barry93], and site-specific models are standard for VLC [@komine04, miramirkhani15]. DCO-OFDM carries data on a biased LED drive, and its clipping noise limits the usable modulation depth [@armstrong09, dimitrov12]. Hospitals have adopted optical links to keep radio emission away from medical devices [@torkestani12, lebas18, zwaag21], and surveys position VLC as an e-health access technology [@guana24, celik22, haas16]. Body sensor networks and their channel models have been reviewed extensively [@movassaghi14]. Orientation and human blockage dominate optical body links [@soltani19, pal24], and adaptive modulation for an optical medical body sensor network has been studied with learned channel estimators [@antaki25]. Hybrid radio–optical scheduling and load balancing address coverage gaps [@basnayaka17, wu20, hamrouni25]. None of these works uses the downlink illumination itself to observe the patient.",
"Two properties of these links matter for sensing. First, the received power on an optical link is a smooth function of geometry. A displacement of the reflecting surface by a fraction of a wavelength has no effect, and only millimetre-scale changes of distance and angle alter the gain. This differs from radio sensing, where phase is sensitive to sub-millimetre motion. Second, the luminaire's drive is shared between illumination and modulation, so every sensing pilot competes with data for the same limited AC swing [@armstrong09, dimitrov12]. Both properties make optical sensing of breathing a low-SNR problem that has to be engineered together with the data link, not added afterwards.",
],
"B. Contactless Sensing and ISAC": [
"Doppler radar measures chest-wall motion without contact and has been reviewed for healthcare monitoring [@li13]. Visible-light sensing reconstructs coarse human posture from shadows cast on ceiling or floor photodiodes [@li15]. ISAC theory characterises the trade-off between communication rate and sensing accuracy as a Pareto frontier bounded by the CRB [@liu22, xiong23]. In VLC, the trade-off is constrained by lighting: dimming support and flicker-free modulation are part of the IEEE 802.15.7 standard [@rajagopal12], and the AC swing available for both data and pilots shrinks as the light is dimmed. Estimating spectral peaks from noisy, gapped signals relies on averaged periodograms [@welch67]. No earlier study has quantified how ceiling luminaires trade data rate against respiratory-rate accuracy under dimming, ambient light and pilot placement.",
"Radio ISAC and optical ISAC differ in what limits sensing. In radio, the transmit power budget and clutter dominate, and the sensing target is a delay or Doppler parameter. In VLC with intensity modulation, the target is a slowly varying gain. The noise is set by the receiver front end and by ambient light, and the waveform must keep a fixed mean for illumination. Pilot placement inside the band therefore becomes a design variable with no direct radio counterpart, because the thermal noise of a transimpedance receiver grows with frequency [@komine04]. Respiratory-rate estimation from noisy multichannel signals also depends on how channels are combined and on how gaps are handled, which is why the estimator is treated as part of the system design here.",
],
"C. Respiration, Falls and Clinical Alarms": [
"The BIDMC dataset provides photoplethysmogram, ECG and impedance respiration with manual breath annotations from 53 ICU recordings and is a standard benchmark for respiratory-rate estimation [@pimentel17]. Fall detection and prevention systems span wearable, ambient and vision sensors [@chaccour17], and the UP-Fall dataset combines wearable, ambient and camera modalities [@martinez19]. Early warning scores weight respiratory rate heavily [@rcp17, smith13], and bedside monitors already generate alarm volumes that burden staff [@drew14]. Continuous wearable monitoring on wards remains mostly in validation studies [@leenen20]. Urgency-aware status updates have been formalised for networked control [@zheng20]. The combination of an ambient sensor and a wearable with explicit alarm-rate guarantees has not been studied.",
"Wearable respiratory monitoring fails in ways that matter clinically. Patches lose skin contact, batteries run out between changes, and radio or optical uplinks are interrupted by body shadowing and by staff. Each failure leaves a gap that is invisible on the central station unless it is explicitly reported. An ambient channel that never detaches but is less precise is therefore a natural complement, provided its estimates come with a quality measure that prevents a poor contactless estimate from overriding a good wearable one.",
],
"D. Conformal and Federated Calibration": [
"Conformal prediction gives distribution-free coverage from exchangeable calibration data [@vovk05, angelopoulos23b], and it has been applied to calibrate wireless predictors [@cohen23]. Online variants adapt under distribution shift [@gibbs21]. Federated conformal predictors aggregate calibration information across clients without sharing raw scores [@lu23], and label-shift-aware versions weight the clients [@plassier23]. Gradient-boosted trees are a strong tabular baseline for event classification [@friedman01]. Statistical comparisons use paired Wilcoxon tests with Holm correction [@wilcoxon45, holm79]. Walking patients and trajectories in earlier simulations come from the CMU motion-capture database [@cmu].",
"Alarm calibration in hospitals is a natural federated problem. Each ward holds its own monitoring data, privacy rules discourage pooling raw recordings, and wards differ in patient mix, lighting and staff routines. Federated conformal methods have so far been evaluated on image and tabular benchmarks with marginal coverage targets [@lu23, plassier23]. Whether a per-ward false-alarm promise can be kept with a few hundred calibration windows per ward, and how pooling compares with local calibration under realistic heterogeneity, has not been examined for monitoring alarms.",
],
}

RW_TABLE = [
["Work", "Sensor", "Signal used", "Communication", "Fusion", "Alarm guarantee"],
["[@li15]", "Ceiling/floor photodiodes", "Shadow pattern (posture)", "No", "No", "No"],
["[@li13]", "Doppler radar", "Chest motion", "No", "No", "No"],
["[@liu22]", "Radio ISAC (general)", "Echo / CRB", "Yes", "No", "No"],
["[@antaki25]", "Body-worn VLC", "Link channel", "Yes", "No", "No"],
["[@lu23]", "—", "Federated calibration", "—", "—", "Coverage, pooled"],
["LUCID-Care", "Ceiling photodiodes + patch", "Body-reflected DCO-OFDM pilots", "Yes (rate–CRB)", "Quality-gated", "Per-ward Mondrian, federated"],
]

NOTATION = [
["Symbol", "Meaning"],
["$h_{\\ell k,t}$", "luminaire $\\ell$ → body → photodiode $k$ gain in frame $t$"],
["$d,\\ \\rho$", "dimming level; pilot fraction of the band"],
["$\\sigma(d)$", "AC (OFDM) optical standard deviation"],
["$N(\\rho)$", "noise density over the pilot band"],
["$\\mathrm{CRB}_k$", "variance bound of the per-frame gain estimate"],
["$x(t)$", "normalised BIDMC respiration; $a$ tidal excursion"],
["$\\mathbf{S}(f)$", "whitened cross-spectral matrix of the pair signals"],
["$q_c,\\ q_w$", "quality of contactless and wearable estimates"],
["$\\alpha$", "target false-alarm rate per window"],
["$\\tau_{w,m}$", "alarm threshold of ward $w$, modality $m$"],
]

METHOD = {
'A': ("Dual-Function Luminaire and Sensing Channel", [
('p', "Each ceiling luminaire $\\ell$ emits an optical power $P_{\\ell}(t)=dP_{\\max}+s_{\\ell}(t)$, where $d\\in(0,1)$ is the dimming level set by the lighting system and $s_\\ell$ is a DCO-OFDM signal. To keep the drive inside the LED's linear range with clipping at $\\kappa$ standard deviations, the AC swing is limited by the nearer of the two range limits:"),
('eq', "\\sigma(d)=\\frac{\\min\\{d,\\,1-d\\}\\,P_{\\max}}{\\kappa}"),
('p', "Light reaches photodiode $k$ next to luminaire $k$ through the static room ($h^{s}_{\\ell k}$) and through the patient's torso, modelled as Lambertian patches $j$ with reflectance $\\rho_b$, area $A_j$, position $\\mathbf{p}_j$ and normal $\\mathbf{n}_j$. The first-order body gain is"),
('eq', "h_{\\ell k}=\\sum_{j}\\frac{(m+1)\\cos^{m}\\phi_{\\ell j}\\cos\\theta_{\\ell j}}{2\\pi r_{\\ell j}^{2}}\\cdot\\frac{\\rho_b A_j A_{\\mathrm{pd}}g_c\\cos\\theta_{jk}\\cos\\psi_{jk}}{\\pi r_{jk}^{2}}"),
('p', "Here $r_{\\ell j}$ and $r_{jk}$ are the two leg lengths, the cosines are the emission, incidence and detection angles, and the detection term vanishes outside the photodiode field of view. Respiration displaces every patch along the torso normal:"),
('eq', "\\mathbf{p}_j(t)=\\bar{\\mathbf{p}}_j(t)+a\\,\\beta_j\\,x(t)\\,\\mathbf{n}_j(t)"),
('p', "$x(t)$ is the band-limited BIDMC impedance waveform scaled to unit peak-to-peak amplitude, $a$ is the tidal excursion of the abdominal wall, and $\\beta_j$ tapers the excursion from abdomen to chest. Linearising (2) around the resting posture gives the respiratory modulation of the gain, which is of the order of $10^{-3}$ of $h$ for millimetre excursions at 2 m:"),
('eq', "\\Delta h_{\\ell k}(t)\\approx a\\,x(t)\\sum_j\\beta_j\\,\\mathbf{n}_j^{\\top}\\nabla_{\\mathbf{p}_j}h_{\\ell k}"),
]),
'B': ("Low-Band Pilots and the Rate–CRB Trade-off", [
('p', "Each luminaire reserves a fraction $\\rho$ of the subcarriers as known pilots, and the pilots of neighbouring luminaires are frequency-interleaved so that each photodiode can separate them. The pilots occupy the lowest band $[0,\\rho B]$, where the LED response is flattest and where the field-effect-transistor term of the receiver noise, which grows with frequency [@komine04], is smallest. Averaging the receiver noise of photodiode $k$ over the pilot band gives the density"),
('eq', "N_k(\\rho)=2q\\left(R\\,d P_{\\max}h^{\\mathrm{dc}}_{k}+I_{\\mathrm{amb}}\\right)+\\frac{8\\pi k_BT\\eta A I_2}{G}\\rho B+\\frac{16\\pi^{2}k_BT\\Gamma\\eta^{2}A^{2}I_3}{g_m}\\rho^{2}B^{2}"),
('p', "where $h^{\\mathrm{dc}}_k$ is the total illumination gain of all luminaires, $I_{\\mathrm{amb}}$ is the ambient-light photocurrent and the remaining constants are those of the Komine–Nakagawa receiver model. Pilot power is $\\rho\\sigma(d)^2$ spread over $\\rho B$, and a frame of length $T_f$ contains $\\rho BT_f$ pilot samples. The per-frame gain estimate $\\hat{h}_{\\ell k,t}=h_{\\ell k,t}+h^{s}_{\\ell k}+\\nu_{\\ell k,t}$ is therefore unbiased, with variance at the CRB:"),
('eq', "\\mathrm{CRB}_k(\\rho,d)=\\mathrm{Var}\\left[\\nu_{\\ell k,t}\\right]=\\frac{N_k(\\rho)}{R^{2}\\sigma(d)^{2}\\rho\\,T_f}"),
('p', "The bound falls as $1/\\rho$ while thermal noise is small, and rises again when the pilot band reaches frequencies where the $\\rho^2$ term dominates. It therefore has an interior minimum. The data rate uses the remaining subcarriers with the per-group bit loading of a 1e-3 symbol-error target:"),
('eq', "R_{\\mathrm{data}}(\\rho,d)=(1-\\rho)\\,B\\,\\frac{\\sum_{g}b_g(d)\\,|\\mathcal{K}_g|}{N_{\\mathrm{FFT}}+L_{\\mathrm{CP}}}"),
('p', "Sweeping $\\rho$ for each dimming and ambient condition traces the achievable rate–CRB–accuracy frontier. Formally, the luminaire designer solves, for a dimming level fixed by the lighting schedule,"),
('eq', "\\max_{\\rho\\in(0,1)}\\ R_{\\mathrm{data}}(\\rho,d)\\quad\\text{s.t.}\\quad \\max_{k}\\mathrm{CRB}_k(\\rho,d)\\le\\epsilon_h"),
('p', "where $\\epsilon_h$ is the estimation accuracy that the respiratory estimator needs. Because $R_{\\mathrm{data}}$ falls monotonically in $\\rho$ and the CRB is unimodal, the solution is the smallest $\\rho$ that meets the accuracy constraint on the decreasing branch of the CRB, or no feasible point if $\\epsilon_h$ is below the minimum CRB at that dimming level. In the latter case sensing must be paused or the lighting raised."),
]),
'C': ("Contactless Respiration, Apnea and Events", [
('p', "Two kinds of gross change corrupt the respiratory signal: posture changes, and staff standing between the luminaires and the bed. Each pair signal is detrended with a 20 s running median, low-passed at 0.3 Hz and divided by its median absolute deviation, $\\tilde{y}_{i,t}$. A frame is flagged when the pooled robust score"),
('eq', "\\chi_t=\\left(\\frac{1}{P}\\sum_{i=1}^{P}\\tilde{y}_{i,t}^{2}\\right)^{1/2}>6"),
('p', "is exceeded, and the flag is dilated by 2 s. Windows of 32 s with more than 5% flagged frames are withheld. On the remaining frames, the band-passed (0.09–0.8 Hz) pair signals are divided by their CRB standard deviation, so that noise is white and of equal power on every pair. Their cross-spectral matrix is estimated from overlapping 16 s Hann segments $\\mathbf{X}_s(f)$:"),
('eq', "\\mathbf{S}(f)=\\frac{1}{S}\\sum_{s=1}^{S}\\mathbf{X}_s(f)\\mathbf{X}_s(f)^{\\mathrm{H}},\\qquad \\lambda(f)=\\lambda_{\\max}\\left(\\mathbf{S}(f)\\right)"),
('p', "Respiration is a single waveform seen through different gains, so it occupies the dominant eigenvector, while whitened noise spreads evenly over all eigenvalues. The respiratory rate and the window quality are"),
('eq', "\\widehat{\\mathrm{RR}}=60\\,\\arg\\max_{f}\\lambda(f),\\qquad q_c=\\frac{\\max_f\\lambda(f)}{\\mathrm{median}_f\\,\\lambda(f)}"),
('p', "An apnea score compares the coherent respiratory-band power $\\bar{\\lambda}_t$ of a 10 s window, the largest eigenvalue of the whitened pair covariance, with the median $\\bar{\\lambda}_{\\mathrm{base}}$ of the preceding 120 s:"),
('eq', "s_{\\mathrm{ap},t}=\\log\\bar{\\lambda}_{\\mathrm{base},t}-\\log\\bar{\\lambda}_{t},\\qquad \\bar{\\lambda}_{\\mathrm{base},t}=\\underset{t-120\\le t'<t}{\\mathrm{median}}\\,\\bar{\\lambda}_{t'}"),
('p', "Bed exits and falls are detected by a gradient-boosted classifier [@friedman01]. Its inputs are the vector of slow pair-level changes between the 10 s before and after the window centre, the norm of that vector and the steepest 1 s slope, each normalised by the CRB."),
]),
'D': ("Wearable Fusion and Federated Conformal Alarms", [
('p', "The wearable patch produces an estimate $\\widehat{\\mathrm{RR}}_w$ with quality $q_w$ computed in the same way from its own impedance channel, when it is present and its uplink delivers more than 95% of the window. Quality thresholds $t_c$ and $t_w$ are chosen on calibration episodes as the smallest values at which 90% of the accepted estimates fall within 2 bpm. The fused estimate takes the better-qualified source:"),
('eq', "\\widehat{\\mathrm{RR}}=\\begin{cases}\\widehat{\\mathrm{RR}}_c,& q_c\\ge t_c\\ \\text{and}\\ (q_w<t_w\\ \\text{or}\\ q_c/t_c\\ge q_w/t_w),\\\\ \\widehat{\\mathrm{RR}}_w,& q_w\\ge t_w\\ \\text{otherwise},\\end{cases}"),
('p', "and no estimate is reported when neither passes. Alarm thresholds come from split conformal calibration [@vovk05]. For a ward $w$ and modality $m$ (Mondrian partition), the non-event calibration scores $s_1,\\dots,s_n$ give the threshold"),
('eq', "\\tau_{w,m}=s_{(\\lceil (n+1)(1-\\alpha)\\rceil)}"),
('p', "which, for exchangeable scores, guarantees"),
('eq', "1-\\alpha\\le\\Pr\\left\\{s_{\\mathrm{new}}\\le\\tau_{w,m}\\right\\}\\le1-\\alpha+\\frac{1}{n+1}"),
('p', "With a few hundred non-event windows per ward, the local quantile is noisy, and the guarantee holds only on average over calibration sets. Federated calibration pools information across wards without moving raw data [@lu23]. In the pooled variant each ward sends a histogram sketch of its scores, and the server returns the pooled $(1-\\alpha)$ quantile. In the normalised variant each ward also sends its median $\\mu_w$ and scale $\\varsigma_w$. The server computes the quantile $q_z$ of the standardised scores, and each ward sets"),
('eq', "\\tau_{w,m}=\\mu_{w,m}+\\varsigma_{w,m}\\,q_z,\\qquad z=\\frac{s-\\mu_{w,m}}{\\varsigma_{w,m}}"),
('p', "which adapts the threshold to the ward's noise level, as in label-shift-aware federated conformal prediction [@plassier23]. The event detector is calibrated the same way, with the window-level level $\\alpha$ derived from a target number of false alarms per hour."),
('p', "The communication cost of federated calibration is small. A histogram sketch with $J$ bins per modality costs each ward $2J$ counts per calibration round, independent of the number of windows, and the server returns one threshold per ward and modality. With $J=512$ bins of width $\\Delta$ on the score axis, the quantile is resolved to $\\Delta$, and the extra false-alarm probability caused by binning is bounded by the score mass of one bin:"),
('eq', "\\left|\\mathrm{FAR}_w^{\\mathrm{sketch}}-\\mathrm{FAR}_w\\right|\\le\\max_{s}\\,\\Pr_w\\left\\{s\\le S<s+\\Delta\\right\\}"),
('p', "This bound is negligible for the 512-bin sketches used here. The pooled variant reveals only aggregate counts, so no ward's individual scores leave the ward."),
]),
'E': ("Properties of the Design", [
('p', "Two properties make the design choices explicit. Write the pilot-band noise (5) as $N_k(\\rho)=c_0+c_1\\rho+c_2\\rho^{2}$, where $c_0$ collects shot noise and $c_1$, $c_2$ the two thermal terms. The CRB (6) is then proportional to $c_0/\\rho+c_1+c_2\\rho$. The linear term does not depend on $\\rho$ in this form, and setting the derivative to zero gives"),
('eq', "\\rho^{\\star}=\\sqrt{\\frac{c_0}{c_2}},\\qquad \\min_{\\rho}\\mathrm{CRB}_k=\\frac{2\\sqrt{c_0c_2}+c_1}{R^{2}\\sigma(d)^{2}T_f}"),
('p', "**Proposition 1.** The CRB-optimal pilot fraction grows with the square root of the shot-noise level and does not depend on the dimming level except through the DC illumination term of $c_0$. With the receiver of Table IV, $\\rho^{\\star}\\approx0.08$ at night and $\\approx0.22$ in daylight. Dimming changes the minimum CRB through $\\sigma(d)^2$ but hardly moves its location. The pilot share can therefore be set from the ambient light alone, and the dimming controller can act without retuning it."),
('p', "The second property concerns calibration. Pooling quantile sketches across $W$ wards guarantees the target only marginally, that is, on average over wards. For ward $w$, whose non-event score distribution $P_w$ differs from the pooled mixture $\\bar{P}$, the false-alarm rate obeys"),
('eq', "\\left|\\mathrm{FAR}_w-\\alpha\\right|\\le d_{\\mathrm{TV}}\\left(P_w,\\bar{P}\\right)+\\frac{1}{n_{\\mathrm{pool}}+1}"),
('p', "whereas local calibration has no heterogeneity term but a sampling term of order $1/\\sqrt{n_w}$, because only the $n_w$ windows of the ward set the quantile. **Proposition 2.** Pooling is preferable whenever the ward-to-mixture total-variation distance is smaller than the sampling error of the local quantile. This is the regime of a hospital with many wards and short calibration periods. Because the false-alarm rate depends only on the tail event above the threshold, the distance in the bound can be replaced by the tail discrepancy at the pooled threshold, and Appendix B checks this sharper condition ward by ward."),
]),
}
