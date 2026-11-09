# LUCID-Care: code for the paper

This code goes with "LUCID-Care: Dual-Function Visible-Light Sensing and Communication for Contactless-Plus-Wearable Patient Monitoring with Federated Conformal Alarm Calibration across Hospital Wards".

Every number, table and figure in the manuscript comes from these scripts. Nothing in the paper was typed in by hand.

## Setup

```
pip install -r requirements.txt
./get_data.sh            # BIDMC PPG and Respiration Dataset (PhysioNet 1.0.0) -> ./bidmc_data.mat
./run_all.sh             # full pipeline -> results/, figs/, out/
```

- **Data path.** Set `BIDMC_MAT=/path/to/bidmc_data.mat` if the file is kept somewhere else.
- **Manuscript build.** This step needs pandoc and the IEEE template saved as `./IEEE_template.docx`.
- **Run time.** The full run takes several hours on two cores. Episode generation and the operating-grid sweep take most of that time.

## Pipeline

| Stage | Script | Output |
|---|---|---|
| Ward geometry and ray-traced luminaire–body–photodiode gains | `rooms.py`, `raytrace.py`, `bvh.py`, `simulate.py`, `ward.py` | used by `episodes.py` |
| Torso patches and BIDMC impedance respiration | `sense.py` | — |
| 8 min sensing episodes: posture, staff, apneas, bed exits and falls | `episodes.py` | `eps/`, `eps_ood/` |
| Dual-function DCO-OFDM model, CRB, downlink rate, artifact gating, eigen-combining | `phy.py`, `isac.py` | — |
| Operating grid, ward deployment, ablation, design sweeps | `exp3.py pareto / wards / ablation / design` | `results/*.json`, `results/wards*.pkl` |
| Link-budget check | `snr3.py` | `results/snr3.json` |
| Fusion, Mondrian conformal alarms (fixed, local, pooled, federated-normalised), events | `analyze3.py` | `results/analysis3.json` |
| α sweep, calibration budget, ROC | `analyze3b.py` | `results/extra3.json` |
| Stress-test wards and record-disjoint split | `analyze3c.py` | `results/external3.json` |
| Wearable availability | `analyze3d.py` | `results/wear3.json` |
| Per-type events and apnea onset latency | `events_type3.py` | `results/events_type3.json` |
| Error by respiratory-rate band and Bland–Altman | `errband3.py` | `results/errband3.json` |
| Ward-level check of Proposition 2 | `appendix3.py` | `results/appendix3.json` |
| Corpus statistics and timing | `corpus_stats3.py` | `results/stats3.json` |
| Figures | `make_figs3.py` | `figs/` |
| Manuscript (OMML equations, numbered citations) | `paper3_text.py`, `refs.py`, `cite.py`, `docbuild.py`, `build_paper3.py` | `out/LUCID_Care_Q1_Paper.docx` |

## Reproducibility notes

- **Seeds.** Episodes are generated from fixed seeds: 7000–7211 for the main corpus and 9000–9105 for the stress-test wards. They regenerate bit-exactly.
- **Noise draws.** Most noise draws use explicit NumPy generators.
  - The one exception is the operating-grid sweep (`exp3.pareto_job`), which seeds its noise from Python's `hash()`.
  - `run_all.sh` fixes `PYTHONHASHSEED=0` so that this sweep repeats exactly.
  - The published grid was produced without that setting, so a rerun reproduces it only up to the noise of the 53 episodes at each grid point.
- **Timing.** The timing values in `results/stats3.json` depend on the machine.
- **No other datasets.** BIDMC is the only external dataset that is required.
