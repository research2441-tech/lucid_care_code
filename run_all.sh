#!/usr/bin/env bash
# Regenerates every number, table and figure of the LUCID-Care paper from the BIDMC recordings.
# Expected run time on a 2-core machine: several hours (episode generation and the operating-grid sweep dominate).
set -e
export PYTHONPATH=.
export PYTHONHASHSEED=0          # exp3.pareto_job seeds its noise from hash(); fix it for repeatable runs
export BIDMC_MAT=${BIDMC_MAT:-./bidmc_data.mat}
mkdir -p results figs out
# 1. Sensing corpus: 212 main episodes (seeds 7000-7211) and 106 stress-test episodes (seeds 9000-9105)
python3 episodes.py
LC_OUT=eps_ood LC_N=106 LC_SEED=9000 LC_DWELL=90 LC_STAFF=9 python3 episodes.py   # LC_STAFF is exclusive: up to 8 visits
python3 corpus_stats3.py
# 2. Experiments on the main corpus
python3 exp3.py pareto            # rate / CRB / accuracy over the operating grid  -> results/pareto.json
python3 exp3.py wards             # ward deployment, wearable, apnea and event scores -> results/wards.pkl
python3 exp3.py ablation          # estimator ablation  -> results/ablation3.json
python3 exp3.py design            # pairs, photodiode area, body gain -> results/design3.json
python3 snr3.py                   # link budget check -> results/snr3.json
# 3. Stress-test wards and wearable-availability variants
LC_OUT=eps_ood LC_WARDS=results/wards_ood.pkl python3 exp3.py wards
for cfg in "0.0 300 20 w0" "0.3 300 20 w30" "0.5 300 20 w50" "0.15 60 20 wfreq"; do
  set -- $cfg
  LC_ABSENT=$1 LC_UP=$2 LC_DOWN=$3 LC_WARDS=results/wards_$4.pkl python3 exp3.py wards
done
# 4. Analyses
python3 analyze3.py               # fusion, Mondrian conformal alarms, events -> results/analysis3.json
python3 analyze3b.py              # alpha sweep, calibration budget, event ROC -> results/extra3.json
python3 analyze3c.py              # stress-test and record-disjoint checks -> results/external3.json
python3 analyze3d.py              # wearable availability -> results/wear3.json
python3 events_type3.py           # per-type events, apnea onset latency -> results/events_type3.json
python3 errband3.py               # error by respiratory-rate band, Bland-Altman -> results/errband3.json
python3 appendix3.py              # ward-level check of Proposition 2 -> results/appendix3.json
# 5. Figures and manuscript (needs the IEEE template saved as ./IEEE_template.docx and pandoc on PATH)
python3 make_figs3.py fig_arch fig_signal fig_pareto fig_fusion fig_design fig_calib fig_snr fig_analytic
[ -f IEEE_template.docx ] && python3 build_paper3.py || echo "IEEE_template.docx not found: manuscript build skipped"
