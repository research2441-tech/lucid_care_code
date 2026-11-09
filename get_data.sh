#!/usr/bin/env bash
# Downloads the BIDMC PPG and Respiration Dataset (PhysioNet, v1.0.0, ODC-By licence) as one MATLAB file.
# If PhysioNet's file layout changes, download bidmc_data.mat manually from https://physionet.org/content/bidmc/1.0.0/
set -e
URL="https://physionet.org/files/bidmc/1.0.0/bidmc_data.mat"
[ -f bidmc_data.mat ] || curl -L --fail -o bidmc_data.mat "$URL"
python3 -c "import scipy.io as s; d=s.loadmat('bidmc_data.mat'); print('bidmc_data.mat OK, records:', d['data'].shape[1])"
