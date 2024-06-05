#!/bin/bash

echo "Processed files:"

ls -t | grep -v "all" | grep 'proj_pcone_n100_.*_tol_high.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n100_tol_high_all.csv
ls -t | grep -v "all" | grep 'proj_pcone_n300_.*_tol_high.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n300_tol_high_all.csv
ls -t | grep -v "all" | grep 'proj_pcone_n500_.*_tol_high.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n500_tol_high_all.csv
ls -t | grep -v "all" | grep 'proj_pcone_n1000_.*_tol_high.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n1000_tol_high_all.csv

ls -t | grep -v "all" | grep 'proj_pcone_n100_.*_tol_low.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n100_tol_low_all.csv
ls -t | grep -v "all" | grep 'proj_pcone_n300_.*_tol_low.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n300_tol_low_all.csv
ls -t | grep -v "all" | grep 'proj_pcone_n500_.*_tol_low.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n500_tol_low_all.csv
ls -t | grep -v "all" | grep 'proj_pcone_n1000_.*_tol_low.*csv'  | head -n 4 | sort | tee /dev/tty | xargs paste -d "," > proj_pcone_n1000_tol_low_all.csv

ls -t | grep "proj_hyper_n10_d1_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename -f  's/.*tol_high\K.*$/\.csv/'
ls -t | grep "proj_hyper_n10_d2_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_high\K.*$/\.csv/'

ls -t | grep "proj_hyper_n20_d1_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_high\K.*$/\.csv/'
ls -t | grep "proj_hyper_n20_d2_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename  -f 's/.*tol_high\K.*$/\.csv/'

ls -t | grep "proj_hyper_n30_d27_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_high\K.*$/\.csv/'
ls -t | grep "proj_hyper_n40_d37_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_high\K.*$/\.csv/'
ls -t | grep "proj_hyper_n50_d47_30_tol_high.*stats.csv" | head -n 1 | tee /dev/tty | rename -f  's/.*tol_high\K.*$/\.csv/'


ls -t | grep "proj_hyper_n10_d1_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename -f  's/.*tol_low\K.*$/\.csv/'
ls -t | grep "proj_hyper_n10_d2_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_low\K.*$/\.csv/'

ls -t | grep "proj_hyper_n20_d1_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_low\K.*$/\.csv/'
ls -t | grep "proj_hyper_n20_d2_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename  -f 's/.*tol_low\K.*$/\.csv/'

ls -t | grep "proj_hyper_n30_d27_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_low\K.*$/\.csv/'
ls -t | grep "proj_hyper_n40_d37_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename -f 's/.*tol_low\K.*$/\.csv/'
ls -t | grep "proj_hyper_n50_d47_30_tol_low.*stats.csv" | head -n 1 | tee /dev/tty | rename -f  's/.*tol_low\K.*$/\.csv/'
