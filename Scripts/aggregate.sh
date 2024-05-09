#!/bin/bash
foldername="simul1"
#awk '(NR == 1) || (FNR > 1)' ../Results/${foldername}/res_j*_norm_summary.csv > ../Results/${foldername}/res_aggregate.csv
#awk '(NR == 1) || (FNR > 1)' ../Results/${foldername}/res_j*_pk1_summary.csv > ../Results/${foldername}/res_aggregate_pk1.csv
#awk '(NR == 1) || (FNR > 1)' ../Results/${foldername}/res_j*_null_summary.csv > ../Results/${foldername}/res_aggregate_null.csv
mlr --csv unsparsify ../Results/${foldername}/res_*summary.csv > ../Results/${foldername}/res_aggregate.csv

echo "stop"
