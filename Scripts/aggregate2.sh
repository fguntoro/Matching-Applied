awk '(NR == 1) || (FNR > 1)' ../Results/res_exp2_j*_pk1_summary.csv > ../Results/res_aggregate_exp2_pk1.csv
awk '(NR == 1) || (FNR > 1)' ../Results/res_exp2_j*_null_summary.csv > ../Results/res_aggregate_exp2_null.csv

echo "stop"
