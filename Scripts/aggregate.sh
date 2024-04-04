awk '(NR == 1) || (FNR > 1)' ../Results/unmatched/res_exp1_j*.csv > ../Results/unmatched/res_aggregate_exp1.csv
awk '(NR == 1) || (FNR > 1)' ../Results/matched/res_exp1_j*.csv > ../Results/matched/res_aggregate_exp1.csv

echo "stop"
