awk '(NR == 1) || (FNR > 1)' ../Results/res_mod2_j*.csv > ../Results/res_final_mod2.csv

echo "stop"
