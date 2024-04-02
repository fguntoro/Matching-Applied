#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=24:mem=50gb
#PBS -J 1-18

module load anaconda3/personal
source activate phd_r

cd /rds/general/user/fg520/home/Matching/Matching-Applied/Scripts/

Rscript main.R ${PBS_ARRAY_INDEX}

echo "stop"
