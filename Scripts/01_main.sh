#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=24:mem=50gb
#PBS -J 1-439

module load anaconda3/personal
source activate phd_r

cd /rds/general/user/fg520/home/Matching/Matching-Applied/Scripts/

Rscript 01_main_fake.R ${PBS_ARRAY_INDEX}

# bash aggregate.sh

echo "stop"
