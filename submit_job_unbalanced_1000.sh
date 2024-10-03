#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haoge.chang@yale.edu
#SBATCH --partition=week
#SBATCH --time=72:00:00
#SBATCH --array=1

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID


cd /vast/palmer/home.grace/hc654/BinaryOutcomePermutationTest



module load R/4.2.0-foss-2020b

module load miniconda
conda deactivate
conda activate r-general

lscpu
Rscript --vanilla 1_simulation_unbalanced_case_1_1000.R $SLURM_ARRAY_TASK_ID
