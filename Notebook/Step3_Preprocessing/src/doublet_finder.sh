#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --partition=himem
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=gdwanglab
#SBATCH --job-name=doublets



ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.1


cd /home/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook1/Step3_Preprocessing

Rscript src/doublet_finder.r