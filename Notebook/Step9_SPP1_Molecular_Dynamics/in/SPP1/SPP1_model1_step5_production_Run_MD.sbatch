#!/usr/bin/env bash
#SBATCH --time=200:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --gres=gpu:2
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --account=gdwanglab
#SBATCH --job-name=replex_md
#SBATCH --output=prod.out  

module purge
module load GCC/10.3.0
module load OpenMPI/4.1.1
module load GROMACS/2021.3-CUDA-11.3.1

cd "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/../../out/step5_production" 

mpirun -np 4 gmx_mpi mdrun \
            -multidir   step5_production_283 \
                        step5_production_303 \
                        step5_production_333 \
                        step5_production_353 \
            -replex 500000 \
            -deffnm remd \
            -nb gpu 