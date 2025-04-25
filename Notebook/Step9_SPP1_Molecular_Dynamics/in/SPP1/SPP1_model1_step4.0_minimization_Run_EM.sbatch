#!/usr/bin/env bash
#SBATCH --time=10:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --gres=gpu
#SBATCH --cpus-per-task=16
#SBATCH --account=gdwanglab
#SBATCH --job-name=SPP1_model1_em_EM
#SBATCH --output=/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/slurm_logs/SPP1_model1_em.out  

###############################################################################
# Auto-generated SBATCH Script on Sat Apr 12 15:55:47 EDT 2025
###############################################################################

# Load required modules
module purge
module load GCC/10.3.0
module load OpenMPI/4.1.1
module load GROMACS/2021.3-CUDA-11.3.1


# Create per-sample output directory if it does not exist
mkdir -p "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/../../out/step4.0_minimization/Run_EM"

# Change to the working directory
cd "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/../../out/step4.0_minimization/Run_EM" || {
    echo "ERROR: Cannot cd to /gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/../../out/step4.0_minimization/Run_EM"
    exit 1
}

########################################
# Construct the grompp command
########################################

GROMPP_CMD=("gmx" "grompp"     -f "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/step4.0_minimization.mdp"     -o "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/../../out/step4.0_minimization/Run_EM/step4.0_minimization_EM.tpr"     -c "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/step3_input.gro"     -r "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/step3_input.gro"
    )

# Check for available checkpoint files in priority order: MD -> NPT -> NVT
for CPT in "" "" ""; do
    if [[ -n "" ]]; then
        GROMPP_CMD+=(-t "")
        break
    fi
done

# Add topology file at the end
GROMPP_CMD+=(-p "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/topol.top")


# Add the index file only if provided
if [[ -n "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/index.ndx" ]]; then
    GROMPP_CMD+=("-n" "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/index.ndx")
fi

# Add the POSRES option only if provided
if [[ -n "" ]]; then
    GROMPP_CMD+=("")
fi

# Add any remaining constant options
GROMPP_CMD+=(-maxwarn -1)

########################################
# Run grompp, optionally piping index selection
########################################
if [[ -n "" ]]; then
    echo "" | "${GROMPP_CMD[@]}"
else
    "${GROMPP_CMD[@]}"
fi



########################################
# Run the actual MD
########################################

MDRUN_CMD=("gmx" "mdrun"     -v \
    -deffnm "/gpfs0/home1/gdwanglab/axk201/personal_projects/Pediatric_glioblastoma_singelcell/Notebook/Step22_Molecular_dynamics/in/SPP1/../../out/step4.0_minimization/Run_EM/step4.0_minimization_EM" \
    -nb gpu 
)

if [[ -n "" ]]; then
    MDRUN_CMD+=("-cpi" "")
fi

"${MDRUN_CMD[@]}"

