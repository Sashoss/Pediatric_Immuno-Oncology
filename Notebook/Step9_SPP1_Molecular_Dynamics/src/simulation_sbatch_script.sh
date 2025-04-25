#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Author: Ambuj Kumar, Nationwide Childrens Hospital
#
# This script accepts all input parameters, validates them (checking that
# required files/directories exist), and then creates a unique SBATCH submission 
# script based on these inputs.
#
# It also creates a log file (with a date/time–stamped filename) in both the 
# "log" and "out" directories. If any validation fails, the error is both logged
# and printed.
# -----------------------------------------------------------------------------

###############################################################################
# Usage/Help Function
###############################################################################
usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Required arguments:
    --init_gro <PATH>       Initial GROMACS coordinate file (*.gro)
    --topol <PATH>          GROMACS topology file (*.top)
    --run_prefix <STRING>   Descriptive prefix for your run (e.g., "step4.0_minimization")
    --mdp <PATH>            GROMACS parameter file (*.mdp)
    --output_dir <PATH>     Directory where simulation outputs will be stored
    --sample <STRING>       Sample identifier (e.g., 'testSample')
    --workdir <PATH>        Working directory for script operations
    --lab <STRING>          Lab/account name for #SBATCH --account

Optional arguments:
    --job_name <STRING>     Job name (default: "defaultJob")
    --run <STRING/NUMBER>   Run ID (default: "001")
    --index <PATH>          GROMACS index file (*.ndx) (optional)
    --index_select <STRING> String for interactive index selection (optional)
    --posres                posres flag (ex. -DPOSRES_ZN)
    -h, --help              Display this help message and exit

Example:
    ./$(basename "$0") \\
        --init_gro system.gro \\
        --topol topol.top \\
        --run_prefix step4.0_minimization \\
        --mdp minim.mdp \\
        --output_dir /path/to/results \\
        --sample testSample \\
        --workdir /path/to/workdir \\
        --lab myLab \\
        --job_name myJob \\
        --run_duration "10:00:00" \\
        --run 001 \\
        --index index.ndx \\
        --index_select "14" \\
        --posres "-DPOSRES_ZN"
EOF
    exit 1
}

###############################################################################
# Initialize variables
###############################################################################
INIT_GRO=""
TOPOL=""
RUN_PREFIX=""
MDP=""
OUTPUT_DIR=""
SAMPLE=""
WORKDIR=""
LAB=""
JOB_NAME="defaultJob"
RUN="001"
INDEX=""
INDEX_SELECT=""
POSRES=""
NVT_CPT=""
NPT_CPT=""
MD_CPT=""
CPI=""
EMAIL=""
RUN_DURATION="10:00:00"


###############################################################################
# Parse arguments
###############################################################################
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        --init_gro)
            INIT_GRO="$2"
            shift 2
            ;;
        --topol)
            TOPOL="$2"
            shift 2
            ;;
        --run_prefix)
            RUN_PREFIX="$2"
            shift 2
            ;;
        --mdp)
            MDP="$2"
            shift 2
            ;;
        --nvt_cpt)
            NVT_CPT="$2"
            shift 2
            ;;
        --npt_cpt)
            NPT_CPT="$2"
            shift 2
            ;;
        --md_cpt)
            MD_CPT="$2"
            shift 2
            ;;
        --cpi)
            CPI="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --workdir)
            WORKDIR="$2"
            shift 2
            ;;
        --lab)
            LAB="$2"
            shift 2
            ;;
        --job_name)
            JOB_NAME="$2"
            shift 2
            ;;
        --run_duration)
            RUN_DURATION="$2"
            shift 2
            ;;
        --run)
            RUN="$2"
            shift 2
            ;;
        --index)
            INDEX="$2"
            shift 2
            ;;
        --index_select)
            INDEX_SELECT="$2"
            shift 2
            ;;
        --email)
            EMAIL="$2"
            shift 2
            ;;
        --posres)
            POSRES="-DPOSRES"
            shift 
            ;;
        *)
            usage
            ;;
    esac
done


###############################################################################
# Validate Required Arguments and Existence of Files/Directories
###############################################################################
if [[ -z "$INIT_GRO" ]]; then
    echo "ERROR: --init_gro is required but missing."
    usage
fi
if [[ ! -f "$INIT_GRO" ]]; then
    echo "ERROR: File not found: $INIT_GRO"
    exit 1
fi

if [[ -z "$TOPOL" ]]; then
    echo "ERROR: --topol is required but missing."
    usage
fi
if [[ ! -f "$TOPOL" ]]; then
    echo "ERROR: File not found: $TOPOL"
    exit 1
fi

if [[ -z "$RUN_PREFIX" ]]; then
    echo "ERROR: RUN_PREFIX not found: $TOPOL"
    usage
fi

if [[ -z "$MDP" ]]; then
    echo "ERROR: --mdp is required but missing."
    usage
    exit 1
fi
if [[ ! -f "$MDP" ]]; then
    echo "ERROR: File not found: $MDP"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --output_dir is required but missing."
    usage
    exit 1
fi

if [[ -z "$SAMPLE" ]]; then
    echo "Sample File not found: $SAMPLE"
    usage
    exit 1
fi

if [[ -z "$WORKDIR" ]]; then
    echo "ERROR: --workdir is required but missing."
    usage
    exit 1
fi

if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: Cannot locate: $WORKDIR"
    exit 1
fi

if [[ -z "$LAB" ]]; then
    echo "ERROR: --lab is required but missing."
    usage
    exit 1
fi
if [[ -z "$EMAIL" ]]; then
    echo "ERROR: --email is required but missing."
    usage
fi


###############################################################################
# Log Directories and Files 
###############################################################################

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "The directory specified by --output_dir does not exist: $OUTPUT_DIR. Creating path."
    mkdir -p ${OUTPUT_DIR}
fi

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p ${OUTPUT_DIR}/Run_${RUN}

###############################################################################
# Adding Timestamp to sbatch script filename 
###############################################################################

SBATCH_SCRIPT="${SAMPLE}_${RUN_PREFIX}_Run_${RUN}.sbatch"
SBATCH_PATH="${WORKDIR}/${SBATCH_SCRIPT}"

###############################################################################
# SBATCH Script 
###############################################################################
cat > "$SBATCH_PATH" <<EOF
#!/usr/bin/env bash
#SBATCH --time=${RUN_DURATION}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --gres=gpu
#SBATCH --cpus-per-task=4
#SBATCH --account=${LAB}
#SBATCH --job-name=${SAMPLE}_${JOB_NAME}_${RUN}
#SBATCH --output=${WORKDIR}/slurm_logs/${SAMPLE}_${JOB_NAME}.out  

###############################################################################
# Auto-generated SBATCH Script on $(date)
###############################################################################

# Load required modules
module purge
module load GCC/10.3.0
module load OpenMPI/4.1.1
module load GROMACS/2021.3-CUDA-11.3.1


# Create per-sample output directory if it does not exist
mkdir -p "${OUTPUT_DIR}/Run_${RUN}"

# Change to the working directory
cd "${OUTPUT_DIR}/Run_${RUN}" || {
    echo "ERROR: Cannot cd to ${OUTPUT_DIR}/Run_${RUN}"
    exit 1
}

########################################
# Construct the grompp command
########################################

GROMPP_CMD=("gmx" "grompp" \
    -f "${MDP}" \
    -o "${OUTPUT_DIR}/Run_${RUN}/${RUN_PREFIX}_${RUN}.tpr" \
    -c "${INIT_GRO}" \
    -r "${INIT_GRO}"
    )

# Check for available checkpoint files in priority order: MD -> NPT -> NVT
for CPT in "${MD_CPT}" "${NPT_CPT}" "${NVT_CPT}"; do
    if [[ -n "$CPT" ]]; then
        GROMPP_CMD+=(-t "$CPT")
        break
    fi
done

# Add topology file at the end
GROMPP_CMD+=(-p "${TOPOL}")


# Add the index file only if provided
if [[ -n "${INDEX}" ]]; then
    GROMPP_CMD+=("-n" "${INDEX}")
fi

# Add the POSRES option only if provided
if [[ -n "${POSRES}" ]]; then
    GROMPP_CMD+=("${POSRES}")
fi

# Add any remaining constant options
GROMPP_CMD+=(-maxwarn -1)

########################################
# Run grompp, optionally piping index selection
########################################
if [[ -n "${INDEX_SELECT}" ]]; then
    echo "${INDEX_SELECT}" | "\${GROMPP_CMD[@]}"
else
    "\${GROMPP_CMD[@]}"
fi



########################################
# Run the actual MD
########################################

MDRUN_CMD=("gmx" "mdrun" \
    -v \\
    -deffnm "${OUTPUT_DIR}/Run_${RUN}/${RUN_PREFIX}_${RUN}" \\
    -nb gpu 
)

if [[ -n "${CPI}" ]]; then
    MDRUN_CMD+=("-cpi" "${CPI}")
fi

"\${MDRUN_CMD[@]}"

EOF

chmod +x "$SBATCH_PATH"

###############################################################################
# Adding Log for the creation of the SBATCH script
###############################################################################

echo "SBATCH script created at: ${SBATCH_PATH}"
