
MMPBSA_RUN_PATH="Path where MMPBSA should be executed"
cd ${MMPBSA_RUN_PATH}

set +e
unset GROUPS

ml purge
ml gmx_MMPBSA/1.6.3

source /export/apps/opt/envs/gmx_MMPBSA-1.6.3/amber.sh

MD_DIR="path where MD simulation results are"  
MMPBSA_DIR="gmx_MMPBSA-master/examples/Protein_protein"
MMPBSA_INPUT="mmpbsa.in"
MMPBSA_OUTPUT="FINAL_RESULTS_MMPBSA"
ANALYSIS_DIR_SUFFIX="_Analysis"

XTC_START=0  # Start time
XTC_END=50000    # End time

# Binding groups
declare -a BINDING_GROUPS=(
    "16 17"  # SPP1-Heavy chain
    "16 18"  # SPP1-Light chain
)

SAMPLES=(
        "23C3_v1 ${MD_DIR}/in/23C3_v1/index_prot_chains.ndx" \
        "Hu23C3_v1 ${MD_DIR}/in/Hu23C3_v1/index_prot_chains.ndx"
)

RUNS=('Run_MD_chain1' 'Run_MD_chain2' 'Run_MD_chain3')
for SAMPLE_ENTRY in "${SAMPLES[@]}"
do
    SAMPLE_NAME=$(echo $SAMPLE_ENTRY | awk '{print $1}')
    SAMPLE_DIR="${MD_DIR}/out/${SAMPLE_NAME}"
    INDEX_FILE_PATH=$(echo $SAMPLE_ENTRY | awk '{print $2}')

    echo "=============================="
    echo "Processing sample: ${SAMPLE_NAME}"
    echo "=============================="

    if [ ! -f "${INDEX_FILE_PATH}" ]; then
        echo "Error: Index file ${INDEX_FILE_PATH} does not exist."
        echo "Please create the index file using gmx make_ndx."
        exit 1
    fi

    for RUN in "${RUNS[@]}"
    do
        echo "----------------------------------------"
        echo "Processing run: ${RUN} for sample: ${SAMPLE_NAME}"
        echo "----------------------------------------"

        ANALYSIS_DIR="${MMPBSA_RUN_PATH}/out/${SAMPLE_NAME}/${RUN}"
        if [ ! -d "${ANALYSIS_DIR}" ]; then
            mkdir -p "${ANALYSIS_DIR}"
        fi

        run_count=${RUN: -1}

        INPUT_TPR="${MD_DIR}/out/${SAMPLE_NAME}/md/${RUN}/md_MD_chain${run_count}.tpr"
        INPUT_XTC="${MD_DIR}/out/${SAMPLE_NAME}/md/${RUN}/centered.xtc"
        TRIMMED_XTC="${MD_DIR}/out/${SAMPLE_NAME}/md/${RUN}/centered_trimmed.xtc"
        if [ ! -f "${INPUT_TPR}" ]; then
            echo "Warning: Input TPR file ${INPUT_TPR} not found. Skipping run."
            continue
        fi

        if [ ! -f "${INPUT_XTC}" ]; then
            echo "Warning: Input XTC file ${INPUT_XTC} not found. Skipping run."
            continue
        fi

        echo "Trimming trajectory for run: ${RUN}"
        echo 0 | gmx trjconv -s "${INPUT_TPR}" -f "${INPUT_XTC}" -n "${INDEX_FILE_PATH}" -o "${TRIMMED_XTC}" -b "${XTC_START}" -e "${XTC_END}"

        GROUPS=("${BINDING_GROUPS[@]}")
        echo "Assigned GROUPS: ${GROUPS[@]}"
        echo "Group comparisons for sample ${SAMPLE_NAME}: ${GROUPS[@]}"

        if [ "${#GROUPS[@]}" -gt 0 ]; then
            for CG_PAIR in "${GROUPS[@]}"
            do
                read -r CG1 CG2 <<< "${CG_PAIR}"  
                OUTPUT_DAT="${ANALYSIS_DIR}/${MMPBSA_OUTPUT}_${CG1}_${CG2}.dat"
                OUTPUT_CSV="${ANALYSIS_DIR}/${MMPBSA_OUTPUT}_${CG1}_${CG2}.csv"

                echo "Running gmx_MMPBSA for group pair: ${CG_PAIR}"
                gmx_MMPBSA -O \
                      -i "${MMPBSA_DIR}/${MMPBSA_INPUT}" \
                      -cs "${INPUT_TPR}" \
                      -ct "${TRIMMED_XTC}" \
                      -ci "${INDEX_FILE_PATH}" \
                      -cg "${CG1}" "${CG2}" \
                      -cp "${MD_DIR}/in/${SAMPLE_NAME}/topol.top" \
                      -o "${OUTPUT_DAT}" \
                      -eo "${OUTPUT_CSV}"
            done
            echo "MMPBSA calculation for run: ${RUN} completed."
        else
            echo "WARNING: GROUPS is empty. Skipping GROUPS loop for sample ${SAMPLE_NAME}, run ${RUN}"
        fi
    done

    echo "All runs for sample: ${SAMPLE_NAME} completed."
done

echo "=============================="
echo "All samples processed successfully."
echo "=============================="

# Clean dir
rm _GMXMMPBSA_COM*
rm COM.prmtop
rm *.log 