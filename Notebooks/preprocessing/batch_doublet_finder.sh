#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --partition=himem
#SBATCH --mem=256G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --account=gdwanglab
#SBATCH --job-name=hgg_seurat



ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.1


cd /home/gdwanglab/axk201/personal_projects/HGG/Notebooks/preprocessing

# Batch DoubletFinder Shell Script



# Input directory containing .rds files
INPUT_DIR=../create_seurat_objects/out/

# Output directory for processed files
OUTPUT_DIR=./out/

# Path to the R script
R_SCRIPT=./doublet_finder.r

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each .rds file in the input directory
for file in "$INPUT_DIR"/*.rds; do
  # Check if there are any .rds files
  if [ ! -e "$file" ]; then
    echo "No .rds files found in $INPUT_DIR"
    exit 1
  fi

  # Get the base filename (without the directory path)
  base_filename=$(basename "$file")

  # Log the current file being processed
  echo "Processing file: $file"

  # Run the R script
  Rscript "$R_SCRIPT" -i "$file" -o "$OUTPUT_DIR"

  # Log success
  echo "Finished processing file: $file"
done

echo "All files processed. Output files are in $OUTPUT_DIR."

