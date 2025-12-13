#!/bin/bash
#SBATCH --job-name=platelet_ml
#SBATCH --output=platelet_ml_%j.log
#SBATCH --error=platelet_ml_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=godziklab

# ==============================================================================
# Platelet ML Pipeline - HPCC Submission Script
# ==============================================================================
# Usage:
#   sbatch scripts/run_ml_pipeline_hpcc.sh
#
# Or run interactively:
#   bash scripts/run_ml_pipeline_hpcc.sh
# ==============================================================================

# Set paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302005"
SCRIPT_DIR="${BASE_DIR}/scripts"
DATA_FILE="${BASE_DIR}/302005_platelet_harmony_integrated.h5ad"

# Load modules (adjust based on your HPCC)
module load R/4.2.0
module load python/3.9  # For h5ad conversion if needed

echo "=========================================="
echo "Platelet ML Pipeline"
echo "=========================================="
echo "Start time: $(date)"
echo "Working directory: ${BASE_DIR}"
echo "Data file: ${DATA_FILE}"
echo ""

# Change to project directory
cd ${BASE_DIR}

# Step 1: Convert h5ad to Seurat (if needed)
RDS_FILE="${BASE_DIR}/302005_platelet_harmony_integrated.rds"
if [ ! -f "${RDS_FILE}" ]; then
    echo "Step 1: Converting h5ad to Seurat format..."
    Rscript ${SCRIPT_DIR}/data/load_h5ad_data.R
else
    echo "Step 1: Seurat file already exists, skipping conversion..."
fi

# Step 2: Run ML comparison pipeline
echo ""
echo "Step 2: Running ML comparison pipeline..."
Rscript -e "
library(Seurat)
source('${SCRIPT_DIR}/models/ML_comparison_pipeline.R')

# Load data
message('Loading Seurat object...')
seurat_obj <- readRDS('${RDS_FILE}')
DefaultAssay(seurat_obj) <- 'RNA'

# Check available columns for classification
message('Available metadata columns:')
print(colnames(seurat_obj@meta.data))

# Run pipeline (adjust outcome_col based on your data)
# You may need to modify the outcome column name
results <- run_ml_pipeline(
  seurat_obj,
  output_dir = '${BASE_DIR}/results/ml_comparison'
)

message('Pipeline completed!')
"

# Step 3: Generate visualizations
echo ""
echo "Step 3: Generating visualizations..."
Rscript -e "
source('${SCRIPT_DIR}/visualization/ML_visualization.R')

# Load results
ml_results <- readRDS('${BASE_DIR}/results/ml_comparison/ml_pipeline_results.rds')

# Generate all plots
generate_all_plots(
  ml_results,
  output_dir = '${BASE_DIR}/results/figures',
  format = 'png'
)

message('Visualizations generated!')
"

echo ""
echo "=========================================="
echo "Pipeline completed!"
echo "End time: $(date)"
echo ""
echo "Results saved to:"
echo "  - ${BASE_DIR}/results/ml_comparison/"
echo "  - ${BASE_DIR}/results/figures/"
echo "=========================================="
