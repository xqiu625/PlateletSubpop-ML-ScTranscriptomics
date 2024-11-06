#!/bin/bash

# Data Download Script for PlateletSubpop-ML-ScTranscriptomics
# This script downloads and processes single-cell RNA sequencing datasets

set -e  # Exit on error
set -u  # Exit on undefined variable

# Create necessary directories
mkdir -p data/{raw,processed}/{covid19,sepsis,sle}
mkdir -p data/tmp
mkdir -p logs

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a logs/download.log
}

# Function to download and process GEO datasets
download_geo_dataset() {
    local accession=$1
    local category=$2
    
    log_message "Downloading ${accession}..."
    cd data/tmp
    wget -r -np -nd -R "index.html*" "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${accession:0:8}nnn/${accession}/suppl/"
    
    if [ -f "${accession}_RAW.tar" ]; then
        tar -xvf "${accession}_RAW.tar"
    fi
    
    find . -name "*.gz" -exec gunzip {} \;
    
    # Move processed files to appropriate directory
    mv * "../raw/${category}/"
    cd ../../
    
    log_message "Completed processing ${accession}"
}

# Function to handle GSE158055 special case
process_gse158055() {
    log_message "Processing GSE158055 (large dataset)..."
    cd data/raw/covid19
    
    # Split the large matrix file
    sed -n '4,1159554303p' GSE158055_covid19_counts.mtx | \
        awk '{print $1"\t"$2"\t"$3}' > GSE158055_covid19_part1
    sed -n '1159554304,2319108602p' GSE158055_covid19_counts.mtx | \
        awk '{print $1"\t"$2"\t"$3}' > GSE158055_covid19_part2
    
    # Add headers to split files
    for part in {1,2}; do
        echo "%%MatrixMarket matrix coordinate real general" > header
        echo "%" >> header
        echo "27943 1462702 1159554300" >> header
        cat header GSE158055_covid19_part${part} > temp
        mv temp GSE158055_covid19_part${part}
    done
    
    # Create separate folders for each part
    for part in {1,2}; do
        mkdir -p GSE158055_part${part}
        mv GSE158055_covid19_part${part} GSE158055_part${part}/matrix.mtx
        cp GSE158055_covid19_features.tsv GSE158055_part${part}/features.tsv
        cp GSE158055_covid19_barcodes.tsv GSE158055_part${part}/barcodes.tsv
        gzip GSE158055_part${part}/*
    done
    
    cd ../../../
    log_message "Completed GSE158055 processing"
}

# Main execution
log_message "Starting data download process..."

# COVID-19 datasets
download_geo_dataset "GSE150728" "covid19"
download_geo_dataset "GSE155673" "covid19"
download_geo_dataset "GSE158055" "covid19"
process_gse158055
download_geo_dataset "GSE151263" "covid19"
download_geo_dataset "GSE163668" "covid19"

# Download E-MTAB-10026
log_message "Downloading E-MTAB-10026..."
cd data/raw/covid19
for i in {1..4}; do
    wget "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10026/E-MTAB-10026.processed.${i}.zip"
    unzip "E-MTAB-10026.processed.${i}.zip"
done
wget "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10026/E-MTAB-10026.sdrf.txt"
cd ../../../

# SLE dataset
download_geo_dataset "GSE142016" "sle"

# Clean up
rm -rf data/tmp

log_message "Data download complete!"
