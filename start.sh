#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Checking for large data assets..."
mkdir -p files_for_plots

# Download the massive database if it doesn't exist yet
if [ ! -f files_for_plots/TCGA_GTEx_plotting_data.bak ]; then
    echo "Downloading TCGA_GTEx_plotting_data..."
    curl -L -o files_for_plots/TCGA_GTEx_plotting_data.bak "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/TCGA_GTEx_plotting_data.bak"
    curl -L -o files_for_plots/TCGA_GTEx_plotting_data.dat "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/TCGA_GTEx_plotting_data.dat"
    curl -L -o files_for_plots/TCGA_GTEx_plotting_data.dir "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/TCGA_GTEx_plotting_data.dir"
fi

# Repeat for other essential files if needed
if [ ! -f files_for_plots/deeploc2_output.bak ]; then
    echo "Downloading deeploc2_output..."
    curl -L -o files_for_plots/deeploc2_output.bak "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/deeploc2_output.bak"
    curl -L -o files_for_plots/deeploc2_output.dat "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/deeploc2_output.dat"
    curl -L -o files_for_plots/deeploc2_output.dir "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/deeploc2_output.dir"
fi

# Repeat for other essential files if needed
if [ ! -f files_for_plots/genes.bak ]; then
    echo "Downloading genes_db..."
    curl -L -o files_for_plots/genes.bak "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/genes.bak"
    curl -L -o files_for_plots/genes.dat "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/genes.dat"
    curl -L -o files_for_plots/genes.dir "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/genes.dir"
fi

# Repeat for other essential files if needed
if [ ! -f files_for_plots/membrane_topology_objects.bak ]; then
    echo "Downloading membrane topology objects..."
    curl -L -o files_for_plots/membrane_topology_objects.bak "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/membrane_topology_objects.bak"
    curl -L -o files_for_plots/membrane_topology_objects.dat "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/membrane_topology_objects.dat"
    curl -L -o files_for_plots/membrane_topology_objects.dir "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/membrane_topology_objects.dir"
fi

# Repeat for other essential files if needed
if [ ! -f files_for_plots/transcripts_to_isoforms_mapping.bak ]; then
    echo "Downloading mapping_df..."
    curl -L -o files_for_plots/transcripts_to_isoforms_mapping.bak "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/transcripts_to_isoforms_mapping.bak"
    curl -L -o files_for_plots/transcripts_to_isoforms_mapping.dat "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/transcripts_to_isoforms_mapping.dat"
    curl -L -o files_for_plots/transcripts_to_isoforms_mapping.dir "https://github.com/apollis44/Transcript_explorer_app/releases/download/Deploy/transcripts_to_isoforms_mapping.dir"
fi

echo "Data assets ready. Starting the application..."
# Replace this with your actual app startup command (e.g., uvicorn, gunicorn, streamlit)
gunicorn app:server