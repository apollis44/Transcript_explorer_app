# Transcript Explorer - Dash Application

This directory contains the interactive frontend for the Transcript Explorer application. It is a [Plotly Dash](https://dash.plotly.com/) web application built with Bootstrap components, designed to dynamically visualize transcript-level data.

## Features

- **Protein Search:** Fast autocompletion search from a precomputed database of human proteins from the Ensembl database.
- **Subcellular Localization:** A custom-styled heatmap displaying the predicted probability scores of different cellular destinations for each transcript isoform. It features a top-aligned, slender colorbar ("Probability") and a soft, high-contrast gradient (royal blue to crimson).
- **Membrane Topology:** Visualizes the predicted structural features (signal peptides, transmembrane, extracellular, and intracellular regions) along the amino acid alignments for each transcript isoform. The legend is aligned vertically on the top right.
- **Expression Analysis:** Visualizes expression distributions across various healthy tissues and cancer types (TCGA and GTEx datasets).

## Codebase Structure

The application consists of the following components:

- **[app.py](file:///c:/Users/antoine.andreoletti/Documents/Transcript_explorer_app/app.py):** The main application file. Defines the Dash layout, sidebar navigation, routing logic, and state management callbacks.
- **[scripts/plots_generation.py](file:///c:/Users/antoine.andreoletti/Documents/Transcript_explorer_app/scripts/plots_generation.py):** Contains modular functions (`create_localization_plot`, `create_topology_plot`, `plot_expression_data`) to generate the interactive Plotly figures.
- **[assets/custom.css](file:///c:/Users/antoine.andreoletti/Documents/Transcript_explorer_app/assets/custom.css):** The styling definition for the application, establishing the premium dark theme, autocompletes, card panels, and custom-colored form switches.
- **`files_for_plots/`:** Database directory containing shelve archives (`genes`, `deeploc2_output`, `TCGA_GTEx_plotting_data`, `transcripts_to_isoforms_mapping`, and `membrane_topology_objects`) containing precomputed datasets.

## Prerequisites

Before running the application, ensure that the databases are available in the release assets of the repository. If they are not available, you can generate them by running the [Transcripts_explorer](https://github.com/apollis44/Transcripts_explorer) pipeline.

## Running the App

To run the application locally:

```bash
./start.sh
```