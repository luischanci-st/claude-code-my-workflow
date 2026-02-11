#################################################################
# - Research Project:
#   Interconnectedness, Performance, Risk, and Response 
#   to Shocks: Evidence from the Banking Network in Chile.
# - Article:
#   Chanci, Kumbhakar and Bobadilla - 
#   Interconnectedness and Banking Performance
#   2025.
#################################################################

=================================================================

OVERVIEW

This replication package contains the complete codebase and instructions
to reproduce the results of the research paper "Interconnectedness
and Banking Performance".

The project includes data processing scripts (R) and estimation 
modules (Julia), featuring a CPU-parallelized version for efficiency.

=================================================================

FOLDER STRUCTURE

The project assumes a single root folder named "Chkb". The internal
structure is organized as follows:

Chkb/
├── main.jl               # MAIN SCRIPT - Run this to execute the full estimation pipeline.
├── Chkb.jl               # Julia Module - Contains all estimation functions/algorithms.
├── dataprocessing.R      # R Script - Processes raw CSV data into .Rda format for Julia.
├── readme.txt            # This file.
├── Data/
│   ├── Raw/              # Place raw CSV/Excel files here (See Data Requirements).
│   └── Processed/        # Output folder for processed .Rda files (Automatically generated).
└── outcomes/
│   ├── includes/         # Estimation results (txt files) will be saved here.

=================================================================

DATA REQUIREMENTS

Before running the code, ensure the following raw files are placed
in the "Chkb/Data/Raw/" folder:

C18_Muestra.csv        # Interbank market transactions (CMF)

Empl_muestra.csv       # Employment data (CMF)

MB2_Muestra.csv        # Balance sheet data (CMF)

APR_Muestra.csv        # Risk-weighted assets and solvency indicators (CMF)

IPC_06072024162200.xlsx # CPI data (Central Bank of Chile)

Note: The "Muestra" files are samples. For full replication, the CMF
must replace these with the complete regulatory datasets.

=================================================================

SOFTWARE/HARDWARE (specifications and requirements)

The following two software are required:

- Julia (v1.9) ; Check https://julialang.org/downloads/ 
- R (v4.5.1 was used)

Required Julia Packages (Automatically installed by main.jl):
    Distributed, DataFrames, LinearAlgebra, Distributions, Optim
    Statistics, StatsBase, GLM, Random, RData, CategoricalArrays
    Dates, StatsModels, RCall, Printf

Required R Packages (Automatically checked by dataprocessing.R):
    tidyverse, readxl, igraph

All the files have been tested on: 

- A Workstation with an Intel Core i9-14900KF Processor, 
  128GB of DDR5 RAM (5600MT/s), and 
  an NVIDIA GeForce RTX 5090 (32GB GDDR7).

=================================================================

REPLICATION INSTRUCTIONS

Step 1: Setup

    Unzip the "Chkb" folder to any location on your computer.

    Ensure the "Data/Raw" folder contains the required CSV/Excel files.

Step 2: Execution

    Open Julia. It could also be from VScode (as I did)

    Open and run the file "main.jl".

Details:

- The script automatically sets the working directory to the "Chkb" folder.
- It checks for and installs missing Julia packages.
- It calls R (via RCall) to execute "dataprocessing.R", generating
the necessary dataframes in "Data/Processed/".
- It loads the "Chkb" module.
- It runs the estimation for Models M1, M2, and M3 using parallel computing (CPU).
- Results are printed to the console and saved in "outcomes/includes/".

=================================================================

OUTPUTS

The code generates the following results in "outcomes/includes/":

- Point estimates for parameters (rho, sigma_v, sigma_u, deltas).
- Standard errors (*_se.txt).
- Significance asterisks (*_pval.txt).

These text files are formatted for direct inclusion in LaTeX tables.

=================================================================

#################################################################
CONTACT

Luis Chanci
luischanci at santotomas dot cl
www.luischanci.com
#################################################################
