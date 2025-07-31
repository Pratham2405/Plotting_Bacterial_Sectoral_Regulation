# Bacterial Growth Analysis: Sectoral Models Comparison

A comprehensive study analyzing the sensitivity of choice of sectoral models on bacterial growth dynamics, comparing PTR, PTRWF, and min() models to understand their limitations and advantages in different environmental conditions.

## Project Summary & Key Objectives
This project investigates bacterial growth laws through computational modeling using three sectoral models of increasing complexity. The study focuses on how bacteria allocate resources between different protein sectors (ribosomes, transporters, metabolic proteins) and how regulatory mechanisms affect growth in both stable and fluctuating environments.

### Key Objectives
- Validate and compare sectoral models (PTR, PTRWF, min()) by gradually increasing complexity.
- Analyze the cost-benefit tradeoffs of regulatory mechanisms in bacterial growth.
- Compare transient responses and growth dynamics across different models.
- Investigate bacterial behavior in fluctuating vs. stable environmental conditions.
- Propose alternative regulation paradigms within the PTRWF framework.

## Environment Setup

### Prerequisites
- Python 3.12 or higher
- Git
  
### Installation

Clone the repository:
```
git clone <repository-url>
cd Plotting_Bacterial_Sectoral_Regulation
```
Create and activate a virtual environment (recommended):
```
python -m venv bacterial-env
source bacterial-env/bin/activate  # On Windows: bacterial-env\Scripts\activate
```
Install dependencies:
Using pip:
```
pip install -r requirements.txt
```
Or using conda:
```
conda env create -f environment.yml
conda activate bacterial-growth-env
```

## Directory Structure
```
Plotting_Bacterial_Sectoral_Regulation/
├── README.md                           # This file - project overview and instructions
├── LICENSE                             # License for the project
├── Project_Report.pdf                  # Complete documentation of the Project with Introduction, Literature Survey and Results
├── requirements.txt                    # Python dependencies for pip
├── environment.yml                     # Conda environment specification
│
├── src/                               # Source code directory
│   ├── PTR_2016/                       # Scripts for replicating PTR model(Jain, 2016)
│   └── PTRWF_2018/                     # Scripts for replicating PTRWF model(Jain, 2018)
│
├── plots/                           # Generated plots and output files
│   ├── PTRWF_2018/                  # PTRWF replication
│   ├── PTR_2016/                    # PTR replication        
│   └── min()/                       # Incorporating min() methodology into PTRWF.
|
└── docs/                             # Additional documentation
    ├── model_equations.md            # Mathematical formulations
    └── parameter_reference.md        # Parameter values and sources
```
