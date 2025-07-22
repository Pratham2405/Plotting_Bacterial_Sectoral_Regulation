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

## Directory Structure
```
bacterial-growth-models/
├── README.md                           # This file - project overview and instructions
├── LICENSE                             # License for the project
├── requirements.txt                    # Python dependencies for pip
├── environment.yml                     # Conda environment specification
│
├── src/                               # Source code directory
│   ├── pop_v_t.py                     # Population dynamics with division simulation
│   ├── mu_v_f.py                      # Growth rate vs food concentration analysis
│   ├── CNB.py                         # Cost-No-Benefit model implementation
│   ├── mu-vs-F_e-w-o-Division.py     # Growth laws without division effects
│   └── n_v_t.py                      # Fluctuating environment analysis
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
