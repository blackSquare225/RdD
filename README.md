# RΔD

## Description
This is a Snakemake-based pipeline designed for fast and efficient identification of deletions in target long-read sequencing data. It automates the process of QC, alignment and deletion calling using modular workflows and reproducible environments.

## Requirements

Before using the pipeline, ensure you have the following installed:
- **[Miniconda/Conda](https://docs.conda.io/en/latest/miniconda.html)** or **Mamba** (recommended for faster package resolution)
- **Git** (for cloning the repository)

## Installation

### Step 1: Clone the Repository
To get started, clone this repository to your local machine:

```bash
git clone https://github.com/blackSquare225/RdD.git
cd RdD
```

### Step 2: Install Snakemake and Set Up the mamba environment

We recommend using Mamba for installing Snakemake as it is faster than Conda. 
First, install Mamba (if not installed) working on base enviroment:
```bash
conda install mamba -c conda-forge
```

Then, create the Snakemake environment using mamba:
```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake
mamba activate snakemake_env
```

Check that Snakemake is correctly installed by running:
```bash
snakemake --version
```

### Step 3: Run the pipeline 
