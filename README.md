# RΔD

## Description
This is a Snakemake-based pipeline that integrates a custom algorithm (DELecter) designed for fast and efficient identification of deletions in target DNA long-read sequencing data. It automates the process of QC, alignment and deletion calling using modular workflows and reproducible environments.

# Author

Alessio Enderti 

- email: **alessioenderti@gmail.com**
- personal github: **https://github.com/endertial/**

For questions or support, reach out via GitHub Issues or email the author

## 📖 Table of Contents
2. [🛠 Installation](#-installation)
3. [📌 Quick Usage](#-quick-usage)
4. [🔍 How the Code Works](#-how-the-script-works)
   - [🔹 Multi-Sample Mode](#-multi-sample-mode)
   - [🔹 Single-Sample Mode](#-single-sample-mode)
   - [🚀 Downstream analysis](#-downstream-analysis)
5. [❗List of Commands](#-list-of-commands)
7. [👨‍💻 Contributing](#-contributing)
8. [📖 Citing](#-citing)


## 🛠 Installation

### Step 1: Install Snakemake and Set Up the mamba environment

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

### Step 2: Clone the Repository
To get started, clone this repository on your local machine:

```bash
git clone https://github.com/blackSquare225/RdD.git
cd RdD
```


## 📌 Quick Usage

### **1️⃣ Multi-Sample Mode**
Processes all fastq.gz files inside barcoded directories in fast_pass/ output ONT.
```bash
bash run.sh --multi -name <run_name> --fastqdir </path/to/fastq_pass/> --threads <num_of_threads>
```

### **2️⃣ Single-Sample Mode**
Processes a single FASTQ file.
```bash
./script.sh --single -name <sample_name> --fastqfile <path/to/fastq.gz/file> --threads <num_of_threads>
```


## 🔍 How the Code Works

### **🔹 Multi-Sample Mode**
ONT (Oxford Nanopore Technology) sequencing allows for multiplexing, meaning multiple samples can be sequenced in a single run. To differentiate them, barcodes are added to each sample before sequencing. The sequencing software then sorts the reads into separate barcoded directories.
When using multi-sample mode, the script processes all barcoded directories inside the fastq_pass/ (`--fastqdir`) folder, which is the standard output directory of an ONT sequencing run. The script scans the `fastq_pass/` directory for all barcode folders (barcode01, barcode02, etc.). For each barcode folder it merges all FASTQ files inside that barcode folder into a single FASTQ file. The merged file is stored in the `data/` directory with a naming convention that includes the run name (`--name`) and barcode. It then updates the configuration file (`config.yaml`) to ensure snakemake run, which processes each sample separately.

In summary, multi-sample mode is ideal for analyzing multiple samples directly from ONT output, without the need to manually reorganize or move files. This ensures an efficient and streamlined workflow for handling multiplexed sequencing runs.

### **🔹 Single-Sample Mode**
In single-sample mode, the script is designed to process a single FASTQ file instead of scanning multiple barcoded directories. This mode is useful when sequencing only one sample or when working with a pre-extracted FASTQ file instead of a full Nanopore sequencing run.
The user must directly specify the path to a single fastq.gz (`--fastqfile`). Additionally, the user provides a sample name (`--name`), which is used in the output configuration.
Since Snakemake requires input files to follow a specific pre-defined format, the script ensures compatibility by creating a symbolic link inside the `data/` directory, renaming the input file accordingly. This allows the pipeline to recognize and process the sample correctly, even if the original file has a different name.

In summary, single-sample mode is ideal for cases where only one sample needs to be analyzed, providing flexibility while maintaining compatibility with the Snakemake workflow.

### **🚀 Downstream analysis**
This pipeline automates the processing of ONT sequencing data by handling both multi-sample and single-sample modes. It aligns sequencing reads using `minimap2 -ax map-ont`, sorts the resulting BAM files, and then processes them using **DELecter**, a custom DEL structural variant (SV) detection tool. DELecter implements a union-find algorithm to cluster reads that belong to the same DEL SV event based on the proximity of their start and end breakpoints enabling precise identification with high confidence.


## ❗ List of Commands

### Options to select run mode:
| Command | Description | 
|---------|-------------|
|  -m, --multi     |  Run in multi-sample mode (requires run name and fastq directory) |
|  -s, --single    |  Run in single-sample mode (requires sample name and fastq file) |
|  -n, --name      |  Name of the run (for multi-sample mode) or sample name (for single-sample mode) |
|  -d, --fastqdir  |  Path to barcoded directories (e.g., /path/to/fastq_pass/) (only for multi-sample mode) |
|  -f, --fastqfile |  Path to the FASTQ file (only for single-sample mode) |
|  -t, --threads   |  Number of threads to use (default: 1) |
|  -h, --help      |  Display this help message |

### Options for deletion variant calling by DELecter
| Command | Description |
|---------|-------------| 
|  --min_size        |  Minimum deletion size (in bp) to take into account. (default: 300) |
|  --max_size        |  Maximum deletion size (in bp) to take into account. (default: 4000) |
|  --mapq            |   Reads with mapping quality lower than this value will be ignored. (default: 25) |
|  --min_support     |  Minimum number of supporting reads for a DEL to be reported. (default 100) |
|  --min_len         |   Minimum DEL length (in bp) to be reported. (default 1000) |
|  --exclude_flag    |   Bitwise sam flag to exclude. (default 3844) |
|  --tolerance       |   Tolerance between breakpoints (in bp) for clustering deletions as the same event. (default: 100) |


## 👨‍💻 Contributing
Feel free to open an issue or submit a pull request if you find a bug or want to improve the code!

## 📖 Citing

Are you using RdD in your works? Please cite:

