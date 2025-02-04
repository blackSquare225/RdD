#!/bin/bash

# Function to print usage information
usage() {
    echo "Usage: bash $0 [-m -n <run_name> -f <fastq_directory>] [-s -n <sample_name> -f <fastq_file>] [-t <threads>]"
    echo ""
    echo "Options:"
    echo "  -m, --multi      Run in multi-sample mode (requires run name and fastq directory)"
    echo "  -s, --single     Run in single-sample mode (requires sample name and fastq file)"
    echo "  -n, --name       Name of the run (for multi-sample mode) or sample name (for single-sample mode)"
    echo "  -f, --fastqdir   Path to barcoded directories (e.g., /path/to/fastq_pass/) (only for multi-sample mode)"
    echo "  -f, --fastqfile  Path to the FASTQ file (only for single-sample mode)"
    echo "  -t, --threads    Number of threads to use (default: 1)"
    echo "  -h, --help       Display this help message"
    echo ""
    echo "Examples:"
    echo "  $0 -m -n run01 -f /path/to/fastq_pass/ -t 8"
    echo "  $0 -s -n sample01 -f /path/to/sample.fastq.gz -t 4"
    exit 1
}

# Check if no arguments were provided
if [ "$#" -eq 0 ]; then
    echo "Error: No arguments provided."
    usage
fi

# Parse command-line arguments
multi_sample=false
single_sample=false
name=""
fastq_file=""
fastqdir=""
threads=1  # Default number of threads

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -m|--multi)
            multi_sample=true
            shift
            ;;
        -s|--single)
            single_sample=true
            shift
            ;;
        -n|--name)
            name="$2"
            shift 2
            ;;
        -f|--fastqdir)
            fastqdir="$2"
            shift 2
            ;;
        -f|--fastqfile)
            fastq_file="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Invalid argument '$1'."
            usage
            ;;
    esac
done

# Ensure mutually exclusive options are enforced
if [ "$multi_sample" = true ] && [ "$single_sample" = true ]; then
    echo "Error: You cannot select both multi-sample and single-sample modes."
    usage
fi

if [ "$multi_sample" = true ]; then
    # Ensure required arguments for multi-sample mode
    if [ -z "$name" ] || [ -z "$fastqdir" ]; then
        echo "Error: Multi-sample mode requires both a run name and a fastq directory."
        usage
    fi
    # Ensure fastqdir exists
    if [ ! -d "$fastqdir" ]; then
        echo "Error: Directory '$fastqdir' does not exist."
        exit 1
    fi
    # Create output directory if it doesn't exist
    output_dir="data"
    mkdir -p "$output_dir"
    echo "Running in multi-sample mode..."
    files=($(find "$fastqdir" -type d -name "barcode*"))
    for dir in "${files[@]}"; do
        bc=$(basename "$dir")
        cat "$dir"/*fastq.gz > "data/${name}_${bc}.fastq.gz"
    done

elif [ "$single_sample" = true ]; then
    # Ensure required arguments for single-sample mode
    if [ -z "$name" ] || [ -z "$fastq_file" ]; then
        echo "Error: Single-sample mode requires both a sample name and a fastq file."
        usage
    fi
    # Ensure fastq file exists
    if [ ! -f "$fastq_file" ]; then
        echo "Error: Fastq file '$fastq_file' does not exist."
        exit 1
    fi
    # Create output directory if it doesn't exist
    output_dir="data"
    mkdir -p "$output_dir"
    echo "Running in single-sample mode for sample: $name"
    cp "$fastq_file" "data/${name}.fastq.gz"
fi

# Create the config directory if it doesn't exist
mkdir -p config
CONFIG_FILE="config/config.yaml"

# Write the static parts of the config file
echo "genome: resources/hg38.fa.gz" > $CONFIG_FILE
echo "threads: $threads" >> $CONFIG_FILE
echo "" >> $CONFIG_FILE

# **Handling Multi-Sample Mode**
if [ "$multi_sample" = true ]; then
    echo "samples:" >> $CONFIG_FILE
    sample_index=1
    for dir in "$fastqdir"/barcode*; do
        if [ -d "$dir" ]; then  # Ensure it's a directory
            bc=$(basename "$dir")
            output_file="data/${name}_${bc}.fastq.gz"
            cat "$dir"/*fastq.gz > "$output_file"
            echo "  sample${sample_index}: $(basename "$output_file" .fastq.gz)" >> $CONFIG_FILE
            ((sample_index++))
        fi
    done

elif [ "$single_sample" = true ]; then
    echo "samples:" >> $CONFIG_FILE
    sample_name=$(basename "$fastq_file" .fastq.gz)
    cp "$fastq_file" "data/${name}.fastq.gz"
    echo "  sample1: ${name}" >> $CONFIG_FILE
fi


echo "Configuration file generated at $CONFIG_FILE"

# Run analysis
echo "Running analysis with $threads threads..."
snakemake runall --cores "$threads" --use-conda --rerun-incomplete


