import pysam
import argparse

def extract_deletions(bam_file, region, output_file, exclude_flag):
    """
    Extract deletions from a BAM file within a specified region, excluding reads with a specific flag.

    Args:
        bam_file (str): Path to the BAM file or URL.
        region (str): Genomic region in the format "chr:START-END".
        output_file (str): Path to the output file.
        exclude_flag (int): The flag value to exclude reads with this flag.
    """
    # Parse the region
    chrom, pos = region.split(':')
    start, end = map(int, pos.split('-'))

    # Open the BAM file (either local or remote URL)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Initialize variables
        unique_reads = set()  # Set to store unique read names
        deletions = []

        # Iterate through the reads in the region
        for read in bam.fetch(chrom, start, end):
            # Exclude reads with the specified flag
            if read.flag == exclude_flag:
                continue

            unique_reads.add(read.query_name)  # Add read name to the set

            # Parse the CIGAR string for deletions
            found_deletion = False
            for cigar_op, length in read.cigar:
                if cigar_op == 2:  # CIGAR operation 2 corresponds to deletion
                    if 500 <= length <= 3000:
                        deletions.append((read.query_name, length))
                        found_deletion = True

        # Count the number of unique reads
        total_reads = len(unique_reads)

        # If no deletions were found, append a single line with NA and 0
        if not deletions:
            deletions.append(('NA', 0))

    # Write the output
    with open(output_file, "w") as out:
        # Write header
        out.write("ReadName\tDeletionLength\tTotalReads\n")
        
        # Write the deletions or the NA line
        for read_name, del_length in deletions:
            out.write(f"{read_name}\t{del_length}\t{total_reads}\n")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Extract deletions from a BAM file within a specified genomic region, excluding reads with a specific flag."
    )
    parser.add_argument(
        "-b", "--bam", required=True, help="Path to the input BAM file or URL."
    )
    parser.add_argument(
        "-r", "--region", required=True, help="Genomic region in the format 'chr:START-END'."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to the output TSV file."
    )
    parser.add_argument(
        "-f", "--exclude_flag", type=int, default=1796, help="Flag to exclude reads with this specific flag. Default is 1796."
    )

    # Parse arguments
    args = parser.parse_args()

    # Call the main extraction function
    extract_deletions(args.bam, args.region, args.output, args.exclude_flag)


if __name__ == "__main__":
    main()
