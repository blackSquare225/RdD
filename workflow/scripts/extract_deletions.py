import pysam
import argparse

def extract_reads_with_deletions(input_bam, output_bam, min_deletion_size, max_deletion_size):
    # Open input BAM file
    bamfile = pysam.AlignmentFile(input_bam, "rb")
    
    # Open output BAM file
    outfile = pysam.AlignmentFile(output_bam, "wb", template=bamfile)

    for read in bamfile.fetch():
        # Get the CIGAR string of the read
        cigartuples = read.cigartuples
        
        if cigartuples is not None:
            # Iterate through the CIGAR operations
            for op, length in cigartuples:
                # Check for deletions (CIGAR operation 2 corresponds to 'D')
                if op == 2 and min_deletion_size <= length <= max_deletion_size:
                    outfile.write(read)
                    break  # No need to check further, we found a qualifying deletion
    
    # Close the BAM files
    bamfile.close()
    outfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads with deletions from a BAM file.")
    parser.add_argument("input_bam", help="Path to the input BAM file")
    parser.add_argument("output_bam", help="Path to the output BAM file")
    parser.add_argument("min_deletion_size", type=int, help="Minimum deletion size")
    parser.add_argument("max_deletion_size", type=int, help="Maximum deletion size")

    args = parser.parse_args()

    extract_reads_with_deletions(args.input_bam, args.output_bam, args.min_deletion_size, args.max_deletion_size)
