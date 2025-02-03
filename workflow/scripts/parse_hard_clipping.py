import pysam
import statistics

# Open the BAM file
bam_file = "run2_barcode24.minimap.srt.hmga2.DEL.srt.bam"  # replace with your BAM file path
bam = pysam.AlignmentFile(bam_file, "rb")

# Initialize a list to store read lengths
read_lengths = []

# Iterate through each read in the BAM file
for read in bam:
    read_lengths.append(read.query_length)

# Close the BAM file
bam.close()

# Calculate the mean read length
if read_lengths:
    mean_read_length = statistics.mean(read_lengths)
    print(f"Mean read length: {mean_read_length}")
else:
    print("No reads found in the BAM file.")


import pysam

# Open the BAM file
bam_file = "run2_barcode24.minimap.srt.hmga2.DEL.srt.bam"  # replace with your BAM file path
bam = pysam.AlignmentFile(bam_file, "rb")

# Initialize a counter for hard-clipped reads
hardclipped_count = 0

# Iterate through each read in the BAM file
for read in bam:
    # Check the CIGAR string for hard clipping
    if any(cigartuples[0] == 5 for cigartuples in read.cigartuples):  # 5 is the code for hard clipping
        hardclipped_count += 1

# Close the BAM file
bam.close()

# Print the number of hard-clipped reads
print(f"Number of hard-clipped reads: {hardclipped_count}")
