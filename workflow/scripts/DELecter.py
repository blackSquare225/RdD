import argparse
import pysam
from collections import defaultdict, Counter
import numpy as np
import csv

def detect_deletions(aligned_bam, min_size=300, max_size=4000, mapq=25, exclude_flag=3844):
    deletions = []
    try:
        with pysam.AlignmentFile(aligned_bam, "rb") as bam:
            for read in bam.fetch():
                if read.flag & exclude_flag:
                    continue
                if read.mapping_quality < mapq:
                    continue
                if not read.cigartuples:
                    continue
                for cigartype, length in read.cigartuples:
                    if cigartype == 2 and min_size <= length <= max_size:
                        start = read.reference_start
                        end = start + length
                        deletions.append((read.reference_name, start, end, length))
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not open file: {aligned_bam}")
    except ValueError as e:
        raise ValueError(f"Error processing BAM file: {e}")
    return deletions


def cluster_deletions(deletions, tolerance=100):
    deletions.sort(key=lambda x: (x[0], x[1], x[2]))
    parent = list(range(len(deletions)))
    rank = [0] * len(deletions)

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            if rank[root_x] > rank[root_y]:
                parent[root_y] = root_x
            elif rank[root_x] < rank[root_y]:
                parent[root_x] = root_y
            else:
                parent[root_y] = root_x
                rank[root_x] += 1

    for i in range(len(deletions)):
        for j in range(i + 1, len(deletions)):
            if deletions[i][0] != deletions[j][0]:
                break
            if abs(deletions[i][1] - deletions[j][1]) <= tolerance and abs(deletions[i][2] - deletions[j][2]) <= tolerance:
                union(i, j)

    clusters = defaultdict(list)
    for i in range(len(deletions)):
        root = find(i)
        clusters[root].append(deletions[i])

    return list(clusters.values())


def summarize_clusters(clusters, bam_file, output_file, mapq=25, exclude_flag=3844, min_len=1000, length_tolerance=0.2, min_read_support=100):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow([
            "Chromosome", "Start", "End", "Length",
            "StdDev_Start", "StdDev_End", "StdDev_Length",
            "Read_Support", "Total_Reads_in_Region", "VAF"
        ])

        for cluster in clusters:
            chrom = cluster[0][0]
            starts = [c[1] for c in cluster]
            ends = [c[2] for c in cluster]
            lengths = [c[3] for c in cluster]
            std_start = np.std(starts)
            std_end = np.std(ends)
            std_length = np.std(lengths)
            read_support = len(cluster)
            
            # Compute medians
            median_start = int(np.median(starts))
            median_end = int(np.median(ends))
            median_length = int(np.median(lengths))
            
            # Step 1: Separate deletions into groups based on length similarity
            deletion_groups = defaultdict(list)
            for deletion in cluster:
                grouped = False
                for key_length in deletion_groups:
                    if abs(deletion[3] - key_length) / key_length <= length_tolerance:
                        deletion_groups[key_length].append(deletion)
                        grouped = True
                        break
                if not grouped:
                    deletion_groups[deletion[3]].append(deletion)
            
            # Step 2: Prioritize the group with at least one deletion ≥ min_len
            large_deletion_groups = {k: v for k, v in deletion_groups.items() if k >= min_len}
            if large_deletion_groups:
                best_group = max(large_deletion_groups.values(), key=len)
            else:
                best_group = max(deletion_groups.values(), key=len)  # If no large deletion, use most frequent size
            
            # Step 3: Choose the best deletion (highest support, closest to median)
            best_group.sort(key=lambda d: (-len(d), abs(d[1] - median_start) + abs(d[2] - median_end)))
            representative_deletion = best_group[0]
            
            rep_start = representative_deletion[1]
            rep_end = representative_deletion[2]
            rep_length = representative_deletion[3]
            
            total_reads_in_region = sum(
                1 for read in bam.fetch(chrom, rep_start, rep_end)
                if not (read.flag & exclude_flag) and read.mapping_quality >= mapq
            )
            


            # Apply filtering criteria before writing to the file
            if read_support > min_read_support and rep_length > min_len:
                writer.writerow([
                    chrom, rep_start, rep_end, rep_length,
                    round(std_start, 2), round(std_end, 2), round(std_length, 2),
                    read_support, total_reads_in_region , round(read_support/total_reads_in_region, 3) 
                ])
    bam.close()

# Main CLI Function
def main():
    parser = argparse.ArgumentParser(
        description="Detect Structural Variant Deletions from BAM files."
    )
    parser.add_argument("-i", "--input_bam", required=True, help="Path to the input BAM file.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output TSV file.")
    parser.add_argument("--min_size", type=int, default=300, help="Minimum deletion size (in bp) to take into account. (default: 300)")
    parser.add_argument("--max_size", type=int, default=4000, help="Maximum deletion size (in bp) to take into account. (default: 4000)")
    parser.add_argument("--mapq", type=int, default=25, help="Reads with mapping quality lower than this value will be ignored. (default: 25)")
    parser.add_argument("--min_support", type=int, default=100, help="Minimum number of supporting reads for a DEL to be reported. (default 100)")
    parser.add_argument("--min_len", type=int, default=1000, help="Minimum DEL length (in bp) to be reported. (default 1000)")
    parser.add_argument("--exclude_flag", type=int, default=3844, help="Bitwise SAM flag to exclude. (default 3844)")
    parser.add_argument("--tolerance", type=int, default=100, help="Tolerance between breakpoints (in bp) for clustering deletions as the same event. (default: 100)")
    
    args = parser.parse_args()

    # Run the workflow
    print("Detecting deletions...")
    deletions = detect_deletions(
        args.input_bam, 
        min_size=args.min_size, 
        max_size=args.max_size, 
        mapq=args.mapq, 
        exclude_flag=args.exclude_flag
    )
    print(f"Detected {len(deletions)} deletions.")
    
    print("Clustering deletions...")
    clusters = cluster_deletions(
        deletions, 
        tolerance=args.tolerance
        )
    print(f"Generated {len(clusters)} clusters.")
    
    print("Summarizing clusters...")
    summarize_clusters(clusters, 
        args.input_bam, 
        args.output_file,
        mapq=args.mapq, 
        exclude_flag=args.exclude_flag, 
        min_len=args.min_len, 
        min_read_support=args.min_support 
        )
    print(f"Results written to {args.output_file}")

if __name__ == "__main__":
    main()
