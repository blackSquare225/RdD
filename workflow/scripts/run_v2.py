import argparse
import pysam
from collections import defaultdict
import numpy as np
import csv

def detect_deletions(aligned_bam, min_size, max_size, min_mapping_quality=25, exclude_flag=3844):
    deletions = []
    try:
        with pysam.AlignmentFile(aligned_bam, "rb") as bam:
            for read in bam.fetch():
                if read.flag & exclude_flag:
                    continue
                if read.mapping_quality < min_mapping_quality:
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


def summarize_clusters(clusters, bam_file, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow([
            "Chromosome", "Start", "End", "Length",
            "StdDev_Start", "StdDev_End", "StdDev_Length",
            "Read_Support", "Total_Reads_in_Region"
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
            min_total_deviation = float("inf")
            representative_deletion = None
            for deletion in cluster:
                start, end, length = deletion[1], deletion[2], deletion[3]
                total_deviation = sum(
                    abs(start - other[1]) + abs(end - other[2]) + abs(length - other[3])
                    for other in cluster
                )
                if total_deviation < min_total_deviation:
                    min_total_deviation = total_deviation
                    representative_deletion = deletion
            rep_start = representative_deletion[1]
            rep_end = representative_deletion[2]
            rep_length = representative_deletion[3]
            total_reads_in_region = sum(
                1 for read in bam.fetch(chrom, rep_start, rep_end)
            )
            writer.writerow([
                chrom, rep_start, rep_end, rep_length,
                round(std_start, 2), round(std_end, 2), round(std_length, 2),
                read_support, total_reads_in_region
            ])
    bam.close()

# Main CLI Function
def main():
    parser = argparse.ArgumentParser(
        description="Detect Structural Variant Deletions from BAM files."
    )
    parser.add_argument("-i", "--input_bam", required=True, help="Path to the input BAM file.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output TSV file.")
    parser.add_argument("--min_size", type=int, default=300, help="Minimum deletion size to detect (default: 300 bp).")
    parser.add_argument("--max_size", type=int, default=4000, help="Maximum deletion size to detect (default: 4000 bp).")
    parser.add_argument("--min_mapping_quality", type=int, default=25, help="Minimum mapping quality (default: 25).")
    parser.add_argument("--exclude_flag", type=int, default=3844, help="SAM flag to exclude reads (default: 3844).")
    parser.add_argument("--tolerance", type=int, default=100, help="Tolerance between breakpoints for clustering deletions as the same event (default: 100 bp).")
    
    args = parser.parse_args()

    # Run the workflow
    print("Detecting deletions...")
    deletions = detect_deletions(
        args.input_bam, args.min_size, args.max_size, args.min_mapping_quality, args.exclude_flag
    )
    print(f"Detected {len(deletions)} deletions.")
    
    print("Clustering deletions...")
    clusters = cluster_deletions(deletions, tolerance=args.tolerance)
    print(f"Generated {len(clusters)} clusters.")
    
    print("Summarizing clusters...")
    summarize_clusters(clusters, args.input_bam, args.output_file)
    print(f"Results written to {args.output_file}")

if __name__ == "__main__":
    main()
