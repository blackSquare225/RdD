rule get_reference_genome:
	output:
		"resources/hg38.fa.gz"
	threads: 1 
	log:
		"logs/get_reference_genome.log"
	benchmark:
		"benckmarks/get_reference_genome.tsv"
	shell:
		"""
		wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
		"""

rule minimap2_index:
	input:
		"resources/hg38.fa.gz"
	output:
		"resources/hg38.mmi"
	threads:
		config["threads"]
	log:
		"logs/minimap2_index.log"
	benchmark:
		"benckmarks/minimap2_index.tsv"
	shell:
		"minimap2 -t {threads} -d {output} {input}"

rule minimap2_align:
	input:
		"resources/hg38.mmi",
		"data/{sample}.fastq.gz",
	output:
		"alignments/{sample}.minimap2.srt.bam"
	threads: 
		config["threads"]
	log:
		"logs/{sample}.minimap2.log"
	benchmark:
		"benchmarks/{sample}.minimap2.tsv"
	conda:
		"../envs/minimap2.yaml"
	shell:
		"""
		minimap2 --MD -ax map-ont -t {threads} {input} 2> {log} | \
		samtools sort -@ {threads} -o {output} 2>> {log}
		"""
		
rule bam_index:	
	input:
		"alignments/{sample}.minimap2.srt.bam"
	output:
		"alignments/{sample}.minimap2.srt.bam.bai"	
	threads:
		config["threads"]
	log:
		"logs/{sample}.bam_index.log"
	benchmark:
		"benchmarks/{sample}.bam_index.tsv"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools index -@ {threads} {input} 2>{log}"
