rule minimap2:
	input:
		config["genome"],
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

rule mosdepth:
	input:
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam.bai",
		bed=config["interval"]
	output:
		"alignments/{sample}.mosdepth.regions.bed.gz"
	params:
		prefix="alignments/{sample}.mosdepth",
	threads: 
		config["threads"]
	log:
		"logs/{sample}.mosdepth.log"
	benchmark:
		"benchmarks/{sample}.mosdepth.tsv"
	conda: 
		"../envs/mosdepth.yaml"
	shell:
		"mosdepth -t {threads} --by {input.bed} {params.prefix} {input.bam} 2>{log}"


#rule write_stats_cov_per_sample:
#	input:
#		hmga2="alignments/{sample}.mosdepth_hmg2a.regions.bed.gz", 
#		DEL="alignments/{sample}.mosdepth_hmg2a_DEL.regions.bed.gz"
#	output:
#		temp("alignments/{sample}_stats_cov.tsv")	
#	threads: 1
#	log:
#		"logs/{sample}.write_stats_cov.log"
#	shell:
#		"""
#		paste <(zcat {input.hmga2}) <(zcat {input.DEL}) | \
#		awk -v sample="{wildcards.sample}" 'OFS=FS="\t"''{{print $1,$2,$3,sample,$4,$8,$8/$4}}' | tr '.' ',' > {output}
#		"""

#rule write_stats_cov_allsamples:
#	input:
#		expand(f"alignments/{{sample}}_stats_cov.tsv", sample=config["samples"].values())
#	output:	
#		"alignments/stats_cov.tsv"
#	threads: 1
#	shell:
#		"""
#		echo -e "CHR\tSTART\tEND\tSAMPLE\tCOV_HMGA2\tCOV_HMGA2del\tRATIO" | cat - <(cat {input}) | sort -k4,4rn > {output}
#		"""