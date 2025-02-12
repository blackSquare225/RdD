#rule mosdepth:
#	input:
#		bam="alignments/{sample}.minimap2.srt.bam",
#		bai="alignments/{sample}.minimap2.srt.bam.bai",
#		bed=config["interval"]
#	output:
#		"alignments/{sample}.mosdepth.regions.bed.gz"
#	params:
#		prefix="alignments/{sample}.mosdepth",
#	threads: 
#		config["threads"]
#	log:
#		"logs/{sample}.mosdepth.log"
#	benchmark:
#		"benchmarks/{sample}.mosdepth.tsv"
#	conda: 
#		"../envs/mosdepth.yaml"
#	shell:
#		"mosdepth -t {threads} --by {input.bed} {params.prefix} {input.bam} 2>{log}"


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


rule call_deletions:
	input:
		script="workflow/scripts/DELecter.py",
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam.bai"
	output:
		tsv="results/{sample}.results.tsv"
	params:
		min_size=config["min_size"],
		max_size=config["max_size"],
		mapq=config["mapq"],
		min_support=config["min_support"],
		min_len=config["min_len"],
		exclude_flag=config["exclude_flag"],
		tolerance=config["tolerance"]
	log:
		"logs/{sample}.DELecter.log"
	benchmark:
		"benchmarks/{sample}.DELecter.tsv"
	conda:
		"../envs/DELecter.yaml"
	shell:
		"""
		python {input.script} \
			--input_bam {input.bam} \
			--output_file {output.tsv} \
			--min_size {params.min_size} \
			--max_size {params.max_size} \
			--mapq {params.mapq} \
			--min_support {params.min_support} \
			--min_len {params.min_len} \
			--exclude_flag {params.exclude_flag} \
			--tolerance {params.tolerance} 2> {log}
		"""
