rule sniffles2:
	input:
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam.bai",
		ref=config["genome"],
	output:
		"results/{sample}.minimap2.sniffles2.vcf",
	threads: 
		config["threads"]
	log:
		"logs/{sample}.minimap2.sniffles2.log"
	benchmark:
		"benchmarks/{sample}.sniffles2.tsv"
	conda:
		"../envs/sniffles2.yaml"
	shell:
		"sniffles --input {input.bam} --output-rnames --reference {input.ref} --vcf {output} 2> {log}"

rule sniffles2_table:
	input:
		"results/{sample}.minimap2.sniffles2.vcf"
	output:
		"results/{sample}.minimap2.sniffles2.tsv",
	threads: 
		config["threads"]
	log:
		"logs/{sample}.sniffles2_table.log"
	benchmark:
		"benchmarks/{sample}.sniffles2_table.tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\t[%DR]\t[%DV]\t%INFO/AF' -o {output} {input} 2>{log}
		"""

rule filter_sniffles2_tsv:
	input:
		"results/{sample}.minimap2.sniffles2.tsv"
	output:
		"results/{sample}.minimap2.sniffles2.filtered.tsv"
	threads:
		config["threads"]
	log:
		"logs/{sample}.filter_sniffles2_tsv"
	benchmark:
		"benchmarks/{sample}.filter_sniffles2_tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"awk '$4 == DEL' {input} > {output} 2> {log}"

rule svim2:
	input:
		ref=config["genome"],
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam"
	output:
		"results/{sample}.svim2.vcf"
	threads: 1
	conda:
		"../envs/svim2.yaml"
	benchmark:
		"benchmarks/{sample}.svim2.tsv"
	log:
		"logs/{sample}.svim2.log"
	params:
		wd="results/{sample}_svim2"
	shell:
		"mkdir -p {params.wd} && svim alignment --min_sv_size 10 --sample {wildcards.sample} --min_mapq 20 --minimum_score 2 --minimum_depth 2 {params.wd} {input.bam} {input.ref} &>{log} && mv {params.wd}/variants.vcf {output}"


rule cuteSV:
	input:
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam.bai",
		ref=config["genome"],
	output:
		"results/{sample}/{sample}.minimap2.cuteSV.vcf.gz"
	params:
		wd="results/{sample}_cutesv",
		vcf="results/{sample}/{sample}.minimap2.cuteSV.vcf"
	threads:
		config["threads"]
	conda:
		"../envs/cutesv.yaml"
	log:
		"logs/{sample}.minimap2.cuteSV.log"
	benchmark:
		"benchmarks/{sample}.minimap2.cuteSV.tsv"
	shell:
		"""
		mkdir -p {params.wd} && \
		cuteSV \
		--min_mapq 20 \
		--min_support 2 \
		--min_size 50 \
		--sample {wildcards.sample} \
		-t {threads} \
		--max_cluster_bias_INS 100 \
		--diff_ratio_merging_INS 0.3 \
		--max_cluster_bias_DEL 100 \
		--diff_ratio_merging_DEL 0.3 \
		{input.bam} \
		{input.ref} \
		{params.vcf} \
		{params.wd} 2>{log} && \
		bgzip {params.vcf}
		"""


rule custom_tool:
	input:
		script="workflow/scripts/run_v2.py",
		bam="alignments/{sample}.minimap2.srt.bam",
		bai="alignments/{sample}.minimap2.srt.bam.bai"
	output:
		tsv="results/{sample}.custom.tsv"
	log:
		"logs/{sample}.custom_tool.log"
	benchmark:
		"benchmarks/{sample}.custom_tool.tsv"
	conda:
		"../envs/custom_tool.yaml"
	shell:
		"python {input.script} --input_bam {input.bam} --output_file {output.tsv}"

