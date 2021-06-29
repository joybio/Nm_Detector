configfile:  "Nm_Detector.yaml"
#__date__ = "2021-5-27"
#__author__ = "Junbo Yang"
#__email__ = "yang_junbo_hi@126.com"
#__license__ = "PKU.jia.group"


#Snakefile
#Snakefile中的每一个rule其实都可以看作是一个简单的shell脚本，通过Snakefile将多个rule组织在一起并按照
#我们定义的顺序来执行。

#from snakemake.utils import makedirs
#from snakemake.utils import listfiles

#import numpy as np
#import os
reads = ["1","2"]
sample = config["samples"]
index = config["hisat2_idx"]
strandness = config["strandness"]
chrom_size = config["chrom_size"]
group = config["group"]
#samples = ["input-1","input-2","OED-1","OED-2"]
rule all:
	## LOCAL ##
#	'''
#	Defines the target rule by specifying all final output files.
#	Additionally, the cluster_logs dir is deleted if
#	snakemake was run locally.
#	'''
	input: 
		#expand(config["results_dir"] + "/{sample}.R1.test.uniq.sorted.3pend",sample=sample)
		config["results_dir"] + "/confident.Nm_sites.xls",
		expand(config["results_dir"] + "/{sample}.R1.test.uniq.3pend.bw",sample=sample)
#	output:
#		expand(config["results_dir"] + "/{sample}.test.uniq.3pend.bed",sample=sample)

#rule cluster:
#    input:
#	script = 'python/dbscan.py',
#	path   = 'data_files/umap/{sample}_umap.csv'
#    output:
#	path = 'output/{sample}'
#    shell:
#	"python {input.script} -data {input.path} -eps '0.3' -min_samples '10' "


#-------------------------------------------------------------------------------
# cutadapt rule1: Dependency packages - cutadapt
#-------------------------------------------------------------------------------
rule cutadapt:
##定义第一条规则。这里的rule可视为snakemake定义的关键字，concat使我们自定义的这一步任务的名称
	input:  raw_R1 = config["input_dir"] + "/{sample}_combined.R1.fastq.gz",
		raw_R2 = config["input_dir"] + "/{sample}_combined.R2.fastq.gz"
## input同样是snakemake的关键字，定义了在这个任务中的输入文件
#"{file}.txt", file=["hello", "world"]) #expand是一个snakemake定义的替换命令
	output: clean_R1 = config["results_dir"] + "/{sample}.trim.R1.fq.gz",
		clean_R2 = config["results_dir"] + "/{sample}.trim.R2.fq.gz"
## output也是snakemake的关键字，定义输出结果的保存文件
	message: "starting cutadaptor ..."
#message。使用message参数可以指定每运行到一个rule时，在终端中给出提示信息，eg.message: "starting mapping ..."。
	params:
#指定程序运行的参数，eg.params: cat="-n",调用方法为{params.cat}。
		error_rate=config['error_rate'],
		minLen=config['min_length'],
		overlap=config['overlap']
#	log:
#用来指定生成的日志文件，eg.log: "logs/concat.log"。
#		config["log_dir"] + "{sample}.trimmed.log"
	shell:
##指定执行方式，主要是三种：shell, run, script。分别对应shell命令,python命令,自己编写的命令
#在run的缩进区域里面可以输入并执行python代码。
#用来执行指定脚本，eg.scripts: "rm_dup.py"
#在shell命令中直接调用config文件中的内容的话，不需要引号，如config[a]而不是config["a"]。
		'''
		cutadapt -j 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
			-A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
			-e {params.error_rate} -O {params.overlap} -m {params.minLen} \
			-o {output.clean_R1} -p {output.clean_R2} {input.raw_R1} {input.raw_R2} \
		'''
#-------------------------------------------------------------------------------
# map rule2: Dependency packages - hisat2
#-------------------------------------------------------------------------------

rule hisat2_mapping:
	input:
		clean_R1 = config["results_dir"] + "/{sample}.trim.R1.fq.gz",
		clean_R2 = config["results_dir"] + "/{sample}.trim.R2.fq.gz" 
##因为这里的变量名是局部变量，仅在每个rule定义的子部分起作用，因此不同子部分可以用相同变量名，但内容取决于'='后面的部分
	output:
		config["results_dir"] + "/{sample}.R1.test.sam"
 #temp()用于标记临时文件，最后会自动删除，以免sam文件过大挤兑储存空间
	log:
		config["log_dir"] + "/{sample}.hisat2.log"
	message:
		"start hisat2 mapping..."
	shell:
		'''
		hisat2 -p 20 -k 1 --rna-strandness {strandness} --pen-noncansplice 1000000 \
		-x {index} \
		-U {input.clean_R1} \
		-S {output} \
		 > {log} 2>&1
		'''
#-------------------------------------------------------------------------------
# extract_unique_mapped_reads rule3: Dependency packages - None
#-------------------------------------------------------------------------------

rule extract_unique_mapped_reads:
	input:
		config["results_dir"] + "/{sample}.R1.test.sam"
	output:
		config["results_dir"] + "/{sample}.R1.test.uniq.sam"
	message:
		"start extract unique mapped reads..."
	shell:
		"cat {input} | grep NH:i:1 > {output}"

#-------------------------------------------------------------------------------
# extract_3pend_sites rule4: Dependency packages - None
#-------------------------------------------------------------------------------

rule extract_3pend_sites:
	input:
		config["results_dir"] + "/{sample}.R1.test.uniq.sam"
	output:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend"
	message:
		"start extract 3pend_sites..."
	shell:
		"python scripts/extract.3pend.py -i {input} -o {output}"

#-------------------------------------------------------------------------------
# count_3pends rule5: Dependency packages - None
#-------------------------------------------------------------------------------

rule count_3pends:
	input:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend"
	output:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend.bed",
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend.chr.bed"
	shell:
		"""
		python scripts/count_3pend.py -i {input} -o {output[0]} -f {output[1]}
		"""

#-------------------------------------------------------------------------------
# filter chromosome by hg38.chrom.sizes rule6: Dependency packages - None
#-------------------------------------------------------------------------------

rule filter_chr_by_UCSC:
	input:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend.chr.bed"
	output:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend.chr.filter.bed"
	shell:
		"""
		python scripts/extract.chr.py -b {input} -o {output}
		"""
#-------------------------------------------------------------------------------
# sort_3pends by chromosome rule7: Dependency packages - None
#-------------------------------------------------------------------------------

rule sort_3pends:
	input:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend.chr.filter.bed"
	output:
		config["results_dir"] + "/{sample}.R1.test.uniq.sorted.3pend.chr.filter.bed"
	shell:
		"""
		sort -k1,1 -k2,2n {input} > {output}
		"""

#-------------------------------------------------------------------------------
# bedToBigWig rule 8: Dependency packages - bedGraphToBigWig
#-------------------------------------------------------------------------------
#prepare bw files for IGV
rule bedToBigWig:
	input:
		bedfile = config["results_dir"] + "/{sample}.R1.test.uniq.sorted.3pend.chr.filter.bed",
		chrom_size = chrom_size
	output:
		config["results_dir"] + "/{sample}.R1.test.uniq.3pend.bw"
	shell:
		"""
		bedGraphToBigWig {input.bedfile} {input.chrom_size} {output}
		"""

#-------------------------------------------------------------------------------
# merge_3pends_matrix rule9: Dependency packages - None
#-------------------------------------------------------------------------------

rule merge_3pends_matrix:
	input:
		expand(config["results_dir"] + "/{sample}.R1.test.uniq.3pend.bed",sample = sample)
	output:
		config["results_dir"] + "/3pend_counts_matrix.xls"
	shell:
		"""
		python scripts/merge_3pend_count_matrix.py {input}
		"""

#-------------------------------------------------------------------------------
# merge_3pends_matrix rule10: Dependency packages - R: DESeq2 dplyr tidyr
#-------------------------------------------------------------------------------

rule Nm_finder:
	input: 
		matrixs = config["results_dir"] + "/3pend_counts_matrix.xls",
		group = group
	output:
		config["results_dir"] + "/pre.Nm_sites.xls"
	shell:
		"""
		Rscript scripts/DEseq2.reads_end.r {input.matrixs} {input.group} {output}
		"""

#-------------------------------------------------------------------------------
# filter_known_sites rule11: Dependency packages - None
#-------------------------------------------------------------------------------

rule filter_known_stop_sites:
	input:
		preNm = config["results_dir"] + "/pre.Nm_sites.xls",
		gtf = config["gtf_file"],
		tRNA = config["tRNA_database"]
	output:
		config["results_dir"] + "/confident.Nm_sites.xls"

	shell:
		"""
		python scripts/filter_known_stop_sites.py -b {input.preNm} -g {input.gtf} -t {input.tRNA} -o {output}
		"""

