# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Sep-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq
# Aim: Snakemake workflow for mapping QC
# ################################### IMPORT ##################################### #


import os
import sys
import re
import glob
from snakemake.utils import read_job_properties
from multiprocessing import cpu_count
library_path = os.path.abspath(workflow.basedir + "/../library/")
sys.path.append(library_path)
import utility
# ################################### GENERAL FUNCTIONS ############################# #


def parse_mapping(WORK_DIR, TITLE):
	"""
	"""
	mapping_bam_List = []
	pre_processed_fastq_path_Dict = {}
	mapping_bam_Dict = {}

	mapping_bam_List = glob.glob(WORK_DIR + "/" + TITLE + "/**/map_star/*.trim_fastp.map_star.bam", recursive=True)
	for each_bam_Path in mapping_bam_List:
		##
		if utility.is_file_exist(each_bam_Path) is True:
			#
			each_sample = os.path.basename(each_bam_Path).split(".trim_fastp", 1)[0]
			mapping_bam_Dict[each_sample] = each_bam_Path
		else:
			#if utility.is_file_exist(each_fastq_path) is True:
			pass
	else:
		##for each_fastq_path in pre_processed_fastq_path_List:
		pass
	return mapping_bam_Dict
# ################################### PARAMTERS FUNCTIONS ########################## #


def build_kraken_layout_kraken_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		kraken_layout_kraken_command = ""
	else:
		#
		kraken_layout_kraken_command = WORK_DIR + "/" + TITLE + "/{sample}/map_star/{sample}.R2.trim_fastp.unmap_star.fastq.gz --paired".format(sample=sample)

	return kraken_layout_kraken_command


def build_sortmerna_layout_sortmerna_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		sortmeRNA_layout_command = ""
	else:
		#
		sortmeRNA_layout_command = "--reads " + WORK_DIR + "/" + TITLE + "/{sample}/map_star/{sample}.R2.trim_fastp.unmap_star.fastq.gz ".format(sample=sample)

	return sortmeRNA_layout_command


def build_dupradar_layout_dupradar_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		dupradar_layout_dupradar_command = " paired=no "
	else:
		#
		dupradar_layout_dupradar_command = " paired=yes "

	return dupradar_layout_dupradar_command

# ################################### CONFIGURATION ################################ #


# +++++++++++++++++++++++++++++++++++++
# MAIN PATH
Bash_Script = os.path.abspath(workflow.basedir + "/../Bash_script")
R_Script_path = os.path.abspath(workflow.basedir + "/../R_script")
Python_Script_path = os.path.abspath(workflow.basedir + "/../python_script")
Script_Path = os.path.abspath(workflow.basedir + "/../template")
# ------------------------------------
# +++++++++++++++++++++++++++++++++++++
# GENERAL
config_general_Dict = config["GENERAL"]
TITLE = config_general_Dict["TITLE"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
# DIRECTORY
config_directory_Dict = config["DIRECTORY"]
WORK_DIR = utility.fix_path(config_directory_Dict["WORK_DIR"])
DATA_DIR = utility.fix_path(config_directory_Dict["DATA_DIR"])
TEMP_DIR = utility.fix_path(config_directory_Dict["TEMP_DIR"])
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
# DATA
config_data_Dict = config["DATA"]
LAYOUT = config_data_Dict["LAYOUT"].lower()
SAMPLE_DELIMITER = config_data_Dict["SAMPLE_DELIMITER"]
PLATFORM = config_data_Dict["PLATFORM"].lower()
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
# REFERENCE
config_reference_Dict = config["REFERENCE"]
HOST = list(config_reference_Dict)[0]
KRAKEN2_DB = config_reference_Dict[HOST]["KRAKEN2_DB"]
STAR_DB = config_reference_Dict[HOST]["STAR_DB"]
STAR_FUSION_DB = config_reference_Dict[HOST]["STAR_FUSION_DB"]
GENOME_FASTA = config_reference_Dict[HOST]["GENOME_FASTA"]
KRONA_DB = config_reference_Dict[HOST]["KRONA_DB"]
SORTMERNA_DB = config_reference_Dict[HOST]["SORTMERNA_DB"]
GENE_BED = config_reference_Dict[HOST]["GENE_BED"]
GENE_GTF = config_reference_Dict[HOST]["GENE_GTF"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
# CLUSTER PRAMETERS
#PROCESSORS = workflow.cores
PROCESSORS = max(2, cpu_count() - 4)
#MEMORY = 14096
MODE = "DEVELOPMENT"
# ------------------------------------
# ################################### WILDCARDS ################################ #


mapping_bam_Dict = parse_mapping(WORK_DIR, TITLE)
#print(mapping_bam_Dict)


mapping_List = []
mapping_qc_List = []

for sample in mapping_bam_Dict:
	##
	mapping_List.append(WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam".format(title=TITLE, sample=sample))
	#
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.flagstat_samtools.txt".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.mosdepth.global.dist.txt".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star_fastqc.html".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.unmap_star.classify_kraken.txt".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.unmap_star.classify_kraken.html".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.unmap_star.classify_sortmeRNA.txt".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.bam_stat.txt".format(title=TITLE, sample=sample))
	mapping_qc_List.append(WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star_dupRadar_drescatter.png".format(title=TITLE, sample=sample))
else:
	##for sample in pre_processed_fastq_path_Dict:
	pass
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		mapping_qc_List
# ################################### PIPELINE RULES ########################### #


rule samtools:
	"""
	"""
	input:
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
	output:
		flagstat_samtools = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.flagstat_samtools.txt",
		idxstats_samtools = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.idxstats_samtools.txt",
	threads: PROCESSORS
	message: "samtools: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	run:
		sample_name = os.path.basename(input.map_bam).split(".bam", 1)[0]
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/samtools
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			#GENERAL_TAG={sample_name}.filter_samtools.markduplicate_picard.index_samtools
			GENERAL_TAG={sample_name}

			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "samtools: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "###################################- FLAGSTAT SAMTOOLS -#######################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			
			printf "%s\\n" "samtools flagstat \\
			--threads {threads} \\
			{input.map_bam} > {output.flagstat_samtools} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			samtools flagstat \\
			--threads {threads} \\
			{input.map_bam} > {output.flagstat_samtools} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.flagstat_samtools.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "samtools: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "###################################- IDXSTATS SAMTOOLS -#######################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			
			printf "%s\\n" "samtools idxstats \\
			{input.map_bam} > {output.idxstats_samtools} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			samtools idxstats \\
			{input.map_bam} > {output.idxstats_samtools} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.idxstats_samtools.log 2>&1
		""")


rule mosdepth:
	"""
	"""
	input:
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
	output:
		mosdepth_summary = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.mosdepth.summary.txt",
		mosdepth_coverage = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.mosdepth.global.dist.txt",
	threads: PROCESSORS
	message: "mosdepth: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	run:
		sample_name = os.path.basename(input.map_bam).split(".bam", 1)[0]
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/mosdepth
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			#GENERAL_TAG={sample_name}.filter_samtools.markduplicate_picard.index_samtools
			GENERAL_TAG={sample_name}

			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "mosdepth: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "###################################- COVERAGE MOSDEPTH -#######################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			
			printf "%s\\n" "mosdepth \\
			--threads {threads} \\
			$REPORT_PATH/{wildcards.sample}.trim_fastp.map_star \\
			{input.map_bam} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			mosdepth \\
			--threads {threads} \\
			$REPORT_PATH/{wildcards.sample}.trim_fastp.map_star \\
			{input.map_bam} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.coverage_mosdepth.log 2>&1

		""")


rule fastqc:
	"""
	"""
	input:
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
	output:
		fastqc_report = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star_fastqc.html",
	threads: PROCESSORS
	message: "fastqc: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	run:
		sample_name = os.path.basename(input.map_bam).split(".bam", 1)[0]
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/fastqc
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			#GENERAL_TAG={sample_name}.filter_samtools.markduplicate_picard.index_samtools
			GENERAL_TAG={sample_name}

			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "fastqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "###################################- COVERAGE MOSDEPTH -#######################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			
			printf "%s\\n" "fastqc \\
			-o $REPORT_PATH \\
			--dir $TEMP_PATH \\
			-f bam \\
			--threads {threads} \\
			{input.map_bam} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			fastqc \\
			-o $REPORT_PATH \\
			--dir $TEMP_PATH \\
			-f bam \\
			--threads {threads} \\
			{input.map_bam} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.qc_fastqc.log 2>&1

		""")


rule kraken:
	"""
	"""
	input:
		unmap_fastq = WORK_DIR + "/{title}/{sample}/map_star/{sample}.R1.trim_fastp.unmap_star.fastq.gz",
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
	output:
		kraken_report = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.unmap_star.classify_kraken.txt",
		krona_report = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.unmap_star.classify_kraken.html",
	threads: PROCESSORS
	message: "kraken: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	params:
		kraken_layout_kraken_command = build_kraken_layout_kraken_command,
	run:
		sample_name = os.path.basename(input.unmap_fastq).split(".fastq.gz", 1)[0]
		sample_name = sample_name.replace(".R1", "")
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/kraken
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			GENERAL_TAG={sample_name}


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "kraken: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "###################################- CLASSIFY KRAKEN -#########################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			
			printf "%s\\n" "kraken2 \\
			--threads {threads} \\
			--db {KRAKEN2_DB} \\
			--report {output.kraken_report} \\
			--gzip-compressed {input.unmap_fastq} \\
			{params.kraken_layout_kraken_command} \\
			1> $REPORT_PATH/${{GENERAL_TAG}}.report_kraken.txt \\
			2>> $LOG_PATH/${{GENERAL_TAG}}classify_kraken.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			kraken2 \\
			--threads {threads} \\
			--db {KRAKEN2_DB} \\
			--report {output.kraken_report} \\
			--gzip-compressed {input.unmap_fastq} \\
			{params.kraken_layout_kraken_command} \\
			1> $REPORT_PATH/${{GENERAL_TAG}}.report_kraken.txt \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "kraken: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "###################################- VISUAL KRONA -###########################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			
			printf "%s\\n" "ktImportTaxonomy \\
			-q 2 \\
			-t 3 \\
			-tax {KRONA_DB} \\
			-o {output.krona_report} \\
			{output.kraken_report} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			ktImportTaxonomy -q 2 -t 3 -tax {KRONA_DB} -o {output.krona_report} {output.kraken_report} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_kraken.visual_krona.log 2>&1

		""")




rule sortmeRNA:
	input:
		unmap_fastq = WORK_DIR + "/{title}/{sample}/map_star/{sample}.R1.trim_fastp.unmap_star.fastq.gz",
	output:
		sortmeRNA_report = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.unmap_star.classify_sortmeRNA.txt",
	threads: PROCESSORS
	message: "sortmeRNA: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	params:
		sortmerna_layout_sortmerna_command = build_sortmerna_layout_sortmerna_command,
	run:
		sample_name = os.path.basename(input.unmap_fastq).split(".fastq.gz", 1)[0]
		sample_name = sample_name.replace(".R1", "")
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/sortmeRNA
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			GENERAL_TAG={sample_name}

			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "kraken: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "###################################- CLASSIFY SORTMERNA -######################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			
			printf "%s\\n" "sortmerna \\
			--threads {threads} \\
			--workdir $TEMP_PATH \\
			--ref {SORTMERNA_DB}/silva-euk-28s-id98.fasta \\
			--ref {SORTMERNA_DB}/rfam-5.8s-database-id98.fasta \\
			--ref {SORTMERNA_DB}/silva-arc-16s-id95.fasta \\
			--ref {SORTMERNA_DB}/silva-euk-18s-id95.fasta \\
			--ref {SORTMERNA_DB}/rfam-5s-database-id98.fasta \\
			--ref {SORTMERNA_DB}/silva-bac-23s-id98.fasta \\
			--ref {SORTMERNA_DB}/silva-bac-16s-id90.fasta \\
			--ref {SORTMERNA_DB}/silva-arc-23s-id98.fasta \\
			--reads {input.unmap_fastq} \\
			{params.sortmerna_layout_sortmerna_command} \\
			--aligned $TEMP_PATH/${{GENERAL_TAG}} --fastx \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/${{GENERAL_TAG}}.log $REPORT_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.txt" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			sortmerna \\
			--threads {threads} \\
			--workdir $TEMP_PATH \\
			--ref {SORTMERNA_DB}/silva-euk-28s-id98.fasta \\
			--ref {SORTMERNA_DB}/rfam-5.8s-database-id98.fasta \\
			--ref {SORTMERNA_DB}/silva-arc-16s-id95.fasta \\
			--ref {SORTMERNA_DB}/silva-euk-18s-id95.fasta \\
			--ref {SORTMERNA_DB}/rfam-5s-database-id98.fasta \\
			--ref {SORTMERNA_DB}/silva-bac-23s-id98.fasta \\
			--ref {SORTMERNA_DB}/silva-bac-16s-id90.fasta \\
			--ref {SORTMERNA_DB}/silva-arc-23s-id98.fasta \\
			--reads {input.unmap_fastq} \\
			{params.sortmerna_layout_sortmerna_command} \\
			--aligned $TEMP_PATH/${{GENERAL_TAG}} --fastx \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			mv $TEMP_PATH/${{GENERAL_TAG}}.log $REPORT_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.txt \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.classify_sortmeRNA.log 2>&1
		""")



rule dupradar:
	"""
	https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html
	"""
	input:
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
	output:
		dupradar_scatterplot = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star_dupRadar_drescatter.png",
		dupradar_barplot = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star_dupRadar_drebp.png",
		dupradar_histogram = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star_dupRadar_ehist.png",
	threads: PROCESSORS
	message: "dupradar: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	params:
		dupradar_layout_dupradar_command = build_dupradar_layout_dupradar_command,
	run:
		sample_name = os.path.basename(input.map_bam).split(".bam", 1)[0]
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/dupradar
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			GENERAL_TAG={sample_name}

			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "###################################- DUPRADAR -################################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1

			printf "%s\\n" "Rscript {R_Script_path}/dupRadar.R \\
			{input.map_bam} \\
			{GENE_GTF} \\
			stranded=yes \\
			{params.dupradar_layout_dupradar_command} \\
			outdir=${{REPORT_PATH}}/ \\
			threads={threads} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			Rscript {R_Script_path}/dupRadar.R \\
			{input.map_bam} \\
			{GENE_GTF} \\
			stranded=yes \\
			{params.dupradar_layout_dupradar_command} \\
			outdir=${{REPORT_PATH}}/ \\
			threads={threads} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.duplicate_dupradar.log 2>&1
		""")


rule rseqc:
	"""
	"""
	input:
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
	output:
		bam_stat = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.bam_stat.txt",
		infer_experiment = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.infer_experiment.txt",
		read_distribution = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.read_distribution.txt",
		read_duplication = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.pos.DupRate.xls",
		junction_saturation = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.junctionSaturation_plot.r",
		junction_annotation = WORK_DIR + "/{title}/{sample}/report/map_star_qc/{sample}.trim_fastp.map_star.rseqc.junction_annotation.txt",

	threads: PROCESSORS
	message: "rseqc: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	run:
		sample_name = os.path.basename(input.map_bam).split(".bam", 1)[0]
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star_qc
			#mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star_qc
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star_qc
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{sample_name}/map_star_qc/rseqc
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#

			GENERAL_TAG={sample_name}


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "###################################- BAM STAT RSEQC -##########################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1

			printf "%s\\n" "bam_stat.py \\
			-i {input.map_bam} > {output.bam_stat} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			bam_stat.py -i {input.map_bam} > {output.bam_stat} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.bam_stat_rseqc.log 2>&1



			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "###################################- INFER EXPERIMENT RSEQC -##################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1

			printf "%s\\n" "infer_experiment.py -r {GENE_BED} -i {input.map_bam} > {output.infer_experiment} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1

			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			infer_experiment.py -r {GENE_BED} -i {input.map_bam} > {output.infer_experiment} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.infer_experiment_rseqc.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "###################################- READ DISTRIBUTION RSEQC -#################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1

			printf "%s\\n" "read_distribution.py -r {GENE_BED} -i {input.map_bam} > {output.read_distribution} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1

			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			read_distribution.py -r {GENE_BED} -i {input.map_bam} > {output.read_distribution} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_distribution_rseqc.log 2>&1



			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "###################################- JUNCTION SATURATION RSEQC -###############" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1

			printf "%s\\n" "junction_saturation.py -r {GENE_BED} -i {input.map_bam} -o $REPORT_PATH/{sample_name}.rseqc \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1

			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			junction_saturation.py -r {GENE_BED} -i {input.map_bam} -o $REPORT_PATH/{sample_name}.rseqc \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_saturation_rseqc.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "###################################- READ DUPLICATION RSEQC -##################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1

			printf "%s\\n" "junction_saturation.py -r {GENE_BED} -i {input.map_bam} -o $REPORT_PATH/{sample_name} > {output.junction_saturation} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1

			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			read_duplication.py -i {input.map_bam} -o $REPORT_PATH/{sample_name}.rseqc \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.read_duplication_rseqc.log 2>&1



			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "rseqc: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "###################################- JUNCTION ANNOTATION RSEQC -###############" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1

			printf "%s\\n" "junction_annotation.py -r {GENE_BED} -i {input.map_bam} -o $REPORT_PATH/{sample_name}.rseqc 2> {output.junction_annotation} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1

			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			junction_annotation.py -r {GENE_BED} -i {input.map_bam} -o $REPORT_PATH/{sample_name}.rseqc 2> {output.junction_annotation} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.junction_annotation_rseqc.log 2>&1


		""")


# ################################### FINITO ################################## #