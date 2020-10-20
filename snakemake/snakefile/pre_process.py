# ################################### INFO ####################################### #
# Author: Amir Shams
# Date: Sep-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq
# Aim: Snakemake workflow for pre_process fastq files
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
# ################################### GENERAL FUNCTIONS ########################## #


def parse_datadir(DATA_DIR, SAMPLE_DELIMITER, LAYOUT, FORWARD_DELIMITER_LIST, REVERSE_DELIMITER_LIST):
	"""
	"""
	fastq_path_List = []
	fastq_path_List = sorted(glob.glob(DATA_DIR + "/" + "*.fastq*", recursive=True))
	fastq_path_List.extend(sorted(glob.glob(DATA_DIR + "/" + "*.fq*", recursive=True)))

	fastq_file_path_Dict = {}
	for each_fastq_path in fastq_path_List:
		##
		each_fastq_basename = os.path.basename(each_fastq_path)
		if LAYOUT != "single":
			#
			reverse_flag = False
			for each_delimiter in REVERSE_DELIMITER_LIST:
				##
				if each_delimiter in each_fastq_basename:
					#
					reverse_flag = True
				else:
					#
					pass
			else:
				##for each_delimiter in REVERSE_DELIMITER_LIST:
				pass
			if reverse_flag is True:
				#
				continue
			else:
				#if reverse_flag is True:
				pass

			for each_delimiter in FORWARD_DELIMITER_LIST:
				##
				if each_delimiter in each_fastq_basename:
					#
					each_fastq_basename = each_fastq_basename.replace(each_delimiter, "")
					break
				else:
					#
					pass
			else:
				##for each_delimiter in FORWARD_DELIMITER_LIST:
				pass
		else:
			#
			pass
		sample_name = "NO_SAMPLE_NAME"
		
		if SAMPLE_DELIMITER != "" and SAMPLE_DELIMITER in each_fastq_basename:
			#
			sample_name = each_fastq_basename.split(SAMPLE_DELIMITER, 1)[0]
		else:
			sample_name = each_fastq_basename.split(".fastq", 1)[0].split(".fq", 1)[0]
		#
		if sample_name not in fastq_file_path_Dict:
			#
			fastq_file_path_Dict[sample_name] = [each_fastq_path]
		else:
			fastq_file_path_Dict[sample_name].append(each_fastq_path)

	else:
		##for each_fastq_path in fastq_path_List:
		pass

	return fastq_file_path_Dict


def get_forward_fastq(wildcards):
	"""
	"""
	fastq_forward_List = fastq_forward_file_path_Dict[wildcards.sample]
	return fastq_forward_List
# ################################### PARAMTERS FUNCTIONS ########################## #


def build_paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc_command = ""
	else:
		#if LAYOUT == "single":
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc_command = "fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads " + str(PROCESSORS) + " $RESULT_PATH/{sample}.R2.trim_fastp.fastq.gz".format(sample=wildcards.sample)
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc_command += " 1>> $LOG_PATH/{sample}.trim_fastp.qc_fastqc.log 2>&1".format(sample=wildcards.sample)
	return paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc_command


def build_paired_trim_fastp_trim_fastp_concat_cat_qc_fastp_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastp_command = ""
	else:
		#if LAYOUT == "single":
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastp_command = "-I $RESULT_PATH/{sample}.R2.trim_fastp.fastq.gz".format(sample=wildcards.sample)

	return paired_trim_fastp_trim_fastp_concat_cat_qc_fastp_command
# ################################### CONFIGURATION ############################## #


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
# LAYOUT
FORWARD_DELIMITER_LIST = ["_1", "_R1"]
REVERSE_DELIMITER_LIST = ["_2", "_R2"]
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

fastq_forward_file_path_Dict = parse_datadir(DATA_DIR, SAMPLE_DELIMITER, LAYOUT, FORWARD_DELIMITER_LIST, REVERSE_DELIMITER_LIST)
#print(fastq_forward_file_path_Dict)

pre_process_List = []

for sample in fastq_forward_file_path_Dict:
	pre_process_List.append(WORK_DIR + "/{title}/{sample}/pre_process/{sample}.R1.trim_fastp.fastq.gz".format(title=TITLE, sample=sample))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		pre_process_List
# ################################### PIPELINE RULES ########################### #


rule trim_fastp:
	"""
	using fastp to trimp technical reads
	"""
	input:
		fastq_forward_path_List = get_forward_fastq,
	output:
		trim_fastq = WORK_DIR + "/{title}/{sample}/pre_process/{sample}.R1.trim_fastp.fastq.gz",
	threads: PROCESSORS
	message: "trim_fastp: {wildcards.sample}"
	resources:
		#mem_mb = MEMORY
	params:
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc = build_paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc_command,
		paired_trim_fastp_trim_fastp_concat_cat_qc_fastp = build_paired_trim_fastp_trim_fastp_concat_cat_qc_fastp_command,
	run:
		for each_fastq_forward_path in input.fastq_forward_path_List:
			##
			each_forward_fastq_name = os.path.basename(each_fastq_forward_path).split(".fastq", 1)[0].split(".fq", 1)[0]
			#
			if LAYOUT != "single":
				#
				each_reverse_fastq_name = each_forward_fastq_name
				each_unstranded_fastq_name = each_forward_fastq_name.split(FORWARD_DELIMITER_LIST[0], 1)[0].split(FORWARD_DELIMITER_LIST[1], 1)[0].split(REVERSE_DELIMITER_LIST[0], 1)[0].split(REVERSE_DELIMITER_LIST[1], 1)[0]
				for fwd, rev in zip(FORWARD_DELIMITER_LIST, REVERSE_DELIMITER_LIST):
					##
					each_reverse_fastq_name = each_reverse_fastq_name.replace(fwd, rev)
				else:
					##
					pass
				each_reverse_fastq_path = each_fastq_forward_path.replace(each_forward_fastq_name, each_reverse_fastq_name)
				#
				paired_trim_fastp_qc_fastqc = "fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads " + str(PROCESSORS) + " " + each_reverse_fastq_path
				paired_trim_fastp_qc_fastqc += " 1>> $LOG_PATH/" + each_unstranded_fastq_name + ".qc_fastqc.log 2>&1"
				#
				paired_trim_fastp_trim_fastp = "-I " + each_reverse_fastq_path + " -O $TEMP_PATH/" + each_reverse_fastq_name + ".trim_fastp.fastq"
				#
				paired_trim_fastp_trim_fastp_qc_fastqc = "fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads " + str(PROCESSORS) + " $TEMP_PATH/" + each_reverse_fastq_name + ".trim_fastp.fastq"
				paired_trim_fastp_trim_fastp_qc_fastqc += " 1>> $LOG_PATH/" + each_unstranded_fastq_name + ".trim_fastp.qc_fastqc.log 2>&1"
				#
				paired_trim_fastp_trim_fastp_concat_cat = "cat $TEMP_PATH/" + each_reverse_fastq_name + ".trim_fastp.fastq >> $RESULT_PATH/{sample}.R2.trim_fastp.fastq ".format(sample=wildcards.sample)
				paired_trim_fastp_trim_fastp_concat_cat += "2>> $LOG_PATH/" + each_unstranded_fastq_name + ".trim_fastp.concat_cat.log"
				
			else:
				#
				each_unstranded_fastq_name = each_forward_fastq_name
				paired_trim_fastp_qc_fastqc = ""
				paired_trim_fastp_trim_fastp = ""
				paired_trim_fastp_trim_fastp_qc_fastqc = ""
				paired_trim_fastp_trim_fastp_concat_cat = ""

			shell("""
				#
				##
				RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/pre_process
				mkdir -p $RESULT_PATH

				REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/pre_process
				mkdir -p $REPORT_PATH

				LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/pre_process
				mkdir -p $LOG_PATH
				
				TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{each_unstranded_fastq_name}/pre_process
				mkdir -p $TEMP_PATH
				
				##
				#


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample} | {each_unstranded_fastq_name}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "###################################- FASTQC -##################################" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads {threads} {each_fastq_forward_path} \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1" \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1

				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "%s\\n" "{paired_trim_fastp_qc_fastqc}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					#
				fi
				
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				#
				##
				start_time="$(date -u +%s)"
				
				fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads {threads} {each_fastq_forward_path} 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
				
				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					#
					##
					start_time="$(date -u +%s)"
					
					{paired_trim_fastp_qc_fastqc}

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.qc_fastqc.log 2>&1
					#
				fi


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample} | {each_unstranded_fastq_name}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "###################################- FASTP -##################################" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "fastp \\
				--thread {threads} \\
				-i {each_fastq_forward_path} -o $TEMP_PATH/{each_forward_fastq_name}.trim_fastp.fastq \\
				{paired_trim_fastp_trim_fastp} \\
				--cut_front \\
				--cut_tail \\
				--cut_mean_quality 30 \\
				--qualified_quality_phred 30 \\
				--unqualified_percent_limit 10 \\
				--length_required 50 \\
				--trim_poly_x \\
				--json $REPORT_PATH/{each_unstranded_fastq_name}.trim_fastp.json \\
				--html $REPORT_PATH/{each_unstranded_fastq_name}.trim_fastp.html \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1" \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				#
				##
				start_time="$(date -u +%s)"


				fastp \\
				--thread {threads} \\
				-i {each_fastq_forward_path} -o $TEMP_PATH/{each_forward_fastq_name}.trim_fastp.fastq \\
				{paired_trim_fastp_trim_fastp} \\
				--cut_front \\
				--cut_tail \\
				--cut_mean_quality 30 \\
				--qualified_quality_phred 30 \\
				--unqualified_percent_limit 10 \\
				--length_required 50 \\
				--trim_poly_x \\
				--json $REPORT_PATH/{each_unstranded_fastq_name}.trim_fastp.json \\
				--html $REPORT_PATH/{each_unstranded_fastq_name}.trim_fastp.html \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1


				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.log 2>&1


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample} | {each_unstranded_fastq_name}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "###################################- FASTQC -##################################" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads {threads} $TEMP_PATH/{each_forward_fastq_name}.trim_fastp.fastq \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1" \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1

				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "{paired_trim_fastp_trim_fastp_qc_fastqc}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					#
				fi

				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				#
				##
				start_time="$(date -u +%s)"

				fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads {threads} $TEMP_PATH/{each_forward_fastq_name}.trim_fastp.fastq \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				
				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1

				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					#
					##
					start_time="$(date -u +%s)"
					
					{paired_trim_fastp_trim_fastp_qc_fastqc}

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.qc_fastqc.log 2>&1
					#
				fi


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample} | {each_unstranded_fastq_name}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "###################################- CONCAT -##################################" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "cat $TEMP_PATH/{each_forward_fastq_name}.trim_fastp.fastq >> $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq \\
				2>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log" \\
				1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				
				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "%s\\n" "{paired_trim_fastp_trim_fastp_concat_cat}" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					#
				fi

				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				#
				##
				start_time="$(date -u +%s)"

				cat $TEMP_PATH/{each_forward_fastq_name}.trim_fastp.fastq >> $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq \\
				2>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
				
				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					#
					##
					start_time="$(date -u +%s)"

					{paired_trim_fastp_trim_fastp_concat_cat}

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{each_unstranded_fastq_name}.trim_fastp.concat_cat.log 2>&1
					#
				fi
			""")

		else:
			##for each_fastq_path in input.fastq_path_List:
			shell("""
				#
				##
				RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/pre_process
				mkdir -p $RESULT_PATH

				REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/pre_process
				mkdir -p $REPORT_PATH

				LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/pre_process
				mkdir -p $LOG_PATH

				TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/pre_process
				mkdir -p $TEMP_PATH
				##
				#


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample}" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "###################################- PIGZ -####################################" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "pigz --force $RESULT_PATH/*.fastq \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1" \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				#
				##
				start_time="$(date -u +%s)"

				pigz --force $RESULT_PATH/*.fastq 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.gzip_pigz.log 2>&1


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample}" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "###################################- FASTQC -##################################" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads {threads} $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq.gz \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1" \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1

				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "{params.paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc}" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					#
				fi

				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				#
				##
				start_time="$(date -u +%s)"


				fastqc -o $REPORT_PATH --dir $TEMP_PATH -f fastq --threads {threads} $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq.gz \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1

				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1

				if [ {LAYOUT} != 'single' ]; then
					#
					printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					#
					##
					start_time="$(date -u +%s)"

					{params.paired_trim_fastp_trim_fastp_concat_cat_qc_fastqc}

					end_time="$(date -u +%s)"
					##
					#
					printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastqc.log 2>&1
					#
				fi


				printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "trim_fastp: {wildcards.sample}" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "###################################- FASTP -##################################" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				
				printf "%s\\n" "fastp \\
				--thread {threads} \\
				-i $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq.gz \\
				{params.paired_trim_fastp_trim_fastp_concat_cat_qc_fastp} \\
				--json $REPORT_PATH/{wildcards.sample}.trim_fastp.qc_fastp.json \\
				--html $REPORT_PATH/{wildcards.sample}.trim_fastp.qc_fastp.html \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1" \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "##" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "#" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				#
				##
				start_time="$(date -u +%s)"

				fastp \\
				--thread {threads} \\
				-i $RESULT_PATH/{wildcards.sample}.R1.trim_fastp.fastq.gz \\
				{params.paired_trim_fastp_trim_fastp_concat_cat_qc_fastp} \\
				--json $REPORT_PATH/{wildcards.sample}.trim_fastp.qc_fastp.json \\
				--html $REPORT_PATH/{wildcards.sample}.trim_fastp.qc_fastp.html \\
				1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1


				end_time="$(date -u +%s)"
				##
				#
				printf "%s\\n" "" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1
				printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/{wildcards.sample}.trim_fastp.qc_fastp.log 2>&1


				if [ {MODE} != 'DEVELOPMENT' ]; then
					rm -rf $TEMP_PATH/*
				fi

			""")

# ################################### FINITO ################################### #