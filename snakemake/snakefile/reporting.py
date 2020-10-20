# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-08-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq
# Aim: Snakemake script to aggregate and report
# ################################### IMPORT ##################################### #


import os
import sys
import re
import glob
import xml.etree.ElementTree as ET
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


def get_bam():
	"""
	"""
	bam_List = []
	bam_List = glob.glob(WORK_DIR + "/" + TITLE + "/**/map_star/*.trim_fastp.map_star.bam", recursive=True)
	return bam_List


def get_bigwig():
	"""
	"""
	bigwig_List = []
	bigwig_List = glob.glob(WORK_DIR + "/" + TITLE + "/**/map_star/*.trim_fastp.map_star.bigwig", recursive=True)
	return bigwig_List
# ################################### PARAMTERS FUNCTIONS ########################## #


def build_featureCounts_layout_featureCounts_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		featureCounts_layout_featureCounts_command = ""
	else:
		#
		featureCounts_layout_featureCounts_command = "-p "

	return featureCounts_layout_featureCounts_command
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

report_List = []
report_List.append(WORK_DIR + "/{title}/report/{title}.igv_session.xml".format(title=TITLE))
report_List.append(WORK_DIR + "/{title}/report/{title}.pre_process_multiqc.html".format(title=TITLE))
report_List.append(WORK_DIR + "/{title}/report/{title}.featureCounts_subread.txt".format(title=TITLE))
# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		report_List
# ################################### PIPELINE RULES ########################### #


rule igv_session:
	"""
	"""
	input:
		bigwig_List = get_bigwig()
	output:
		igv_session = WORK_DIR + "/{title}/report/{title}.igv_session.xml",
	threads: PROCESSORS
	message: "igv_session: {wildcards.title}"
	#resources:
	#	mem_mb = MEMORY
	run:
		igv_session_String = "<?xml version='1.0' encoding='UTF-8' standalone='no'?>"
		igv_session_String += "<Session genome='hg19' hasGeneTrack='true' hasSequenceTrack='true' locus='All' version='8'>".format(igv_session=output.igv_session)
		Resource_String = "<Resources>"
		Panel_String = ""
		for each_bigwig in input.bigwig_List:
			##
			bigwig_name = os.path.basename(each_bigwig)
			sample_name = os.path.basename(bigwig_name).split(".bigwig", 1)[0]
			
			Resource_String += "<Resource path='/Volumes{bigwig_path}'/>".format(bigwig_path=each_bigwig)
			Panel_String += "<Panel height='1075' name='DataPanel' width='2543'>"
			Panel_String += "<Track attributeKey='{bigwig_name}' altColor='0,0,178' autoScale='true' clazz='org.broad.igv.track.DataSourceTrack' color='0,0,178' fontSize='15' id='/Volumes{bigwig_path}' name='{bigwig_name}' renderer='BAR_CHART' visible='true' windowFunction='mean'>".format(bigwig_path=each_bigwig, bigwig_name=bigwig_name)
			Panel_String += "<DataRange baseline='0.0' drawBaseline='true' flipAxis='false' maximum='' minimum='' type='LINEAR'/>"
			Panel_String += "</Track>"
			Panel_String += "</Panel>"
		else:
			##for each_bigwig in bigwig_List:
			pass
		Resource_String += "</Resources>"

		Panel_String += "<Panel height='277' name='FeaturePanel' width='1858'>"
		Panel_String += "<Track attributeKey='Reference sequence' clazz='org.broad.igv.track.SequenceTrack' fontSize='10' id='Reference sequence' name='Reference sequence' visible='true'/>"
		Panel_String += "<Track attributeKey='RefSeq Genes' clazz='org.broad.igv.track.FeatureTrack' color='0,0,178' colorScale='ContinuousColorScale;0.0;444.0;255,255,255;0,0,178' fontSize='10' height='35' id='hg19_genes' name='RefSeq Genes' visible='true'/>"
		Panel_String += "</Panel>"

		Panel_String += "<PanelLayout dividerFractions='0.789010989010989'/>"
		igv_session_String += Resource_String
		igv_session_String += Panel_String
		igv_session_String += "<HiddenAttributes>"
		igv_session_String += "<Attribute name='DATA FILE'/>"
		igv_session_String += "<Attribute name='DATA TYPE'/>"
		igv_session_String += "<Attribute name='NAME'/>"
		igv_session_String += "</HiddenAttributes>"
		igv_session_String += "</Session>"
		igv_session_Xml = ET.fromstring(igv_session_String)
		
		o = open(output.igv_session, 'wb')
		o.write(ET.tostring(igv_session_Xml, "utf-8"))
		#o.write(igv_session_String)
		o.close()
		# #################################################
		#
		shell("""
			sed -i "s/$USER\///g" {output.igv_session}
		""")


rule multiqc:
	"""
	"""
	output:
		pre_process_multiqc_report = WORK_DIR + "/{title}/report/{title}.pre_process_multiqc.html",
		mapping_multiqc_report = WORK_DIR + "/{title}/report/{title}.map_star_multiqc.html",
		mapping_qc_multiqc_report = WORK_DIR + "/{title}/report/{title}.map_star_qc_multiqc.html",
	threads: PROCESSORS
	message: "multiqc: {wildcards.title}"
	#resources:
	#	mem_mb = MEMORY
	run:
		shell("""

			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/report
			mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/report
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/report/log
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/report
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#
			GENERAL_TAG={TITLE}

			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "map_star: {wildcards.title}" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "###################################- PRE_PROCESS MULTIQC -####################################" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "multiqc --force \\
			--exclude general_stats \\
			--filename {WORK_DIR}/{wildcards.title}/report/${{GENERAL_TAG}}.pre_process_multiqc \\
			-d {WORK_DIR}/{wildcards.title}/**/report/pre_process \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			multiqc --force \\
			--exclude general_stats \\
			--filename {WORK_DIR}/{wildcards.title}/report/${{GENERAL_TAG}}.pre_process_multiqc \\
			-d {WORK_DIR}/{wildcards.title}/**/report/pre_process \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.pre_process_multiqc.log 2>&1

			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "map_star: {wildcards.title}" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "###################################- MAP_STAR MULTIQC -####################################" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "multiqc --force \\
			--exclude general_stats \\
			--filename {WORK_DIR}/{wildcards.title}/report/${{GENERAL_TAG}}.map_star_multiqc \\
			-d {WORK_DIR}/{wildcards.title}/**/report/map_star \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			multiqc --force \\
			--exclude general_stats \\
			--filename {WORK_DIR}/{wildcards.title}/report/${{GENERAL_TAG}}.map_star_multiqc \\
			-d {WORK_DIR}/{wildcards.title}/**/report/map_star \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_multiqc.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "map_star: {wildcards.title}" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "###################################- MAP_STAR_QC MULTIQC -####################################" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "multiqc --force \\
			--exclude general_stats \\
			--filename {WORK_DIR}/{wildcards.title}/report/${{GENERAL_TAG}}.map_star_qc_multiqc \\
			-d {WORK_DIR}/{wildcards.title}/**/report/map_star_qc \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			multiqc --force \\
			--exclude general_stats --ignore *bam_stat* --ignore *flagstat* \\
			--filename {WORK_DIR}/{wildcards.title}/report/${{GENERAL_TAG}}.map_star_qc_multiqc \\
			-d {WORK_DIR}/{wildcards.title}/**/report/map_star_qc \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.map_star_qc_multiqc.log 2>&1

		""")


rule featureCounts:
	"""
	"""
	input:
		bam_List = get_bam()
	output:
		featureCounts_table = WORK_DIR + "/{title}/report/{title}.featureCounts_subread.txt",
	threads: PROCESSORS
	message: "featureCounts: {wildcards.title}"
	#resources:
	#	mem_mb = MEMORY
	params:
		featureCounts_layout_featureCounts_command = build_featureCounts_layout_featureCounts_command,
	run:
		bam_String = " ".join(input.bam_List)
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/report
			mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/report
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/report/log
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/report
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#
			declare -a shell_bam_List=({input.bam_List})

			GENERAL_TAG={TITLE}
			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "map_star: {wildcards.title}" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "###################################- FEATURECOUNTS SUBREAD -####################################" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1

			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "featureCounts -a {GENE_GTF} \\
			{params.featureCounts_layout_featureCounts_command} \\
			-T {threads} \\
			-o $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread \\
			{bam_String} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1

			printf "%s\\n" "rm $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.summary" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "sed 1d $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread | cut -f1,7- > $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.tmp" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			for index in "${{!shell_bam_List[@]}}"
			do
				bam_path=${{shell_bam_List[$index]}}
				sample_name=$(basename $bam_path)
				printf "%s%s%s\\n" 'sed -i ' " s|$bam_path|$sample_name|g " '$REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.tmp' 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			done
			printf "%s\\n" "mv $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.tmp > {output.featureCounts_table}" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1

			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			featureCounts -a {GENE_GTF} \\
			{params.featureCounts_layout_featureCounts_command} \\
			-T {threads} \\
			-o $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread \\
			{bam_String} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			rm $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.summary
			

			sed 1d $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread | cut -f1,7- > $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.tmp

			for index in "${{!shell_bam_List[@]}}"
			do
				bam_path=${{shell_bam_List[$index]}}
				sample_name=$(basename $bam_path)
				sed -i "s|${{bam_path}}|${{sample_name}}|g" $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.tmp
			done

			mv $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread.tmp {output.featureCounts_table}
			rm $REPORT_PATH/${{GENERAL_TAG}}.featureCounts_subread
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1

		""")















