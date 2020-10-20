# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Sep-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq
# Aim: Snakemake workflow for mapping
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


def parse_pre_process(WORK_DIR, TITLE):
	"""
	"""
	trim_fastq_path_List = []
	trim_fastq_path_Dict = {}
	trim_fastq_path_List = sorted(glob.glob(WORK_DIR + "/" + TITLE + "/**/pre_process/*.R1.trim_fastp.fastq.gz", recursive=True))
	for each_fastq_path in trim_fastq_path_List:
		##
		if utility.is_file_exist(each_fastq_path) is True:
			#
			each_sample = os.path.basename(each_fastq_path).split(".R1", 1)[0]
			trim_fastq_path_Dict[each_sample] = each_fastq_path
		else:
			#if utility.is_file_exist(each_fastq_path) is True:
			pass
	else:
		##for each_fastq_path in trim_fastq_path_List:
		pass
	return trim_fastq_path_Dict


def get_mapped_bam(wildcards):
	"""
	"""
	mapped_bam_List = []
	dehost_fastq_List = []
	for sample in pre_processed_fastq_path_Dict:
		##
		mapped_bam_List.append(WORK_DIR + "/" + TITLE + "/{sample}/mapping/{sample}.trimmed_fastp.mapped_bowtie2.bam".format(sample=sample))
		dehost_fastq_List.append(WORK_DIR + "/" + TITLE + "/{sample}/mapping/{sample}.R1.trimmed_fastp.dehost_kraken2.fastq.gz".format(sample=sample))
	else:
		##for sample, sample_Dict in sample_metadata_Dict.items():
		pass
	return mapped_bam_List


def get_dehost_fastq(wildcards):
	"""
	"""
	dehost_fastq_List = []
	for sample in pre_processed_fastq_path_Dict:
		##
		dehost_fastq_List.append(WORK_DIR + "/" + TITLE + "/{sample}/mapping/{sample}.R1.trimmed_fastp.dehost_kraken2.fastq.gz".format(sample=sample))
	else:
		##for sample, sample_Dict in sample_metadata_Dict.items():
		pass
	return dehost_fastq_List
# ################################### PARAMTERS FUNCTIONS ########################## #


def build_map_star_layout_star_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		map_star_layout_star_command = ""
	else:
		#
		map_star_layout_star_command = WORK_DIR + "/" + TITLE + "/{sample}/pre_process/{sample}.R2.trim_fastp.fastq.gz".format(sample=wildcards.sample)
		map_star_layout_star_command += " --outFilterIntronMotifs RemoveNoncanonical"

	return map_star_layout_star_command


def build_map_star_layout_samtools_filter_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		map_star_layout_samtools_filter_command = " -F 260 "
	else:
		#
		map_star_layout_samtools_filter_command = " -f 2 -F 256 "

	return map_star_layout_samtools_filter_command


def build_map_star_layout_star_unmapped_R1_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		map_star_layout_star_unmapped_R1_command = "mv $TEMP_PATH/{sample}.Unmapped.out.mate1 $RESULT_PATH/{sample}.R1.trim_fastp.unmap_star.fastq".format(sample=wildcards.sample)
	else:
		#
		map_star_layout_star_unmapped_R1_command = "mv $TEMP_PATH/{sample}.Unmapped.out.mate1 $RESULT_PATH/{sample}.R1.trim_fastp.unmap_star.fastq".format(sample=wildcards.sample)

	return map_star_layout_star_unmapped_R1_command


def build_map_star_layout_star_unmapped_R2_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		map_star_layout_star_unmapped_R2_command = ""
	else:
		#
		map_star_layout_star_unmapped_R2_command = "mv $TEMP_PATH/{sample}.Unmapped.out.mate2 $RESULT_PATH/{sample}.R2.trim_fastp.unmap_star.fastq".format(sample=wildcards.sample)

	return map_star_layout_star_unmapped_R2_command


def build_map_star_layout_star_samtools_fastq_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		map_star_layout_star_samtools_fastq_command = "samtools fastq " + WORK_DIR + "/" + TITLE + "/{sample}/map_star/{sample}.trim_fastp.map_star.bam ".format(sample=wildcards.sample)
		map_star_layout_star_samtools_fastq_command += "> $RESULT_PATH/{sample}.R1.trim_fastp.map_star.fastq".format(sample=wildcards.sample)
	else:
		#
		map_star_layout_star_samtools_fastq_command = "samtools fastq " + WORK_DIR + "/" + TITLE + "/{sample}/map_star/{sample}.trim_fastp.map_star.bam ".format(sample=wildcards.sample)
		map_star_layout_star_samtools_fastq_command += "-1 $RESULT_PATH/{sample}.R1.trim_fastp.map_star.fastq -2 $RESULT_PATH/{sample}.R2.trim_fastp.map_star.fastq".format(sample=wildcards.sample)

	return map_star_layout_star_samtools_fastq_command


def build_map_star_fusion_layout_fastq_command(wildcards):
	"""
	"""
	map_star_fusion_layout_fastq_command = ""
	if LAYOUT == "single":
		#
		map_star_fusion_layout_fastq_command = "--left_fq $RESULT_PATH/{sample}.R1.trim_fastp.map_star.fastq.gz ".format(sample=wildcards.sample)
	else:
		#
		map_star_fusion_layout_fastq_command = "--left_fq $RESULT_PATH/{sample}.R1.trim_fastp.map_star.fastq.gz --right_fq $RESULT_PATH/{sample}.R2.trim_fastp.map_star.fastq.gz ".format(sample=wildcards.sample)

	return map_star_fusion_layout_fastq_command


def build_map_star_layout_featureCounts_command(wildcards):
	"""
	"""
	if LAYOUT == "single":
		#
		map_star_layout_featureCounts_command = ""
	else:
		#
		map_star_layout_featureCounts_command = "-p "

	return map_star_layout_featureCounts_command

# ################################### CONFIGURATION ############################ #


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
EFFECTIVE_GENOME_SIZE = config_reference_Dict[HOST]["EFFECTIVE_GENOME_SIZE"]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
# CLUSTER PRAMETERS
#PROCESSORS = workflow.cores
PROCESSORS = max(2, cpu_count() - 4)
#MEMORY = 14096
MODE = "DEVELOPMENT"
# ------------------------------------
# ################################### WILDCARDS ################################ #

pre_process_fastq_path_Dict = parse_pre_process(WORK_DIR, TITLE)
#print(pre_processed_fastq_path_Dict)

pre_process_List = []
mapping_List = []

for sample in pre_process_fastq_path_Dict:
	##
	pre_process_List.append(WORK_DIR + "/{title}/{sample}/pre_process/{sample}.R1.trim_fastp.fastq.gz".format(title=TITLE, sample=sample))
	#
	mapping_List.append(WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam".format(title=TITLE, sample=sample))
	
else:
	##for sample in pre_process_fastq_path_Dict:
	pass

# ################################### PIPELINE FLOW ############################ #


rule End_Point:
	input:
		mapping_List
# ################################### PIPELINE RULES ########################### #

rule map_star:
	"""
	map_star
	"""
	input:
		trim_fastq = WORK_DIR + "/{title}/{sample}/pre_process/{sample}.R1.trim_fastp.fastq.gz",
	output:
		map_fastq = WORK_DIR + "/{title}/{sample}/map_star/{sample}.R1.trim_fastp.map_star.fastq.gz",
		map_bam = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.bam",
		map_count = WORK_DIR + "/{title}/{sample}/map_star/{sample}.trim_fastp.map_star.raw_count.tab",
	threads: PROCESSORS
	message: "map_star: {wildcards.sample}"
	#resources:
	#	mem_mb = MEMORY
	params:
		map_star_layout_star = build_map_star_layout_star_command,
		map_star_layout_samtools_filter = build_map_star_layout_samtools_filter_command,
		map_star_layout_star_unmapped_R1 = build_map_star_layout_star_unmapped_R1_command,
		map_star_layout_star_unmapped_R2 = build_map_star_layout_star_unmapped_R2_command,
		map_star_layout_star_samtools_fastq = build_map_star_layout_star_samtools_fastq_command,
		map_star_fusion_layout_fastq = build_map_star_fusion_layout_fastq_command,
		map_star_layout_featureCounts_command = build_map_star_layout_featureCounts_command,
	run:
		fastq_name = os.path.basename(input.trim_fastq).split(".fastq", 1)[0]
		sample_name = fastq_name.split(".R1", 1)[0]
		shell("""
			#
			##
			RESULT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/map_star
			mkdir -p $RESULT_PATH

			REPORT_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/report/map_star
			mkdir -p $REPORT_PATH

			LOG_PATH={WORK_DIR}/{TITLE}/{wildcards.sample}/log/map_star
			mkdir -p $LOG_PATH

			TEMP_PATH={TEMP_DIR}/{TITLE}/{wildcards.sample}/{fastq_name}/map_star
			mkdir -p $TEMP_PATH
			rm -rf $TEMP_PATH/*
			##
			#
			GENERAL_TAG={sample_name}.trim_fastp.map_star

			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "###################################- STAR -####################################" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "STAR \\
			--runThreadN {threads} \\
			--genomeDir {STAR_DB} \\
			--twopassMode Basic \\
			--readFilesIn {input.trim_fastq} {params.map_star_layout_star} \\
			--readFilesCommand zcat  \\
			--outSAMtype BAM SortedByCoordinate  \\
			--outTmpDir=$TEMP_PATH/{sample_name}  \\
			--outFileNamePrefix $TEMP_PATH/{sample_name}. \\
			--outSAMstrandField intronMotif \\
			--outSAMmode Full \\
			--outSAMattrRGline ID:{wildcards.sample} PL:illumina LB:{wildcards.sample} SM:{wildcards.sample} \\
			--chimSegmentMin 12  \\
			--chimJunctionOverhangMin 8  \\
			--chimOutJunctionFormat 1 \\
			--alignSJDBoverhangMin 10 \\
			--alignMatesGapMax 100000  \\
			--alignIntronMax 100000  \\
			--alignSJstitchMismatchNmax 5 -1 5 5  \\
			--chimMultimapScoreRange 3  \\
			--chimScoreJunctionNonGTAG -4  \\
			--chimMultimapNmax 20  \\
			--chimNonchimScoreDropMin 10  \\
			--outSAMattributes All  \\
			--limitBAMsortRAM 500000000 \\
			--quantMode GeneCounts  \\
			--outReadsUnmapped Fastx  \\
			--outSAMunmapped Within 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/{sample_name}.ReadsPerGene.out.tab {output.map_count}" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/{sample_name}.SJ.out.tab $RESULT_PATH/{sample_name}.SJ.out.tab" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/{sample_name}.Chimeric.out.junction $RESULT_PATH/${{GENERAL_TAG}}.chimeric_junction.bed" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "{params.map_star_layout_star_unmapped_R1}" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "{params.map_star_layout_star_unmapped_R2}" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/{sample_name}.Log.out $LOG_PATH/${{GENERAL_TAG}}.initial_report_star.log" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/{sample_name}.Log.progress.out $LOG_PATH/${{GENERAL_TAG}}.progress_report_star.log" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "cp $TEMP_PATH/{sample_name}.Log.final.out $LOG_PATH/${{GENERAL_TAG}}.final_report_star.log" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "mv $TEMP_PATH/{sample_name}.Log.final.out $REPORT_PATH/${{GENERAL_TAG}}.Log.final.out" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "cp {output.map_count} $REPORT_PATH/${{GENERAL_TAG}}.ReadsPerGene.out.tab" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "python {Python_Script_path}/STAR_SJtab2JunctionsBed.py -f $RESULT_PATH/${{GENERAL_TAG}}.SJ.out.tab > $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "LC_COLLATE=C sort -k1,1 -k2,2n $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed > $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed.sorted" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "bgzip -c $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed.sorted > $RESULT_PATH/${{GENERAL_TAG}}.splice_junction.bed.gz" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "tabix -f -p bed $RESULT_PATH/${{GENERAL_TAG}}.splice_junction.bed.gz" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			STAR \\
			--runThreadN {threads} \\
			--genomeDir {STAR_DB} \\
			--twopassMode Basic \\
			--readFilesIn {input.trim_fastq} {params.map_star_layout_star} \\
			--readFilesCommand zcat  \\
			--outSAMtype BAM SortedByCoordinate  \\
			--outTmpDir=$TEMP_PATH/{sample_name}  \\
			--outFileNamePrefix $TEMP_PATH/{sample_name}. \\
			--outSAMstrandField intronMotif \\
			--outSAMmode Full \\
			--outSAMattrRGline ID:{wildcards.sample} PL:illumina LB:{wildcards.sample} SM:{wildcards.sample} \\
			--chimSegmentMin 12  \\
			--chimJunctionOverhangMin 8  \\
			--chimOutJunctionFormat 1 \\
			--alignSJDBoverhangMin 10 \\
			--alignMatesGapMax 100000  \\
			--alignIntronMax 100000  \\
			--alignSJstitchMismatchNmax 5 -1 5 5  \\
			--chimMultimapScoreRange 3  \\
			--chimScoreJunctionNonGTAG -4  \\
			--chimMultimapNmax 20  \\
			--chimNonchimScoreDropMin 10  \\
			--outSAMattributes All  \\
			--limitBAMsortRAM 500000000 \\
			--quantMode GeneCounts  \\
			--outReadsUnmapped Fastx  \\
			--outSAMunmapped Within 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			mv $TEMP_PATH/{sample_name}.ReadsPerGene.out.tab {output.map_count}
			mv $TEMP_PATH/{sample_name}.SJ.out.tab $RESULT_PATH/${{GENERAL_TAG}}.SJ.out.tab
			mv $TEMP_PATH/{sample_name}.Chimeric.out.junction $RESULT_PATH/${{GENERAL_TAG}}.chimeric_junction.bed
			{params.map_star_layout_star_unmapped_R1}
			{params.map_star_layout_star_unmapped_R2}
			mv $TEMP_PATH/{sample_name}.Log.out $LOG_PATH/${{GENERAL_TAG}}.initial_report_star.log
			mv $TEMP_PATH/{sample_name}.Log.progress.out $LOG_PATH/${{GENERAL_TAG}}.progress_report_star.log
			cp $TEMP_PATH/{sample_name}.Log.final.out $LOG_PATH/${{GENERAL_TAG}}.final_report_star.log
			mv $TEMP_PATH/{sample_name}.Log.final.out $REPORT_PATH/${{GENERAL_TAG}}.Log.final.out
			cp {output.map_count} $REPORT_PATH/${{GENERAL_TAG}}.ReadsPerGene.out.tab

			python {Python_Script_path}/STAR_SJtab2JunctionsBed.py -f $RESULT_PATH/${{GENERAL_TAG}}.SJ.out.tab > $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed
			LC_COLLATE=C sort -k1,1 -k2,2n $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed > $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed.sorted
			bgzip -c $TEMP_PATH/${{GENERAL_TAG}}.SJ.out.bed.sorted > $RESULT_PATH/${{GENERAL_TAG}}.splice_junction.bed.gz
			tabix -f -p bed $RESULT_PATH/${{GENERAL_TAG}}.splice_junction.bed.gz

			
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.log 2>&1


			AWK_ALIGNMENT_FILTER="awk \'BEGIN{{OFS=FS}}{{if ( \$3 != \\\"chrUn\\\" && \$3 !~ /chrUn/ && \$3 !~ /random/ ) print \$0}}\'"

			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "###################################- FILTER SAMTOOLS -#########################" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "samtools view --threads {threads} -Sh {params.map_star_layout_samtools_filter} $TEMP_PATH/{sample_name}.Aligned.sortedByCoord.out.bam \\
			| $AWK_ALIGNMENT_FILTER \\
			| samtools sort --threads {threads} -O bam -T $TEMP_PATH/{sample_name} -o $TEMP_PATH/${{GENERAL_TAG}}.tmp - \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			samtools view --threads {threads} -Sh {params.map_star_layout_samtools_filter} $TEMP_PATH/{sample_name}.Aligned.sortedByCoord.out.bam \\
			| awk 'BEGIN{{OFS=FS}}{{if ( $3 != \"chrUn\" && $3 !~ /chrUn/ && $3 !~ /random/ ) print $0}}' \\
			| samtools sort --threads {threads} -O bam -T $TEMP_PATH/{sample_name} -o $TEMP_PATH/${{GENERAL_TAG}}.tmp - 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1

			
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "###################################- MARKDUPLICATE PICARD -####################" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "picard MarkDuplicates INPUT=$TEMP_PATH/${{GENERAL_TAG}}.tmp OUTPUT={output.map_bam} ASSUME_SORTED=true \\
			REMOVE_DUPLICATES=false METRICS_FILE=$REPORT_PATH/${{GENERAL_TAG}}.markduplicate_picard.txt VALIDATION_STRINGENCY=LENIENT \\
			TMP_DIR=$TEMP_PATH 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			picard MarkDuplicates INPUT=$TEMP_PATH/${{GENERAL_TAG}}.tmp OUTPUT={output.map_bam} ASSUME_SORTED=true \\
			REMOVE_DUPLICATES=false METRICS_FILE=$REPORT_PATH/${{GENERAL_TAG}}.markduplicate_picard.txt VALIDATION_STRINGENCY=LENIENT \\
			TMP_DIR=$TEMP_PATH 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "###################################- FEATURECOUNTS SUBREAD -###################" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "featureCounts -a {GENE_GTF} \\
			{params.map_star_layout_featureCounts_command} \\
			-T {threads} \\
			-o $RESULT_PATH/${{GENERAL_TAG}}.featureCounts \\
			{output.map_bam} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "mv $RESULT_PATH/${{GENERAL_TAG}}.featureCounts.summary $REPORT_PATH/${{GENERAL_TAG}}.featureCounts.summary" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
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
			{params.map_star_layout_featureCounts_command} \\
			-T {threads} \\
			-o $RESULT_PATH/${{GENERAL_TAG}}.featureCounts \\
			{output.map_bam} \\
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

			mv $RESULT_PATH/${{GENERAL_TAG}}.featureCounts.summary $REPORT_PATH/${{GENERAL_TAG}}.featureCounts.summary
			
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" 1>> $LOG_PATH/${{GENERAL_TAG}}.featureCounts_subread.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" 1> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "###################################- INDEX SAMTOOLS -####################" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "samtools index -@ {threads} -b {output.map_bam} 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			samtools index -@ {threads} -b {output.map_bam} 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1

			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "###################################- FASTQ SAMTOOLS -##########################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "{params.map_star_layout_star_samtools_fastq} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			{params.map_star_layout_star_samtools_fastq} \\
			2>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log

			
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "###################################- PIGZ -###############################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "pigz --force $RESULT_PATH/*.fastq 2>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			pigz --force $RESULT_PATH/*.fastq 2>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log
			
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.index_samtools.fastq_samtools.gzip_pigz.log 2>&1


			printf "%s\\n" "###################################- JOB INFO -################################" \\
			1> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "map_star: {wildcards.sample}" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "###################################- DEEPTOOLS BIGWIG -########################" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			
			printf "%s\\n" "bamCoverage --bam {output.map_bam} --outFileName $RESULT_PATH/${{GENERAL_TAG}}.bigwig --binSize 5 --normalizeUsing RPGC --extendReads 200 \\
			--outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {EFFECTIVE_GENOME_SIZE}" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1

			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "##" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "#" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "EXECUTING...." 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			#
			##
			start_time="$(date -u +%s)"

			bamCoverage --bam {output.map_bam} --outFileName $RESULT_PATH/${{GENERAL_TAG}}.bigwig --binSize 5 --normalizeUsing RPGC --extendReads 200 \\
			--outFileFormat bigwig --numberOfProcessors {threads} --effectiveGenomeSize {EFFECTIVE_GENOME_SIZE} \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			
			end_time="$(date -u +%s)"
			##
			#
			printf "%s\\n" "" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "DONE!!!!" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "ELAPSED TIME: %s seconds\\n" "$(($end_time-$start_time))" 1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1
			printf "%s\\n" "------------------------------------------------------------------------------" \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.filter_samtools.markduplicate_picard.deeptools_bigwig.log 2>&1

			#STAR FUSION
			printf "%s\\n" "#!/bin/bash" > $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "source ~/.bashrc" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "conda activate RNA-seq_conda" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "STAR-Fusion \\
			--chimeric_junction $RESULT_PATH/${{GENERAL_TAG}}.chimeric_junction.bed \\
			--genome_lib_dir {STAR_FUSION_DB} \\
			--FusionInspector inspect \\
			--examine_coding_effect \\
			{params.map_star_fusion_layout_fastq} \\
			--output_dir $RESULT_PATH/${{GENERAL_TAG}}.star_fusion \\
			--CPU 50 \\
			1>> $LOG_PATH/${{GENERAL_TAG}}.star_fusion.log 2>&1" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh
			printf "%s\\n" "echo 'Pipeline execution successfully finished at: '\$(date)" >> $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh

			
			#cd $LOG_PATH
			#sbatch --mem=100G --cpus-per-task=52 --partition=norm --time=1-00:00:00 $RESULT_PATH/${{GENERAL_TAG}}.star_fusion.sh \\
			#> $LOG_PATH/${{GENERAL_TAG}}.star_fusion.log 2>&1


			if [ {MODE} != 'DEVELOPMENT' ]; then
				rm -rf $TEMP_PATH/*
			fi


		""")

# ################################### FINITO ################################### #