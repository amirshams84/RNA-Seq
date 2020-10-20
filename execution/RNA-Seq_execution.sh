#! /bin/bash
# ################################### INFO ####################################### #
# Author: Amir Shams
# Date: Sep-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq pipeline
# Aim: Bash script to build config json file and execute RNA-Seq pipeline
################################################################################## #
set -o pipefail
set -e
################################################################################## #
#
# Input Assertion
args=("$@")
if [ $# -eq 0 ]; then
    echo "No JSON file provided!!!"
    echo "Aborting!!"
    exit 1

fi
################################################################################## #
#
# Function to parse json file for basic/critical information
# remember the list is going to be like this ["element, "&&", "element", "&&", "element"]
declare -a json_List
function parse_json {
	
	index=0
	IFS="&&"
	export PYTHONIOENCODING=utf8
	for line in $(cat $1 | python -c 'import os,sys,json; \
		data = json.load(sys.stdin); \
		print(\
		os.path.expanduser(data["DIRECTORY"]["WORK_DIR"])+"&&"+\
		os.path.expanduser(data["DIRECTORY"]["TEMP_DIR"])+"&&"+\
		str(data["GENERAL"]["TITLE"])
		)')
	do
		#echo $index
		#echo $line
		json_List[$index]=$line
		index=$(($index+1))
		
	done
	IFS=" "
	
};
################################################################################## #
#
# Function which get the absoulte path of a file
function read_link() {
    local path=$1
    if [ -d $path ] ; then
        local abspath=$(cd $path; pwd -P)
    else
        local prefix=$(cd $(dirname -- $path) ; pwd -P)
        local suffix=$(basename $path)
        local abspath="$prefix/$suffix"
    fi
    if [ -e $abspath ] ; then
        echo $abspath
    else
        echo "$1 is not accessible"
		echo "make sure the path is correct(Please use absoulte path)!"
		echo "Aborting!!"
    	exit 1
    fi
};
################################################################################## #
#

SAMPLE_SHEET_FILE_PATH=${args[0]}
SAMPLE_SHEET_FILE_PATH_DIRECTORY=${args[1]}

#build json file path
EXEC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BUILD_JSON_SCRIPT_PATH=${EXEC_DIR}/../snakemake/python_script
REFERENCE_SHEET_FILE_PATH=${EXEC_DIR}/../template/RNA-Seq_reference_metadata.csv
#
python ${BUILD_JSON_SCRIPT_PATH}/build_config_json.py $SAMPLE_SHEET_FILE_PATH $REFERENCE_SHEET_FILE_PATH $SAMPLE_SHEET_FILE_PATH_DIRECTORY
################################################################################## #
#

for each_config in $(find $SAMPLE_SHEET_FILE_PATH_DIRECTORY -maxdepth 1 -mindepth 1 -name "[!.]*_metadata.json")
do
	##
	JSON_CONFIG_FILE=$each_config
	JSON_NAME=$(basename $JSON_CONFIG_FILE)
	JOB_NAME=${JSON_NAME/_metadata.json}
	#
	parse_json $JSON_CONFIG_FILE
	#
	WORK_DIR=${json_List[0]}
	WORK_DIR=$(read_link $WORK_DIR)
	mkdir -p $WORK_DIR
	#
	TITLE=${json_List[4]}
	#echo ${TITLE}
	mkdir -p ${WORK_DIR}/${TITLE}
	#
	TEMP_DIR=${json_List[2]}
	#echo ${TEMP_DIR}
	mkdir -p $TEMP_DIR
	#
	LOG_DIR=${WORK_DIR}/${TITLE}/cluster_log
	mkdir -p $LOG_DIR
	#
	mkdir -p $LOG_DIR/$JOB_NAME
	#
	SNAKEMAKE_DIR=${EXEC_DIR}/../snakemake
	SNAKE_FILE_DIR=$SNAKEMAKE_DIR/snakefile
	CLUSTER_CONFIG=${EXEC_DIR}/../snakemake/cluster_config
	####################################################
	#
	#EXECUTION_MODE=${args[1]}
	#EXECUTION_MODE='DEVELOPMENT'
	EXECUTION_MODE='PRODUCTION'

	if [ "$EXECUTION_MODE" == "DEVELOPMENT" ]
	then

		PROCESSORS=10
		MEMORY=10000

	elif [ "$EXECUTION_MODE" == "PRODUCTION" ]
	then

		CLUSTER_CONFIG_FILE=${CLUSTER_CONFIG}/cluster_production.yaml

	fi
	####################################################
	#
	EXECUTION_SCRIPT=${JSON_CONFIG_FILE/_metadata.json/_execution.sh}

	printf "%s\\n" "#!/bin/bash" > $EXECUTION_SCRIPT
	printf "%s\\n" "" >> $EXECUTION_SCRIPT
	printf "%s\\n" "source /data/${USER}/conda/etc/profile.d/conda.sh" >> $EXECUTION_SCRIPT
	printf "%s\\n" "" >> $EXECUTION_SCRIPT
	printf "%s\\n" "conda activate RNA-Seq_conda" >> $EXECUTION_SCRIPT
	printf "%s\\n" "" >> $EXECUTION_SCRIPT
	printf "%s%s\\n" "echo 'RNA-Seq Pipeline execution initiated at: '" "\$(date)" >> $EXECUTION_SCRIPT
	printf "%s\\n" "" >> $EXECUTION_SCRIPT
	printf "%s\\n" "cd $LOG_DIR/$JOB_NAME " >> $EXECUTION_SCRIPT

	if [ "$EXECUTION_MODE" = "DEVELOPMENT" ]
	then
		# #####################################################################################################################################################################
		#
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --cores --unlock" >> $EXECUTION_SCRIPT 

		printf "%s\\n" "" >> $EXECUTION_SCRIPT
		
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --cores" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT
		
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/mapping.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --cores" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT
		
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/mapping_qc.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --cores" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT

		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/reporting.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --cores" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT

	else
		# #####################################################################################################################################################################
		#
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --cores --unlock" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT
		
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --cores \
		--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster=\"sbatch --cpus-per-task={cluster.core} --mem={cluster.memory} \
		--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
		--output=${LOG_DIR}/${JOB_NAME}/{cluster.output} --error=${LOG_DIR}/${JOB_NAME}/{cluster.error} {cluster.extra}\"" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT
		
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/mapping.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --cores \
		--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster=\"sbatch --cpus-per-task={cluster.core} --mem={cluster.memory} \
		--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
		--output=${LOG_DIR}/${JOB_NAME}/{cluster.output} --error=${LOG_DIR}/${JOB_NAME}/{cluster.error} {cluster.extra}\"" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT
		
		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/mapping_qc.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --cores \
		--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster=\"sbatch --cpus-per-task={cluster.core} --mem={cluster.memory} \
		--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
		--output=${LOG_DIR}/${JOB_NAME}/{cluster.output} --error=${LOG_DIR}/${JOB_NAME}/{cluster.error} {cluster.extra}\"" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT


		printf "%s\\n" "snakemake --snakefile $SNAKE_FILE_DIR/reporting.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --cores \
		--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster=\"sbatch --cpus-per-task={cluster.core} --mem={cluster.memory} \
		--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
		--output=${LOG_DIR}/${JOB_NAME}/{cluster.output} --error=${LOG_DIR}/${JOB_NAME}/{cluster.error} {cluster.extra}\"" >> $EXECUTION_SCRIPT 
		
		printf "%s\\n" "" >> $EXECUTION_SCRIPT
	fi


	printf "%s%s\\n" "echo 'Pipeline execution successfully finished at: '" "\$(date)" >> $EXECUTION_SCRIPT
	printf "%s%s\\n" "# ################################### FINITO ##################################### #" >> $EXECUTION_SCRIPT

	####################################################
	#Submit job to slurm
	EXECUTION_LOG=${JSON_CONFIG_FILE/_metadata.json/_execution.log}
	printf "%s\\n" "$JOB_NAME Submitted to slurm with the following jobid:"
	sbatch --mem=10G --cpus-per-task=5 --partition=norm --time=2-00:00:00 $EXECUTION_SCRIPT

done
# ################################### FINITO ##################################### #