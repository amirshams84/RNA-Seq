#! /bin/bash
# ################################### INFO ####################################### #
# Author: Amir Shams
# Date: Sep-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq snakemake pipeline
# Aim: Bash script to run RNA-Seq pipeline
################################################################################## #
set -o pipefail
set -e

################################################################################## #
#
# Please edit this based on github repository cloned path

PIPELINE_DIRECTORY_PATH=/data/shamsaddinisha/Development/RNA-Seq
################################################################################## #
#
#PATH assertion
BUILD_JSON_PATH=${PIPELINE_DIRECTORY_PATH}/snakemake/python_script/build_config_json.py    
if [ ! -f $BUILD_JSON_PATH ]; then
	printf "%s%s\n" "PIPELINE_DIRECTORY_PATH:" $PIPELINE_DIRECTORY_PATH
	echo "PIPELINE_DIRECTORY_PATH not set properly, please make sute to enter the path properly."
	echo "ABORTING!!!!"
	exit 1
fi
################################################################################## #
#
#Input assertion
args=("$@")
if [ $# -eq 0 ]; then
	echo "No CSV metadata file provided, the pipeline requires sample metadata file in csv format."
	echo "ABORTING!!!!"
	exit 1
fi
################################################################################## #
#
#function which get the absoulte path of a file
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
SAMPLE_SHEET_FILE_PATH_DIRECTORY=$(read_link $SAMPLE_SHEET_FILE_PATH)
SAMPLE_SHEET_FILE_PATH_DIRECTORY=$(dirname $SAMPLE_SHEET_FILE_PATH_DIRECTORY)
#
RUN_DIRECTORY_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#
PIPELINE_EXECUTION_PATH=${PIPELINE_DIRECTORY_PATH}/execution/RNA-Seq_execution.sh

################################################################################## #
#
bash -c "$PIPELINE_EXECUTION_PATH $SAMPLE_SHEET_FILE_PATH $SAMPLE_SHEET_FILE_PATH_DIRECTORY"
# ################################### FINITO ##################################### #