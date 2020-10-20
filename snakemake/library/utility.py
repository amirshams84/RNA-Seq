# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: July-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq
# Aim: python script for utility task
# ################################### IMPORT ##################################### #


import os
import sys
import pandas
from collections import OrderedDict
# ################################### FUNCTIONS ################################## #


def is_file_exist(file_Path):
	"""
	"""
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		return False


def is_path_readable(the_Path):
	"""
	"""
	if os.path.exists(the_Path) and os.access(the_Path, os.R_OK):
		#
		return True
	else:
		return False


def is_path_writable(the_Path):
	"""
	"""
	if os.path.exists(the_Path) and os.access(the_Path, os.R_OK) and os.access(the_Path, os.W_OK):
		#
		return True
	else:
		return False


def fix_path(the_Path):
	"""
	"""
	if the_Path[-1] == "/":
		#
		the_Path = the_Path[:-1]
	else:
		#
		pass
	return the_Path


def is_dataframe_empty(samplesheet_DF):
	"""
	"""
	if samplesheet_DF.empty is True or len(samplesheet_DF.index) == 0:
		#
		return True
	else:
		return False


def parse_samplesheet_to_dict(samplesheet_file_Path, sample_column=None):
	"""
	parse input samplesheet with all parameters into dict for future access
	"""
	if is_file_exist(samplesheet_file_Path) is True:
		#
		pass
	else:
		#if is_file_exist(samplesheet_file_Path) is True:
		print("FATAL ERROR:samplesheet file path is not accessible or the premission to read is not granted!!!\nSamplesheet file Path: " + samplesheet_file_Path + "\n")
		raise Exception("ABORTING!!!")
		sys.exit(2)
	#
	samplesheet_Dict = OrderedDict()
	#
	if sample_column is None:
		#set the default sample_id as index
		sample_column = "sample_id"
	else:
		#set the samplesheet index
		pass
	#
	samplesheet_DF = pandas.read_csv(
		samplesheet_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter=",",
		index_col=None
	)
	#
	if is_dataframe_empty(samplesheet_DF) is True:
		#dataframe is empty
		print("FATAL ERROR: Samplesheet file is empty!!!\nSamplesheet file Path: " + samplesheet_file_Path + "\n")
		raise Exception("ABORTING!!!")
		sys.exit(2)
	else:
		#
		pass
		
	#
	transposed_samplesheet_DF = samplesheet_DF.transpose()
	samplesheet_Dict = transposed_samplesheet_DF.to_dict()
	#
	return samplesheet_Dict


def build_sample_metadata_dict(samplesheet_Dict):
	"""
	read parse samplesheet dict into organized format
	"""
	sample_metadata_Dict = OrderedDict()
	for index, sample_Dict in samplesheet_Dict.items():
		##
		#print(index)
		sample_id = sample_Dict["sample_id"]
		#print(sample_id)
		if sample_id not in sample_metadata_Dict:
			#the first entery
			iterator = 0
			sample_metadata_Dict[sample_id] = OrderedDict()
			sample_metadata_Dict[sample_id][iterator] = {}
			sample_metadata_Dict[sample_id][iterator]["single_end"] = sample_Dict["single_end"]
			sample_metadata_Dict[sample_id][iterator]["is_sra"] = sample_Dict["is_sra"]
			sample_metadata_Dict[sample_id][iterator]["is_ftp"] = sample_Dict["is_ftp"]
			sample_metadata_Dict[sample_id][iterator]["host"] = sample_Dict["host"]
			sample_metadata_Dict[sample_id][iterator]["platform"] = sample_Dict["platform"]
			sample_metadata_Dict[sample_id][iterator]["fastq_1"] = sample_Dict["fastq_1"]
			sample_metadata_Dict[sample_id][iterator]["fastq_2"] = sample_Dict["fastq_2"]
		else:
			#if sample_id not in sample_metadata_Dict:
			#print(list(sample_metadata_Dict[sample_id].keys())[-1])
			iterator = list(sample_metadata_Dict[sample_id].keys())[-1] + 1
			sample_metadata_Dict[sample_id][iterator] = {}
			sample_metadata_Dict[sample_id][iterator]["single_end"] = sample_Dict["single_end"]
			sample_metadata_Dict[sample_id][iterator]["is_sra"] = sample_Dict["is_sra"]
			sample_metadata_Dict[sample_id][iterator]["is_ftp"] = sample_Dict["is_ftp"]
			sample_metadata_Dict[sample_id][iterator]["host"] = sample_Dict["host"]
			sample_metadata_Dict[sample_id][iterator]["platform"] = sample_Dict["platform"]
			sample_metadata_Dict[sample_id][iterator]["fastq_1"] = sample_Dict["fastq_1"]
			sample_metadata_Dict[sample_id][iterator]["fastq_2"] = sample_Dict["fastq_2"]
	else:
		##for sample, sample_Dict in samplesheet_Dict.items():
		pass
	return sample_metadata_Dict







