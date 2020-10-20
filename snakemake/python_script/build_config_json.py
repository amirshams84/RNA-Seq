# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Sep-20-2020
# Email: amir.shams84@gmail.com
# Project: RNA-Seq pipeline
# Aim: python script to process RNA-Seq_metadata
# ################################### IMPORT ##################################### #


import os
import sys
import pandas
import numpy
import json
import collections
# ################################### FUNCTIONS ################################## #


def is_file_exist(file_Path):
	"""
	"""
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		return False


def is_dataframe_empty(samplesheet_DF):
	"""
	"""
	if samplesheet_DF.empty is True or len(samplesheet_DF.index) == 0:
		#
		return True
	else:
		return False


def parse_reference_sheet_to_dict(reference_sheet_file_Path):
	"""
	"""
	if is_file_exist(reference_sheet_file_Path) is True:
		#
		pass
	else:
		#if is_file_exist(samplesheet_file_Path) is True:
		print("FATAL ERROR:Reference_sheet file path is not accessible or the premission to read is not granted!!!\nReference_sheet file Path: " + reference_sheet_file_Path + "\n")
		raise Exception("ABORTING!!!")
		sys.exit(2)
	#
	reference_sheet_DF = pandas.read_csv(
		reference_sheet_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter=",",
		index_col="ref_id"
	)
	#
	if is_dataframe_empty(reference_sheet_DF) is True:
		#dataframe is empty
		print("FATAL ERROR: reference_sheet is empty!!!\nReference_sheet file Path: " + reference_sheet_file_Path + "\n")
		raise Exception("ABORTING!!!")
		sys.exit(2)
	else:
		#
		pass
	#
	reference_sheet_DF.replace(numpy.nan, '', regex=True)
		
	#
	transposed_reference_sheet_DF = reference_sheet_DF.transpose()
	reference_sheet_Dict = transposed_reference_sheet_DF.to_dict()
	return reference_sheet_Dict


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
	samplesheet_Dict = {}
	#
	if sample_column is None:
		#set the default sample_id as index
		sample_column = "title"
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

	samplesheet_DF.replace(numpy.nan, '', regex=True)
		
	#
	transposed_samplesheet_DF = samplesheet_DF.transpose()
	samplesheet_Dict = transposed_samplesheet_DF.to_dict()
	#
	return samplesheet_Dict


def build_metadata_json(samplesheet_Dict, reference_sheet_Dict, json_directory_Path):
	"""
	"""
	for index, sample_Dict in samplesheet_Dict.items():
		##
		sample_metadata_Dict = collections.OrderedDict()

		input_path = sample_Dict["input"]
		output_path = sample_Dict["output"]
		reference_ID = sample_Dict["ref_id"].lower()
		sample_delimiter = str(sample_Dict["sample_delimiter"]).lower()
		platform = str(sample_Dict["platform"]).lower()
		if str(sample_Dict["single_end"]) == "1":
			#
			layout = "single"
		else:
			#
			layout = "paired"
		#
		sample_metadata_Dict["GENERAL"] = {"TITLE": sample_Dict["title"]}

		sample_metadata_Dict["DIRECTORY"] = {
			"DATA_DIR": input_path,
			"WORK_DIR": output_path,
			"TEMP_DIR": output_path + "/" + sample_Dict["title"] + "/temp"
		}

		sample_metadata_Dict["DATA"] = {
			"LAYOUT": layout,
			"PLATFORM": platform,
			"SAMPLE_DELIMITER": sample_delimiter
		}
		sample_metadata_Dict["REFERENCE"] = {}

		sample_metadata_Dict["REFERENCE"] = {
			reference_ID: {
				"KRAKEN2_DB": reference_sheet_Dict[reference_ID]["kraken2_index"],
				"GENE_GTF": reference_sheet_Dict[reference_ID]["genome_gtf"],
				"STAR_DB": reference_sheet_Dict[reference_ID]["star_index"],
				"KRONA_DB": reference_sheet_Dict[reference_ID]["krona_index"],
				"SORTMERNA_DB": reference_sheet_Dict[reference_ID]["sortmerna_index"],
				"STAR_FUSION_DB": reference_sheet_Dict[reference_ID]["star_fusion_index"],
				"GENE_BED": reference_sheet_Dict[reference_ID]["genome_bed"],
				"GENOME_FASTA": reference_sheet_Dict[reference_ID]["genome_fasta"],
				"EFFECTIVE_GENOME_SIZE": reference_sheet_Dict[reference_ID]["effective_genome_size"],
			}
		}
		
		with open(json_directory_Path + "/RNA-Seq_pipeline_" + sample_Dict["title"] + "_" + reference_ID + "_metadata.json", 'w') as metadata_json:
			json.dump(sample_metadata_Dict, metadata_json, indent=2)
	else:
		##for index, sample_Dict in samplesheet_Dict.items():
		pass

	return True
# ################################### EXEC ####################################### #

reference_sheet_Dict = parse_reference_sheet_to_dict(sys.argv[2])
#print(reference_sheet_Dict)
samplesheet_Dict = parse_samplesheet_to_dict(sys.argv[1])
build_metadata_json(samplesheet_Dict, reference_sheet_Dict, sys.argv[3])
# ################################### FINITO ##################################### #
