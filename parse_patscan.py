#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO	
import pandas as pd

__version__ = "0.1.0"
__author__ = "Pavlo Hrab"

def parse_args():
	parser = argparse.ArgumentParser(description="This script takes a genbank file and a patscan file and returns a csv file with the location of the patscan sequences in the genbank file",
	epilog=f"Version {__version__} - (C) {__author__}")
	parser.add_argument("-g", "--genbank_file", help="The genbank file")
	parser.add_argument("-p", "--patscan_file", help="The patscan file")
	parser.add_argument("-o", "--output_file", help="The output file")
	parser.add_argument("--write_faa", action="store_true", help="Write a fasta file with the protein sequences")
	parser.add_argument("--write_fna", action="store_true", help="Write a fasta file with the nucleotide sequences")
	parser.add_argument("-v","--version", action="version", version=f"{__version__}")
	return parser.parse_args()


def is_in_feature(genomic_coordinates, feature):
	start, end = genomic_coordinates
	cds_start, cds_end = feature.location.start, feature.location.end

	return cds_start <= start <= cds_end and cds_start <= end <= cds_end

def read_patscan(fl):
	result_dict = {}
	# Split the data into lines
	with open (fl, "r") as file:
		lines = file.readlines()
	# Process each line
	for index,line in enumerate(lines):
		line = line.strip()
		# Split the line into sequence and numbers
		seq_info, sequence = line.split(':')[1:]
		
		# Extract the sequence name
		seq_name = seq_info.split('[')[0].strip()
		
		# Extract the numbers from []
		numbers = [int(num) for num in seq_info.split('[')[1].split(']')[0].split(',')]
		
		# Add the information to the dictionary
		result_dict[f"Seq_{index}"] = {"sequence": sequence.strip(), "location": tuple(numbers)}
	return result_dict

def initialize_dict():
	data = {
			"ID": None,
			"Start": None,
			"Stop": None,
			"Sequence": None,
			"Intergenic": None,
			"In_feature": None,
			"Feature_start": None,
			"Feature_stop": None,
			"Feature_type": None,
			"Locus_tag": None,
			"Protein_id": None,
			"Feature_product": None,
			"Feature_translation": None,
			"Feature_sequence": None
		}
	return data

def update_dict(data, k,v, patscan_dic=None, feature=None, record=None, intergenic=False):
	if intergenic:
		data.update({
				"ID": k,
				"Start": patscan_dic[k]["location"][0],
				"Stop": patscan_dic[k]["location"][1],
				"Sequence": patscan_dic[k]["sequence"],
				"Intergenic" : True,
				"In_feature": False
			})
	else:
		data.update({
				"ID": k,
				"Start": v["location"][0],
				"Stop": v["location"][1],
				"Sequence": v["sequence"],
				"Intergenic" : False,
				"In_feature": True,
				"Feature_start": feature.location.start,
				"Feature_stop": feature.location.end,
				"Feature_type": feature.type,
				"Feature_sequence": feature.location.extract(record.seq)})
		if "locus_tag" in feature.qualifiers:
			data.update({"Locus_tag": feature.qualifiers["locus_tag"][0]})
		if "protein_id" in feature.qualifiers:
			data.update({"Protein_id": feature.qualifiers["protein_id"][0]})
		if "product" in feature.qualifiers:
			data.update({"Feature_product": feature.qualifiers["product"][0]})
		if "translation" in feature.qualifiers:
			data.update({"Feature_translation": feature.qualifiers["translation"][0]})
	return data

def add_intergenic(patscan_dic, intergenic, data_list):
	for k in intergenic:
		data = initialize_dict()
		data = update_dict(data, k,v=None, patscan_dic=patscan_dic, intergenic=True)
		data_list.append(data)
	return data_list

def get_sequence_location(records, patscan_dic):
	seen = []
	df = pd.DataFrame( columns=["ID", "Start", "Stop", "Sequence", "Intergenic", "In_feature", "Feature_start", "Feature_stop", "Feature_type", "Locus_tag","Protein_id", "Feature_product", "Feature_translation", "Feature_sequence"])
	data_list = []
	for record in records:
		for feature in record.features:
			for k,v in patscan_dic.items():
				data = initialize_dict()
				if is_in_feature(v["location"], feature):
					if feature.type != "source" and feature.type == "CDS":
						data = update_dict(data, k,v, feature=feature, record=record, intergenic=False)
						seen.append(k)
						data_list.append(data)
						continue
					elif feature.type == "gene":
						continue
					elif feature.type == "source":
						continue
					else:
						data = update_dict(data, k,v, feature=feature, record=record, intergenic=False)
						seen.append(k)
						data_list.append(data)
						continue
	# Get the intergenic sequences
	intergenic = [x for x in patscan_dic if x not in seen]
	data_list = add_intergenic(patscan_dic, intergenic, data_list)
	

	df = pd.DataFrame(data_list)
	return df

def write_faa_file(df, output_file):
	filtered_df = df[df['Feature_translation'].notna()]
	unique_seqs = filtered_df['Feature_translation'].unique().tolist()
	with open(f"{output_file}.faa", 'w') as fasta_file:
		for seq in unique_seqs:
			# Assuming you have a 'Sequence' column in your DataFrame
			locus_tag = filtered_df[filtered_df['Feature_translation'] == seq]['Locus_tag'].values[0]
			id_seq = filtered_df[filtered_df['Feature_translation'] == seq]['ID'].values[0]
			
			# Write to the fasta file
			fasta_file.write(f'>{locus_tag}__{id_seq}\n{seq}\n')

def write_fna_file(df, output_file):
	filtered_df = df[df['Feature_sequence'].notna()]
	unique_seqs = filtered_df['Feature_sequence'].unique().tolist()
	with open(f"{output_file}.fna", 'w') as fasta_file:
		for seq in unique_seqs:
			# Assuming you have a 'Sequence' column in your DataFrame
			locus_tag = filtered_df[filtered_df['Feature_sequence'] == seq]['Locus_tag'].values[0]
			id_seq = filtered_df[filtered_df['Feature_sequence'] == seq]['ID'].values[0]
			
			# Write to the fasta file
			fasta_file.write(f'>{locus_tag}__{id_seq}\n{seq}\n')

def transform_to_absolute_paths(args):
	args.genbank_file = os.path.abspath(args.genbank_file)
	args.patscan_file = os.path.abspath(args.patscan_file)
	args.output_file = os.path.abspath(args.output_file)
	return args

def main():
	args = parse_args()
	args = transform_to_absolute_paths(args)
	print(f"Working on {os.path.basename(args.genbank_file)} and {os.path.basename(args.patscan_file)}")
	records = SeqIO.parse(args.genbank_file, "genbank")
	patscan_dic = read_patscan(args.patscan_file)
	df = get_sequence_location(records, patscan_dic)
	if args.write_faa:
		write_faa_file(df, args.output_file)
	if args.write_fna:
		write_fna_file(df, args.output_file)
	df.to_csv(f"{args.output_file}.csv", index=False)

if __name__ == "__main__":
	main()