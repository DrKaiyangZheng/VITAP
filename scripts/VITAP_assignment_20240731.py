import csv
import argparse
import re
import os
import glob
import zipfile
import shutil
import subprocess
import statistics
import random
import numpy as np
import pandas as pd
import networkx as nx
from sys import argv, exit
import time
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor
from threading import Lock
from datetime import datetime
from shutil import copy
import threading
from datetime import datetime
from collections import defaultdict, Counter
from functools import partial
from multiprocessing import Pool
from itertools import combinations
from networkx.algorithms import bipartite
pd.options.mode.chained_assignment = None

def deplication_check(fasta_file):
	# Store Genome_IDs from the FASTA file
	genome_ids = []
	with open(fasta_file, 'r') as file:
		for line in file:
			if line.startswith('>'):
				# Extract content from '>' to the first whitespace
				genome_id = line[1:].split()[0]
				genome_ids.append(genome_id)
	# Find duplicate Genome_IDs
	duplicated_ids = [item for item, count in Counter(genome_ids).items() if count > 1]

	# If there are duplicate Genome_IDs
	if duplicated_ids:
		print("[WARNING] Duplicated Genome_IDs found in the input FASTA file:")
		for dup_id in duplicated_ids:
			print(f"-> {dup_id}")
		# Terminate the program
		exit(1)
	else:
		return


def extract_short_sequences(fasta_file, protein_file):
	fasta_ids = set([record.id for record in SeqIO.parse(fasta_file, "fasta")])
	
	protein_ids = set()
	for record in SeqIO.parse(protein_file, "fasta"):
		protein_id = record.id.rsplit('_', 1)[0]
		protein_ids.add(protein_id)
	
	short_ids = fasta_ids - protein_ids
	
	if not short_ids:
		return []
	
	short_sequences = [record for record in SeqIO.parse(fasta_file, "fasta") if record.id in short_ids]
	return short_sequences
		
def generate_short_gff(short_sequences, short_gff_file):
	with open(short_gff_file, "w") as gff:
		for record in short_sequences:
			gff.write(f"{record.id}\tSeqkit_translate\tCDS\t1\t{len(record)}\tNaN\tNaN\tNaN\tNaN\n")

def orf_count(gff_file):
	# reading GFF
	gff_df = pd.read_csv(gff_file, sep="\t", header=None, names=["id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
	# counting number of ORFs for each genome id
	orf_num_df = gff_df.groupby('id').size().reset_index(name='ORF_number')
	orf_num_df.set_index('id', inplace=True)
	return orf_num_df           

def taxonomy_assigning(blast_results_file, ictv_file, genome_cutoff_file, taxon_level, genome_length_file, orf_count_df, output_file, graph_output_file):
	# reading BLAST file
	blast_results = pd.read_csv(blast_results_file, sep="\t", names=["qseqid", "sseqid", f"bitscore"])

	# reading ICTV file
	ictv_data = pd.read_csv(ictv_file)
	
	# extract qseqid and sseqid of genome ids from BLAST result
	blast_results["qseqid_genome_id"] = blast_results["qseqid"].apply(lambda x: x.rsplit('_', 1)[0])
	blast_results["sseqid_genome_id"] = blast_results["sseqid"].apply(lambda x: x.split(".")[0])

	ictv_data.set_index("Virus GENBANK accession", inplace=True)

	blast_results = blast_results.join(ictv_data[taxon_level], on="sseqid_genome_id").rename(columns={taxon_level: "sseqid_taxa"})
	
	# Step 1: label self_hit
	blast_results['hit_type'] = blast_results.apply(lambda x: 'self_hit' if x['qseqid'] == x['sseqid'] else 'other', axis=1)

	# Step 2: find top_hit's sseqid_taxa of each qseqid
	# avoid to label self_hit as top_hit
	top_hits_taxa = blast_results[blast_results['hit_type'] == 'other'].groupby('qseqid')['sseqid_taxa'].first()

	# Step 3: Update the hit_type to top_hit for rows where the sseqid_taxa is the same as that of the top_hit
	blast_results['hit_type'] = blast_results.apply(lambda x: 'top_hit' if x['hit_type'] != 'self_hit' and x['sseqid_taxa'] == top_hits_taxa.get(x['qseqid'], None) else x['hit_type'], axis=1)

	# Step 4: Calculate bitscore_calculation_weight1: Based on classification by top hit
	blast_results['bitscore_calculation_weight1'] = blast_results.apply(lambda x: 1.2 if x['hit_type'] == 'top_hit' else (1 if x['hit_type'] == 'self_hit' else 0.8), axis=1)

	# Step 5: Calculate bitscore_calculation_weight2: Based on the most frequent classification
	qseqid_sseqid_taxa_counts = blast_results.groupby(['qseqid', 'sseqid_taxa']).size().reset_index(name='qseqid_sseqid_taxa_count')
	qseqid_total_counts = blast_results['qseqid'].value_counts().to_dict()
	qseqid_sseqid_taxa_counts['sseqid_taxa_percentage'] = qseqid_sseqid_taxa_counts.apply(lambda x: x['qseqid_sseqid_taxa_count'] / qseqid_total_counts[x['qseqid']], axis=1)
	qseqid_sseqid_taxa_counts['bitscore_calculation_weight2'] = qseqid_sseqid_taxa_counts['sseqid_taxa_percentage'].apply(lambda x: 1.2 if x > 0.5 else 1)
	blast_results = blast_results.merge(qseqid_sseqid_taxa_counts[['qseqid', 'sseqid_taxa', 'bitscore_calculation_weight2']], on=['qseqid', 'sseqid_taxa'])
	blast_results.loc[blast_results['hit_type'] == 'self_hit', 'bitscore_calculation_weight2'] = 1

	# calculate taxon_bitscore
	blast_results[f"{taxon_level}_bitscore"] = blast_results['bitscore'] * blast_results['bitscore_calculation_weight1'] * blast_results['bitscore_calculation_weight2']

	# Count the number of ORF hits for each genome in each taxonomic unit
	qseq_genome2taxon_df = blast_results[['qseqid', 'qseqid_genome_id', 'sseqid_taxa']].drop_duplicates()
	genome2taxon_df = qseq_genome2taxon_df.groupby(['qseqid_genome_id', 'sseqid_taxa']).size().reset_index(name='orf_occurance_taxon_count')

	# Extract the required columns
	genome_taxon2taxon_df = blast_results[["qseqid_genome_id", "sseqid_taxa", f"{taxon_level}_bitscore"]]
	genome_taxon2taxon_df = genome_taxon2taxon_df.merge(genome2taxon_df, on=['qseqid_genome_id', 'sseqid_taxa'])

	# Group by and perform sum and count
	grouped = genome_taxon2taxon_df.groupby(['qseqid_genome_id', 'sseqid_taxa', 'orf_occurance_taxon_count'])
	sum_df = grouped[f"{taxon_level}_bitscore"].sum().reset_index(name='sum_bitscore')
	
	# Calculate the occurrence frequency of each group
	count_df = grouped.size().reset_index(name='count')
	
	# Combine the results of sum and count
	grouped = pd.merge(sum_df, count_df, on=['qseqid_genome_id', 'sseqid_taxa', 'orf_occurance_taxon_count'])
	
	# Calculate the total sum of counts corresponding to each qseqid_taxa, compute the proportion, and add it to a new column
	grouped['total_count_per_qseqid_taxa'] = grouped.groupby('qseqid_genome_id')['count'].transform('sum')
	grouped['perc_sseqid_taxa'] = 10 * (grouped['count'] / grouped['total_count_per_qseqid_taxa'])
	
	genome_length_df = pd.read_csv(genome_length_file, sep="\t")
	genome_length_df.set_index('#id', inplace=True)
	grouped = grouped.merge(genome_length_df, left_on='qseqid_genome_id', right_index=True)
	grouped = grouped.merge(orf_count_df, left_on='qseqid_genome_id', right_index=True)

	# calculate taxon_score
	grouped = grouped.reset_index()  # Reset Index
	grouped[f'{taxon_level}_taxon_score'] = (grouped['sum_bitscore'] / grouped['count']) * ((grouped['orf_occurance_taxon_count'] / grouped['ORF_number']) ** 2) * grouped['perc_sseqid_taxa']
	
	# Extract the columns sseqid_taxa, taxon_score from the grouped data
	taxon2taxon_df = grouped[['qseqid_genome_id', 'sseqid_taxa', f'{taxon_level}_taxon_score']]

	# read thresholds of each Genome
	df_thresholds2 = pd.read_csv(genome_cutoff_file, sep="\t")
	
	# Set the threshold with taxon_level as the index
	df_thresholds2.set_index(taxon_level, inplace=True)

	# Merge the threshold dataframe with taxon2taxon_df based on the column sseqid_taxa
	taxon2taxon_df = taxon2taxon_df.join(df_thresholds2, on='sseqid_taxa')
	taxon2taxon_df[f'{taxon_level}_weight'] = taxon2taxon_df[f'{taxon_level}_taxon_score'] / taxon2taxon_df[f'{taxon_level}_score_cut-off']

	# keep necessary columns
	taxonomy_assigned_graph = taxon2taxon_df[['qseqid_genome_id', 'sseqid_taxa', f'{taxon_level}_weight']].rename(columns={'qseqid_genome_id': 'Genome_ID', 'sseqid_taxa': taxon_level})    
	taxonomy_assigned_graph.to_csv(output_file, sep="\t", index=False)


def multipartite_taxonomic_graph_generator(result_dir, vmr_mapping_file, taxon_categories):
	def process_row(row):
		max_average = -1
		max_count = 0
		max_lineage = None

		# If the weight is greater than or equal to 0.6, calculate the average weight from the current point to the end and the corresponding joint taxonomic level
		for i, weight in enumerate(weights):
			if row[weight] >= 0.6:
				average = row[weights[i:]].mean()
				lineage = ";".join(row[lineages[i:]])
				return pd.Series([average, lineage])

		# If all weights are less than 0.9, find the combination with the highest average weight
		total_max_average = -1
		for i in range(len(weights), 2, -1):  # Terminates when the number of remaining taxonomic level is 3
			current_weights = row[weights[-i:]]  # Recursing step by step from Species
			if not current_weights.isna().any():
				current_average = current_weights.mean()
				if current_average > total_max_average:
					total_max_average = current_average

		# Repeat the cycle to find combinations that meet the new criteria
		for i in range(len(weights), 2, -1):  # 修Terminates when the number of remaining taxonomic level is 3
			current_weights = row[weights[-i:]]
			if not current_weights.isna().any():
				current_average = current_weights.mean()
				current_count = i
				if current_average >= 0.95 * total_max_average and current_count >= 3:
					if current_count > max_count:
						max_average = current_average
						max_count = current_count
						max_lineage = ";".join(row[lineages[-i:]])

		if max_average != -1:
			return pd.Series([max_average, max_lineage])
		else:
			return pd.Series(['NaN', 'NaN'])

	def pad_lineage(lineage):
		elements = lineage.split(';')
		return ';'.join(['-'] * (8 - len(elements)) + elements)

	def count_dashes(lineage):
		return lineage.count('-')
		
	def replace_first_two_elements(lineage):
		elements = lineage.split(';')
		if elements[0] != '-' or elements[1] != '-':
			elements[0] = '-'
			elements[1] = '-'
		return ';'.join(elements)   

	print('[INFO] Importing taxonomic bipartite graphs and ICTV taxonomic hierarchy')
	species_annotation_file = os.path.join(result_dir, "Species_bipartite.graph")
	species_df = pd.read_csv(species_annotation_file, sep='\t')

	vmr_df = pd.read_csv(vmr_mapping_file)
	merged_df = species_df.merge(vmr_df[['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']], left_on='Species', right_on='Species', how='left')

	for taxon in taxon_categories[1:]:
		weight_column_name = f'{taxon}_weight'
		taxon_column_index = merged_df.columns.get_loc(taxon) + 1
		merged_df.insert(taxon_column_index, weight_column_name, None)

	for taxa in taxon_categories[1:]:
		unit_annotation_file = os.path.join(result_dir, f"{taxa}_bipartite.graph")
		unit_df = pd.read_csv(unit_annotation_file, sep='\t')
		merged_df = merged_df.merge(unit_df[['Genome_ID', taxa, f'{taxa}_weight']], on=['Genome_ID', taxa], how='left', suffixes=('', f'_from_{taxa}'))
		merged_df[f'{taxa}_weight'] = merged_df[f'{taxa}_weight_from_{taxa}']
		merged_df.drop(columns=[f'{taxa}_weight_from_{taxa}'], inplace=True)

	weights = ['Species_weight', 'Genus_weight', 'Family_weight', 'Order_weight', 'Class_weight', 'Phylum_weight', 'Kingdom_weight', 'Realm_weight']
	lineages = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']

	print('[INFO] Calculating average taxonomic score for target genomes: \n  -> finding the lowest effective taxonomic level\n  -> determining the optimal taxonomic hierarchy termination point')
	merged_df[['lineage_score', 'lineage']] = merged_df.apply(lambda row: process_row(row), axis=1)

	merged_df = merged_df[merged_df['lineage'].str.count(';') >= 3]

	#merged_df.to_csv('merged_df.csv', header=True, index=False)
	merged_df['lineage'] = merged_df['lineage'].apply(pad_lineage)
	all_lineage_df = merged_df[['Genome_ID', 'lineage', 'lineage_score']].drop_duplicates()

	print('[INFO] Determine best-fit taxonomic hierarchy')
	all_lineage_df['dash_count'] = all_lineage_df['lineage'].apply(count_dashes)
	high_score_df = all_lineage_df[all_lineage_df['lineage_score'] >= 0.9]
	best_high_score_indices = high_score_df.groupby('Genome_ID')[['dash_count', 'lineage_score']].idxmin()
	best_high_score_df = high_score_df.loc[best_high_score_indices['dash_count']]
	low_score_df = all_lineage_df[all_lineage_df['lineage_score'] < 0.9]
	low_score_df = low_score_df[~low_score_df['Genome_ID'].isin(best_high_score_df['Genome_ID'])]
	best_low_score_df = low_score_df.loc[low_score_df.groupby('Genome_ID')['lineage_score'].idxmax()]
	best_lineage_df = pd.concat([best_high_score_df, best_low_score_df]).drop_duplicates()
	best_lineage_df = best_lineage_df[['Genome_ID', 'lineage', 'lineage_score']]

	print('[INFO] Assigning confidence level to taxonomic hierarchy')
	conditions = [
		best_lineage_df['lineage_score'] >= 0.9,
		(best_lineage_df['lineage_score'] < 0.9) & (best_lineage_df['lineage_score'] >= 0.1),
	]
	choices = ['High-confidence', 'Medium-confidence']
	best_lineage_df['Confidence_level'] = np.select(conditions, choices, default='Low-confidence')
	best_lineage_df.loc[best_lineage_df['Confidence_level'].isin(['Low-confidence', 'Medium-confidence']), 'lineage'] = best_lineage_df.loc[best_lineage_df['Confidence_level'].isin(['Low-confidence', 'Medium-confidence']), 'lineage'].apply(replace_first_two_elements)
	#best_lineage_df.to_csv('best_lineage_df.csv', header=True, index=False)
	
	return best_lineage_df, all_lineage_df

def assignment():        
	# ======================================================================
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description = "The VITAP (Viral Taxonomic Assignment Pipeline, main program) v.1.5, Copyright 2024 Kaiyang Zheng, Viral/microbial diversity Lab. Ocean University of China", epilog='*******************************************************************\nExample usage: VITAP assignment -i target_genome.fasta -d ./DB_20230510/ -o VITAP_result\n*******************************************************************\n')
	args_parser.add_argument('-i', '--fasta', required=True, help = 'Input genome sequence file in FASTA format.')
	args_parser.add_argument('-d', '--db', required=True, help = 'VITAP database folder, which contains 10 required files: VMR genome mapping csv file, VMR protein DIAMOND database, eight lineage weight files (Species/Genus/Family/Order/Class/Phylum/Kingdom/Realm). ')
	args_parser.add_argument('-p', '--cpu', required=False, default = 1, help='Number of threads available for DIAMOND')
	args_parser.add_argument('-o', '--out', required=True, help = 'Result folder name')
	args_parser.add_argument("--low_conf", action = 'store_true', help="Exporting taxonomic assignments with low confidence, which will be discarded by default [lineage score < 0.1].")  
	args_parser = args_parser.parse_args()         
	# ======================================================================
	input_fasta = args_parser.fasta
	db_folder = args_parser.db + '/'
	threads = str(args_parser.cpu)
	result_folder = args_parser.out + '/'
	log_file = os.path.join(result_folder, "VITAP_run.log")
	os.makedirs(result_folder, exist_ok=True)
	taxon_categories = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
	start_time_wall_clock = time.time()
	start_time_cpu = time.process_time()

	print("===== The VITAP (Viral Taxonomic Assignment Pipeline, main program) v.1.5, Copyright 2024 Kaiyang Zheng, Viral/microbial diversity Lab. Ocean University of China ===== ")
	print("[INFO] A number of temporary files will be generated in your result directory, please keep them or delete them followed by the guidlines.")

	# ===== check duplicate ID =====
	deplication_check(input_fasta)

	# ===== ORFs prediction =====
	print(f'===== ORF calling for {input_fasta} =====') 
	merge_fasta = os.path.join(result_folder, "merge_genome.fasta")
	merge_gff = os.path.join(result_folder, "merge_genome.gff")
	input_faa = os.path.join(result_folder, "merge_genome.faa")
	select_genome = os.path.join(result_folder, "ICTV_selected_genomes.fasta")

	VMR_fasta_files = glob.glob(db_folder + '*.fasta')
	for fasta_file in VMR_fasta_files:
		VMR_fasta_file = fasta_file

	if not Path(input_faa).is_file():
		log_file = os.path.join(result_folder, "VITAP_run.log")
		print("[INFO] Ten ICTV reference genomes were randomly selected and merged to your genome set! ")
		with open(log_file, "a") as log_file:
			subprocess.run(["seqkit", "sample", "-n", "5", "-o", select_genome, VMR_fasta_file], stdout=log_file, stderr=log_file)
			with open(merge_fasta, 'w') as merged_file:
				with open(input_fasta, 'r') as ifile, open(select_genome, 'r') as sg:
					merged_file.write(ifile.read())
					merged_file.write(sg.read())
			subprocess.run(["prodigal-gv", "-p", "meta", "-i", merge_fasta, "-a", input_faa, "-f", "gff", "-o", merge_gff], stdout=log_file, stderr=log_file)
	else:
		print(f"[INFO] Using exited file : {input_faa}")           
	# extract short sequence
	short_sequences = extract_short_sequences(merge_fasta,input_faa)
	if short_sequences: 
		print(f"[INFO] Processing short sequences in {merge_fasta} ignored by prodigal.")  
		# save short sequence to file
		short_genome_file = os.path.join(result_folder, f"target_short_genome.fasta")
		SeqIO.write(short_sequences, short_genome_file, "fasta")
		# run seqkit command
		short_faa_file = os.path.join(result_folder, f"target_short_genome.faa")
		subprocess.run(["seqkit", "translate", "-f", "6", "-F", "--clean", "-o", short_faa_file, short_genome_file])
		# append reading frame of short sequence to file
		with open(input_faa, "a") as final_output_file, open(short_faa_file, "r") as short_output_file:
			final_output_file.write(short_output_file.read())
		os.remove(short_genome_file)
		os.remove(short_faa_file)   
	else:
		print("[INFO] Genome and ORF files have consistent non-redundant FASTA IDs. √")
		
	# ===== Calculate the length of sequences =====
	length_file = os.path.join(result_folder, f"merged_genome_length.tsv")
	if not Path(length_file).is_file():
		print(f"[INFO] Length statitic for merged_genome.fasta")
		with open("VITAP_run.log", "w") as log_file:
			subprocess.run(["seqkit", "fx2tab", "-l", "-n", "-i", "-H", "-o", length_file, merge_fasta])
	else:
		print(f"[INFO] The merge_genome_length.tsv exists, skipping.") 

	# ===== Calculate the number of ORFs =====
	print("[INFO] Statistic of the number of ORF per genome")
	orf_count_df = orf_count(merge_gff)

	# ===== align with VMR database using BLAST =====   
	print('===== Aligning to ICTV reference database based on BLAST-algorithm =====')

	# Initialize db_dmnd_file to None
	db_dmnd_file = None
	# Read CSV files
	VMR_csv_files = glob.glob(db_folder + '*.csv')
	for csv_file in VMR_csv_files:
		VMR_csv_file = csv_file

	# Read DMND files
	db_dmnd_files = glob.glob(db_folder + '*.dmnd')
	for dmnd_file in db_dmnd_files:
		db_dmnd_file = dmnd_file  # This will get the last dmnd file if there are multiple

	target_blast_fp = os.path.join(result_folder, f"target_ICTV_blastp.align")

	if not Path(target_blast_fp).is_file(): 
		print(f"[INFO] Please delete the output file ({target_blast_fp}) if this step aborted!")
		print(f"[INFO] Please keep the output file ({target_blast_fp}) if this step finished!")
		log_file = os.path.join(result_folder, "VITAP_run.log")
		with open(log_file, "a") as log_file:
			# Check if db_dmnd_file has been defined before running subprocess
			if db_dmnd_file:
				subprocess.run(["diamond", "blastp", "-p", threads, "-q", input_faa, "-d", db_dmnd_file, "--sensitive", "-f", "6", "qseqid", "sseqid", "bitscore", "-o", target_blast_fp, "-k", "1000", "--max-hsps", "1", "-e", "1e-3"], stdout=log_file, stderr=log_file)
			else:
				print("[WARNING] Error: Diamond database is not defined. Please check if dmnd files exist in the specified folder (Please do not forget to add '.' before relative path).")
	else:
		print(f"[INFO] Using exited file : {target_blast_fp}")

	# ===== calculate score =====
	print('===== Locating contig on taxonomic framework and generating genome community network =====') 
	for taxa in taxon_categories:   
		genome_cutoff_file = os.path.join(db_folder, f"{taxa}_genome.threshold")
		taxa_annot_file = os.path.join(result_folder, f"{taxa}_bipartite.graph")
		blast_graph_file = os.path.join(result_folder, "genome2genome.nwk")
		if not Path(taxa_annot_file).is_file():
			print(f"[INFO] Generating {taxa} graph for target genomes")
			taxonomy_assigning(target_blast_fp, VMR_csv_file, genome_cutoff_file, taxa, length_file, orf_count_df, taxa_annot_file, taxa_annot_file)
		else:
			print(f"[INFO] Using exited file: {taxa_annot_file} for {taxa} graph")
			
	# ===== complete the taxonomic classification and correction ===== 
	print("===== Determining viral lineages based on multi-partite graph ===== ")
	best_lineage_df, all_lineage_df = multipartite_taxonomic_graph_generator(result_folder, VMR_csv_file, taxon_categories)

	if not args_parser.low_conf:
		best_lineage_df.loc[best_lineage_df['Confidence_level'] == 'Low-confidence', 'lineage'] = '-;-;-;-;-;-;-;-'

	print("===== Exporting results ===== ")
	best_lineage_df_outfile = os.path.join(result_folder, f"best_determined_lineages.tsv")
	all_lineage_df_outfile = os.path.join(result_folder, f"all_lineages.tsv")

	best_lineage_df.to_csv(best_lineage_df_outfile, sep = '\t', index = False)
	all_lineage_df.to_csv(all_lineage_df_outfile, sep = '\t', index = False)
	
	end_time_wall_clock = time.time()
	end_time_cpu = time.process_time()
	elapsed_time_seconds_wall_clock = end_time_wall_clock - start_time_wall_clock
	elapsed_time_seconds_cpu = end_time_cpu - start_time_cpu
	elapsed_time_hours_wall_clock = elapsed_time_seconds_wall_clock / 3600
	elapsed_time_hours_cpu = elapsed_time_seconds_cpu / 3600
	print(f"[INFO] Time-consuming (wall clock): {elapsed_time_hours_wall_clock:.1f} hours")
	print(f"[INFO] Time-consuming (CPU): {elapsed_time_hours_cpu:.1f} hours")
	
	print('===== All done ===== ')
	
if __name__ == "__main__":
	assignment()    

