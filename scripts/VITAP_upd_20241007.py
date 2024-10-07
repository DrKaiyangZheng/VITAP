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
import pandas as pd
from sys import argv, exit
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
from collections import defaultdict
from functools import partial
from multiprocessing import Pool

def clean_virus_id(virus_id):
	if ':' in virus_id:
		virus_id = virus_id.split(':')[1].strip()
	return virus_id.strip()

def fill_empty_cells(row, header):
	filled_row = [row[0]]
	for idx, cell in enumerate(row[1:], 1):
		if cell != '':
			filled_row.append(cell)
		else:
			for next_cell in row[idx+1:]:
				if next_cell != '':
					filled_row.append(f"[{header[idx]}]_{next_cell}")
					break
			else:
				filled_row.append("")
	return filled_row

def extract_start_end_sites(virus_id):
	match = re.match(r'(\w+)\s*\((\d+)\.(\d+)\)', virus_id)
	if match:
		start_end = f"{match.group(2)}~{match.group(3)}"
		virus_id = match.group(1)
	else:
		start_end = "full_length"
	return virus_id, start_end
	
def download_and_process_genome(row):
	global completed_downloads
	virus_id = row[0]
	start_end_sites = row[-1]
	if virus_id in downloaded_ids:
		with counter_lock:
			progress_bar.update(1)
		return

	output_file = os.path.join(output_folder, f"{virus_id}.fasta")

	# download fasta using efetch
	for attempt in range(10):
		with open("VITAP_VMR_update.log", "a") as log_file:
			subprocess.run(["efetch", "-id", virus_id, "-format", "fasta", "-db", "nuccore"], stdout=open(output_file, "w"), stderr=log_file)

		if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
			break

		if attempt == 9:
			with open("VITAP_VMR_update.log", "a") as update_log:
				update_log.write(f"Fail to download genome with accession {virus_id}\n")

	# update progress bar
	with counter_lock:
		completed_downloads += 1
		progress_bar.update(1)

	return row 

def zip_folder(folder, output_zip):
	with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
		for root, _, files in os.walk(folder):
			for file in files:
				zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), folder))     

def user_input_func():
	global user_input
	while True:
		user_input = input("=====\n[WARNING] Enter 'D' to delete the VMR_Genome folder, 'R' to retain it, or 'C' to compress it into a zip file (D/R/C). The default operation 'C' will be automatically condect in 5 minutes: \n=====\n")
		user_input = user_input.lower()
		if user_input in ('d', 'r', 'c'):
			break
		else:
			print("[INFO] Invalid input. Please enter 'D', 'R', or 'C'.")

def delete_subfolders(folder):
	for item in os.listdir(folder):
		item_path = os.path.join(folder, item)
		if os.path.isdir(item_path):
			shutil.rmtree(item_path)

def remove_invalid_lines(file_path):
	valid_bases = set('ATCGRYKMSWBDHVN')
	with open(file_path, 'r') as file:
		lines = file.readlines()

	valid_lines = []
	for line in lines:
		if line.startswith('>') or all(c.upper() in valid_bases for c in line.strip()):
			valid_lines.append(line)

	with open(file_path, 'w') as file:
		file.writelines(valid_lines)

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
	# read GFF
	gff_df = pd.read_csv(gff_file, sep="\t", header=None, names=["id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
	# keep the substring before the last "." in the first column as a new column
	gff_df['Genome_id'] = gff_df['id'].apply(lambda x: x.split('.')[0])
	# count the number of ORFs for each Genome_id
	orf_num_df = gff_df.groupby('Genome_id').size().reset_index(name='ORF_number')
	orf_num_df.set_index('Genome_id', inplace=True)
	return orf_num_df

def taxon_cutoff(blast_results_file, ictv_file, taxon_level, genome_length_file, orf_count_df, taxon_threshold_output):  

	blast_results = pd.read_csv(blast_results_file, sep="\t", names=["qseqid", "sseqid", "bitscore"])

	ictv_data = pd.read_csv(ictv_file)

	blast_results["qseqid_genome_id"] = blast_results["qseqid"].apply(lambda x: x.split(".")[0])
	blast_results["sseqid_genome_id"] = blast_results["sseqid"].apply(lambda x: x.split(".")[0])

	ictv_data.set_index("Virus GENBANK accession", inplace=True)

	blast_results = blast_results.join(ictv_data[taxon_level], on="qseqid_genome_id").rename(columns={taxon_level: "qseqid_taxa"})
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
	
	# count the number of ORF hits for each genome in each taxonomic unit
	qseq_genome2taxon_df = blast_results[['qseqid', 'qseqid_genome_id', 'sseqid_taxa']].drop_duplicates()
	genome2taxon_df = qseq_genome2taxon_df.groupby(['qseqid_genome_id', 'sseqid_taxa']).size().reset_index(name='orf_occurance_taxon_count')

	# extract required columns
	genome_taxon2taxon_df = blast_results[["qseqid_genome_id", "qseqid_taxa", "sseqid_taxa", f"{taxon_level}_bitscore"]]
	del blast_results
	genome_taxon2taxon_df = genome_taxon2taxon_df.merge(genome2taxon_df, on=['qseqid_genome_id', 'sseqid_taxa'])

	# group by and perform sum and count
	grouped = genome_taxon2taxon_df.groupby(['qseqid_genome_id', 'qseqid_taxa', 'sseqid_taxa', 'orf_occurance_taxon_count'])
	sum_df = grouped[f"{taxon_level}_bitscore"].sum().reset_index(name='sum_bitscore')
	
	# calculate the occurrence frequency of each group
	count_df = grouped.size().reset_index(name='count')
	
	# combine the results of sum and count
	grouped = pd.merge(sum_df, count_df, on=['qseqid_genome_id', 'qseqid_taxa', 'sseqid_taxa', 'orf_occurance_taxon_count'])
	
	# calculate the total sum of counts corresponding to each qseqid_taxa, compute the proportion, and add it to a new column
	grouped['total_count_per_qseqid_taxa'] = grouped.groupby('qseqid_genome_id')['count'].transform('sum')
	grouped['perc_sseqid_taxa'] = 10 * (grouped['count'] / grouped['total_count_per_qseqid_taxa'])

	genome_length_df = pd.read_csv(genome_length_file, sep="\t")
	genome_length_df['id'] = genome_length_df['#id'].apply(lambda x: x.split('.')[0])
	genome_length_df.set_index('id', inplace=True)
	grouped = grouped.merge(genome_length_df, left_on='qseqid_genome_id', right_index=True)
	grouped = grouped.merge(orf_count_df, left_on='qseqid_genome_id', right_index=True)
	grouped = grouped.reset_index()  # reset index
	grouped['taxon_score'] = (grouped['sum_bitscore'] / grouped['count']) * ((grouped['orf_occurance_taxon_count'] / grouped['ORF_number']) ** 2) * grouped['perc_sseqid_taxa']

	taxon2taxon_df = grouped[['qseqid_taxa', 'sseqid_taxa', 'taxon_score']]
	
	# calculate taxon_score_cut-off
	taxon_score_thresholds = []
	for taxa in taxon2taxon_df['qseqid_taxa'].unique():
		same_taxa_scores = taxon2taxon_df.loc[(taxon2taxon_df['qseqid_taxa'] == taxa) & (taxon2taxon_df['sseqid_taxa'] == taxa), 'taxon_score']
		diff_taxa_scores = taxon2taxon_df.loc[(taxon2taxon_df['qseqid_taxa'] == taxa) & (taxon2taxon_df['sseqid_taxa'] != taxa), 'taxon_score']		
		if not diff_taxa_scores.empty:  # 如果存在$qseqid_taxa!=$sseqid_taxa的行
			taxon_score_cut_off = (min(same_taxa_scores) + max(diff_taxa_scores)) / 2
		else:
			taxon_score_cut_off = min(same_taxa_scores) * 3 / 4

		taxon_score_thresholds.append({f'{taxon_level}': taxa, f'{taxon_level}_score_cut-off': taxon_score_cut_off})

	taxon_score_thresholds_df = pd.DataFrame(taxon_score_thresholds)
	taxon_score_thresholds_df.to_csv(taxon_threshold_output, sep="\t", index=False)
		
def delete_temp_files(path):
	temp_files = Path(path).glob('*.temp')

	for temp_file in temp_files:
		temp_file.unlink()

#=======================================================================
args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description = "The VITAP (Viral Taxonomic Assignment Pipeline, update program) v.1.5, Copyright 2024 Kaiyang Zheng, Viral/microbial diversity Lab. Ocean University of China", epilog='*******************************************************************\nExample usage: VITAP upd --vmr raw_VMR.csv -o VMR_reformat.csv -d VMR-MSL\n*******************************************************************\n')
args_parser.add_argument('--vmr', required=True, help = 'Input raw ICTV VMR source file in csv format.')
args_parser.add_argument('-o', '--out', required=True, help = 'Reformat ICTV VMR source file for VITAP utilization')
args_parser.add_argument('-d', '--db', required=True, help = 'You need to assign a name for a new DB')
args_parser = args_parser.parse_args()         

# ===== Loading initial VMR table =====
input_file = args_parser.vmr
output_file = args_parser.out
output_folder = "VMR_Genome"
db_name = args_parser.db
os.makedirs(output_folder, exist_ok=True)
VMR_csv_file = output_file

#delete empty files first
for root, dirs, files in os.walk(output_folder):
	for name in files:
		file_path = os.path.join(root, name)
		if os.path.getsize(file_path) == 0:
			os.remove(file_path)

# ===== Reformatting VMR table =====
with open(input_file, "r", encoding="utf-8") as infile, open(output_file, "w", encoding="utf-8", newline='') as outfile:
	reader = csv.reader(infile)
	writer = csv.writer(outfile)

	# Writing table header
	header = next(reader)
	header.append('Start/End site')
	writer.writerow(header)

	# Processing data rows
	for row in reader:
		virus_ids = row[0].split(";")
		for virus_id in virus_ids:
			cleaned_id = clean_virus_id(virus_id)
			if cleaned_id:
				cleaned_id, start_end_sites = extract_start_end_sites(cleaned_id)
				new_row = [cleaned_id] + row[1:] + [start_end_sites]
				filled_row = fill_empty_cells(new_row, header)
				writer.writerow(filled_row)

# Y/N confirmation
while True:
	user_input = input(f"=====\n[WARNING] You need to double check the VMR dataframe generated in last step [{output_file}], and ensure all start/end information is correct. Please pay special attention to these categories: Peduoviridae (eg. AE006468), Belpaoviridae (eg. LK928904), all GTA-viriform (Bartogtaviriformidae and Rhodogtaviriformidae). The related metedata of these categories is incorrect in ICTV_VMR-MSL38_210426. Please make sure that all mistakes have been fixed and saved, and enter 'N' to abort, then rerun the program. If no mistakes existed yet, enter 'Y' to continue (Y/N): \n=====\n")
	if user_input.lower() == 'y':
		break
	elif user_input.lower() == 'n':
		print("[INFO] Exiting the program.")
		exit()

# ===== Downloading genome FASTA =====
VMR_csv_file = output_file

counter_lock = Lock()

completed_downloads = 0

downloaded_ids = {file[:-6] for file in os.listdir(output_folder) if file.endswith('.fasta')}      

with open(VMR_csv_file, "r", encoding="utf-8") as f:
	reader = csv.reader(f)

	next(reader)
	rows = list(reader)
	total_rows = len(rows)

	# efetch allows parallel downloads using up to 3 threads
	with ThreadPoolExecutor(max_workers=3) as executor:
		progress_bar = tqdm(total=total_rows, desc="Downloading genomes")
		futures = [executor.submit(download_and_process_genome, row) for row in rows]
		for future in as_completed(futures):
			future.result()
		progress_bar.close()

# Processing integrated viral sequences
for row in rows:
	virus_id = row[0]
	start_end_sites = row[-1]

	if start_end_sites == "full_length":
		continue

	start, end = map(int, start_end_sites.split('~'))
	start -= 1

	input_fasta = os.path.join(output_folder, f"{virus_id}.fasta")
	output_fasta = os.path.join(output_folder, f"{virus_id}_segment.fasta")

	if start is not None and end is not None:
		#os.system(f"seqkit subseq --quiet --id-regexp '^(\\S+)\.\s?' --chr {virus_id} -r {start}:{end} {input_fasta} | sed 's/>.* {virus_id}/>{virus_id}/g' > {output_fasta}")
		#os.system(f"seqkit subseq --quiet --id-regexp '^(\\\\S+)\\.\\s?' --chr {virus_id} -r {start}:{end} {input_fasta} | sed 's/>.* {virus_id}/>{virus_id}/g' > {output_fasta}")
		command = (
			f"seqkit subseq --quiet --id-regexp '^(\\\\S+)\\.\\s?' --chr {virus_id} "
			f"-r {start}:{end} {input_fasta} | "
			f"sed 's/>.* {virus_id}/>{virus_id}/g' > {output_fasta}"
		)

		with open("VITAP_VMR_update.log", "a") as log_file:
			subprocess.run(
				command,
				shell=True,
				stdout=log_file,
				stderr=subprocess.STDOUT
			)
		
		fai_path = os.path.join(output_folder, '*.fai')
		os.remove(input_fasta)
		for fai in glob.glob(fai_path):
			os.remove(fai)

print("[INFO] All files successfully downloaded and processed")

# ===== Get current date and generate new folder name =====
today = datetime.today().strftime('%Y%m%d')
#updated_DB_folder = f"DB_{today}"
updated_DB_folder = f"DB_{db_name}"
os.makedirs(updated_DB_folder, exist_ok=True)

# ===== Merge sequences into DB_(date) folder =====
db_genome_file = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.fasta")
if not Path(db_genome_file).is_file():
	print(f"[INFO] Merging file to VMR_genome_{db_name}.fasta...")
	with open(db_genome_file, "w") as db_genome:
		for fasta_file in glob.glob(os.path.join(output_folder, "*.fasta")):
			with open(fasta_file, "r") as single_fasta:
				db_genome.write(single_fasta.read())
		remove_invalid_lines(db_genome_file)        
else:
	print(f"[INFO] The VMR_genome_{db_name}.fasta exists, skipping.")
	
# ===== Statistical sequence length =====
length_file = os.path.join(updated_DB_folder, f"VMR_genome_length_{db_name}.tsv")
if not Path(length_file).is_file():
	print(f"[INFO] Length statitic for VMR_genome_{db_name}.fasta")
	with open("VITAP_VMR_update.log", "w") as log_file:
		subprocess.run(["seqkit", "fx2tab", "-l", "-n", "-i", "-H", "-o", length_file, db_genome_file])
else:
	print(f"[INFO] The VMR_genome_length_{db_name}.tsv exists, skipping.")			
	
# ===== Prodigal =====
db_prot_file = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.faa")
db_gff_file = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.gff")
if not Path(db_prot_file).is_file() or not Path(db_gff_file).is_file():
	print(f"[INFO] ORF calling for VMR_genome_{db_name}.fasta")
	with open("VITAP_VMR_update.log", "w") as log_file:
		subprocess.run(["prodigal", "-p", "meta", "-f", "gff", "-i", db_genome_file, "-a", db_prot_file, "-o", db_gff_file], stdout=log_file, stderr=log_file)
else:
	print(f"[INFO] VMR_genome_{db_name}.faa and VMR_genome_{db_name}.gff exist, skipping.")          
	

# ===== Short sequence extraction and end-to-end reading frame translation =====
print("[INFO] Processing short sequences ignored by prodigal.")
short_sequences = extract_short_sequences(db_genome_file, db_prot_file) 

if short_sequences:    
	short_genome_file = os.path.join(updated_DB_folder, f"VMR_short_genome_{db_name}.fasta")
	SeqIO.write(short_sequences, short_genome_file, "fasta")
	short_faa_file = os.path.join(updated_DB_folder, f"VMR_short_genome_{db_name}.faa")
	subprocess.run(["seqkit", "translate", "-f", "6", "-F", "--clean", "-o", short_faa_file, short_genome_file])
	with open(db_prot_file, "a") as final_output_file, open(short_faa_file, "r") as short_output_file:
		final_output_file.write(short_output_file.read())
	short_gff_file = os.path.join(updated_DB_folder, f"VMR_short_genome_{db_name}.gff")
	generate_short_gff(short_sequences, short_gff_file)
	with open(db_gff_file, "a") as final_gff_file, open(short_gff_file, "r") as short_gff_output_file:
		final_gff_file.write(short_gff_output_file.read())
	os.remove(short_genome_file)
	os.remove(short_faa_file)
	os.remove(short_gff_file)  
	del short_sequences 
else:
	print("[INFO] Genome and ORF files have consistent on non-redundant FASTA IDs. √")
	del short_sequences      

# ===== Cleaning GFF =====
print("[INFO] Cleaning the GFF file")

with open(db_gff_file, "r") as infile:
	lines = infile.readlines()
filtered_lines = [line for line in lines if not line.startswith('#')]

with open(db_gff_file, "w") as outfile:
	outfile.writelines(filtered_lines)
	
# ===== Statistical total number of ORFs =====
print("[INFO] Statistic of the number of ORF per genome")
db_gff_file = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.gff")
orf_count_df = orf_count(db_gff_file)	

# ===== Generate DB_VMR file =====
db_VMR_path = os.path.join(updated_DB_folder, f"VMR_taxonomy_map_{db_name}.csv")
print(f"[INFO] Moving {VMR_csv_file} to {updated_DB_folder} as {db_VMR_path}")
shutil.copy(VMR_csv_file, db_VMR_path)

# ===== self-Diamond  =====  
blast_fp = os.path.join(updated_DB_folder, f"Self_BLAST_{db_name}.align")
blast_db = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.dmnd")
if not Path(blast_db).is_file():    
	print("[INFO] Building ICTV reference protein databse.")
	subprocess.run(["diamond", "makedb", "--in", db_prot_file, "-d", blast_db])
	if not Path(blast_fp).is_file():		
		print("[INFO] Self-aligning of ICTV reference proteins.")
		with open("VITAP_VMR_update.log", "a") as log_file:
			subprocess.run(["diamond", "blastp", "-q", db_prot_file, "-d", blast_db, "-f", "6", "qseqid", "sseqid", "bitscore", "-o", blast_fp, "-k", "100", "--max-hsps", "1", "-e", "1e-3" ], stdout=log_file, stderr=log_file)
	else:
		print(f"[INFO] {blast_fp} exists, self-aligning was finished, skipping.")      
else:
	print("[INFO] ICTV reference protein databse exists.")
	if not Path(blast_fp).is_file():
		print("[INFO] Self-aligning of ICTV reference proteins.")
		log_file = open("VITAP_VMR_update.log", "w")
		subprocess.run(["diamond", "blastp", "-q", db_prot_file, "-d", blast_db, "-f", "6", "qseqid", "sseqid", "bitscore", "-o", blast_fp, "-k", "100", "--max-hsps", "1", "-e", "1e-3" ], stdout=log_file, stderr=log_file)
		log_file.close()
	else:
		print(f"[INFO] {blast_fp} exists, self-aligning was finished, skipping.")      

# ===== Assigning classification information to qseqid in Diamond alignment results =====
taxon_categories = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
for taxa in taxon_categories:
	taxon_threshold_output = os.path.join(updated_DB_folder, f'{taxa}_genome.threshold')	
	if not Path(taxon_threshold_output).is_file():
		print(f'[INFO] Calculating best-fit taxonomic threshold for {taxa}')	
		taxon_cutoff(blast_fp, VMR_csv_file, taxa, length_file, orf_count_df, taxon_threshold_output)	
	else:
		print(f"[INFO] {taxon_threshold_output} exists, skipping.")
						
delete_temp_files(updated_DB_folder)
print("[INFO] All updating steps finished.")

