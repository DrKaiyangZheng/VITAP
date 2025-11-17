import argparse
import os
import glob
import subprocess

import numpy as np
import pandas as pd
from sys import argv, exit
import time
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
from collections import Counter
import polars as pl
from typing import List, Tuple, Dict, Any, Optional

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

def orf_count(gff_file: str) -> pl.DataFrame:
    """
    Reads a GFF file and counts the number of entries (ORFs) for each genome ID.
    """
    # Define column names
    column_names = [
        "id", "source", "type", "start", "end", 
        "score", "strand", "phase", "attributes"
    ]

    # reading GFF
    gff_df = pl.read_csv(
        gff_file,
        comment_prefix="#",
        separator="\t",
        has_header=False,
        new_columns=column_names
    )

    # counting number of ORFs for each genome id
    # We group by 'id' and use pl.len() to get the size of each group.
    orf_num_df = gff_df.group_by("id").agg(
        pl.len().alias("ORF_number")
    )

    # Polars does not have a DataFrame index like pandas.
    # The 'id' column is already a key column in the result.
    # The .set_index() step is not needed.
    
    return orf_num_df         

def taxonomy_assigning(
    blast_results_file: str, 
    ictv_file: str, 
    genome_cutoff_file: str, 
    taxon_level: str, 
    genome_length_file: str, 
    gff_file: str,  # Assume this is a Polars DataFrame
    output_file: str,  # This parameter was unused in the original code
):
	
    orf_count_df = orf_count(gff_file)

    # Reading BLAST file
    blast_results = pl.read_csv(
        blast_results_file, 
        separator="\t", 
        has_header=False, 
        new_columns=["qseqid", "sseqid", "bitscore"]
    )

    # Reading ICTV file
    ictv_data = pl.read_csv(ictv_file)
    
    # Extract qseqid and sseqid of genome ids
    # .rsplit('_', 1)[0] -> Get everything before the last '_'
    # .split(".")[0] -> Get everything before the first '.'
    blast_results = blast_results.with_columns(
        pl.col("qseqid").str.extract(r"^(.*)_[^_]*$").alias("qseqid_genome_id"),
        pl.col("sseqid").str.split(".").list.first().alias("sseqid_genome_id")
    )

    # Polars joins on columns, not indexes.
    # Join with ictv_data to get 'sseqid_taxa'
    blast_results = blast_results.join(
        ictv_data.select(["Virus GENBANK accession", taxon_level]), 
        left_on="sseqid_genome_id", 
        right_on="Virus GENBANK accession"
    ).rename({taxon_level: "sseqid_taxa"})
    
    # Step 1: Label self_hit
    blast_results = blast_results.with_columns(
        pl.when(pl.col("qseqid") == pl.col("sseqid"))
          .then(pl.lit("self_hit"))
          .otherwise(pl.lit("other"))
          .alias("hit_type")
    )

    # Step 2: Find top_hit's sseqid_taxa of each qseqid
    top_hits_taxa_df = blast_results.filter(
        pl.col("hit_type") == "other"
    ).group_by("qseqid").agg(
        pl.first("sseqid_taxa").alias("top_hit_taxa")
    )

    # Join the top_hit_taxa back to the main dataframe
    blast_results = blast_results.join(
        top_hits_taxa_df, on="qseqid", how="left"
    )

    # Step 3: Update hit_type to top_hit
    blast_results = blast_results.with_columns(
        pl.when(
            (pl.col("hit_type") != "self_hit") & 
            (pl.col("sseqid_taxa") == pl.col("top_hit_taxa"))
        )
        .then(pl.lit("top_hit"))
        .otherwise(pl.col("hit_type"))
        .alias("hit_type")
    )

    # Step 4: Calculate bitscore_calculation_weight1
    blast_results = blast_results.with_columns(
        pl.when(pl.col("hit_type") == "top_hit").then(pl.lit(1.2))
          .when(pl.col("hit_type") == "self_hit").then(pl.lit(1.0))
          .otherwise(pl.lit(0.8))
          .alias("bitscore_calculation_weight1")
    )

    # Step 5: Calculate bitscore_calculation_weight2
    
    # Count occurrences of each (qseqid, sseqid_taxa) pair
    qseqid_sseqid_taxa_counts = blast_results.group_by(
        ["qseqid", "sseqid_taxa"]
    ).agg(
        pl.len().alias("qseqid_sseqid_taxa_count")
    )
    
    # Get total counts per qseqid
    qseqid_total_counts = blast_results.group_by("qseqid").agg(
        pl.len().alias("qseqid_total_count")
    )
    
    # Join to calculate percentage and weight2
    qseqid_sseqid_taxa_counts = qseqid_sseqid_taxa_counts.join(
        qseqid_total_counts, on="qseqid"
    ).with_columns(
        (pl.col("qseqid_sseqid_taxa_count") / pl.col("qseqid_total_count"))
        .alias("sseqid_taxa_percentage")
    ).with_columns(
        pl.when(pl.col("sseqid_taxa_percentage") > 0.5)
          .then(pl.lit(1.2))
          .otherwise(pl.lit(1.0))
          .alias("bitscore_calculation_weight2")
    )

    # Join weight2 back to blast_results
    blast_results = blast_results.join(
        qseqid_sseqid_taxa_counts.select(
            ["qseqid", "sseqid_taxa", "bitscore_calculation_weight2"]
        ),
        on=["qseqid", "sseqid_taxa"]
    )

    # Update weight2 for self_hit
    blast_results = blast_results.with_columns(
        pl.when(pl.col("hit_type") == "self_hit")
          .then(pl.lit(1.0))
          .otherwise(pl.col("bitscore_calculation_weight2"))
          .alias("bitscore_calculation_weight2")
    )

    # Calculate taxon_bitscore
    blast_results = blast_results.with_columns(
        (pl.col("bitscore") * pl.col("bitscore_calculation_weight1") * pl.col("bitscore_calculation_weight2"))
        .alias(f"{taxon_level}_bitscore")
    )

    # Count the number of ORF hits for each genome in each taxonomic unit
    qseq_genome2taxon_df = blast_results.select(
        ["qseqid", "qseqid_genome_id", "sseqid_taxa"]
    ).unique()
    
    genome2taxon_df = qseq_genome2taxon_df.group_by(
        ["qseqid_genome_id", "sseqid_taxa"]
    ).agg(
        pl.len().alias("orf_occurance_taxon_count")
    )

    # Extract required columns and merge
    genome_taxon2taxon_df = blast_results.select(
        ["qseqid_genome_id", "sseqid_taxa", f"{taxon_level}_bitscore"]
    ).join(
        genome2taxon_df, on=["qseqid_genome_id", "sseqid_taxa"]
    )

    # Group by and perform sum and count in one operation
    grouped = genome_taxon2taxon_df.group_by(
        ["qseqid_genome_id", "sseqid_taxa", "orf_occurance_taxon_count"]
    ).agg(
        pl.col(f"{taxon_level}_bitscore").sum().alias("sum_bitscore"),
        pl.len().alias("count")
    )

    # Calculate total_count_per_qseqid_taxa using a window function
    # and then calculate perc_sseqid_taxa
    grouped = grouped.with_columns(
        pl.col("count").sum().over("qseqid_genome_id").alias("total_count_per_qseqid_taxa")
    ).with_columns(
        (10 * (pl.col("count") / pl.col("total_count_per_qseqid_taxa"))).alias("perc_sseqid_taxa")
    )

    # Load genome_length and join
    genome_length_df = pl.read_csv(genome_length_file, separator="\t")
    
    # Join with genome_length and orf_count (which is already a Polars DF)
    grouped = grouped.join(
        genome_length_df, left_on="qseqid_genome_id", right_on="#id"
    ).join(
        orf_count_df, left_on="qseqid_genome_id", right_on="id"
    )

    # Calculate taxon_score
    grouped = grouped.with_columns(
        ((pl.col("sum_bitscore") / pl.col("count")) * ((pl.col("orf_occurance_taxon_count") / pl.col("ORF_number")) ** 2) * pl.col("perc_sseqid_taxa"))
        .alias(f"{taxon_level}_taxon_score")
    )

    # Extract the columns sseqid_taxa, taxon_score
    taxon2taxon_df = grouped.select(
        ["qseqid_genome_id", "sseqid_taxa", f"{taxon_level}_taxon_score"]
    )

    # Read thresholds and join
    df_thresholds2 = pl.read_csv(genome_cutoff_file, separator="\t")
    
    taxon2taxon_df = taxon2taxon_df.join(
        df_thresholds2, left_on="sseqid_taxa", right_on=taxon_level
    ).with_columns(
        (pl.col(f"{taxon_level}_taxon_score") / pl.col(f"{taxon_level}_score_cut-off"))
        .alias(f"{taxon_level}_weight")
    )

    # Keep necessary columns and write to file
    taxonomy_assigned_graph = taxon2taxon_df.select(
        pl.col("qseqid_genome_id").alias("Genome_ID"),
        pl.col("sseqid_taxa").alias(taxon_level),
        pl.col(f"{taxon_level}_weight")
    )
    
    taxonomy_assigned_graph.write_csv(output_file, separator="\t")

def _process_taxonomy_lists(
    w_list: List[Optional[float]], 
    l_list: List[Optional[str]]
) -> Dict[str, Any]:
    """
    This is the optimized logic from your `process_row` function.
    It operates on simple lists, making it much faster than using
    pandas Series objects for every row.
    """
    n_weights = len(w_list)

    # --- Logic 1: Check for >= 0.6 ---
    for i in range(n_weights):
        weight = w_list[i]
        if weight is not None and weight >= 0.6:
            # Calculate mean of the tail slice
            valid_weights = [w for w in w_list[i:] if w is not None]
            if not valid_weights:
                continue
            
            average = sum(valid_weights) / len(valid_weights)
            
            # Join the tail lineage
            valid_lineage = [str(e) for e in l_list[i:] if e is not None]
            lineage = ";".join(valid_lineage)
            
            return {"average": average, "lineage": lineage}

    # --- Logic 2: No weight >= 0.6 (Combined loops) ---
    tail_results = [] # Store (average, count)
    total_max_average = -1.0

    # Iterate from count = n_weights down to 3
    for i_count in range(n_weights, 2, -1):
        current_weights_list = w_list[-i_count:]
        
        # Check for NaNs (None)
        if None in current_weights_list:
            continue
            
        # We know they are all valid floats now
        current_average = sum(current_weights_list) / len(current_weights_list)
        tail_results.append((current_average, i_count))
        
        if current_average > total_max_average:
            total_max_average = current_average

    if total_max_average == -1.0:
        return {"average": None, "lineage": None} # No valid segment

    # Find best match from our saved results
    threshold = 0.95 * total_max_average
    
    # tail_results is already sorted by count descending (e.g., 7, 6, 5...)
    for avg, count in tail_results:
        if avg >= threshold:
            lineage_list = l_list[-count:]
            max_lineage = ";".join([str(e) for e in lineage_list if e is not None])
            return {"average": avg, "lineage": max_lineage}

    # No segment met the criteria
    return {"average": None, "lineage": None}


# --- Part 2: Main Polars Function ---

def multipartite_taxonomic_graph_generator_polars(
    result_dir: str, 
    vmr_mapping_file: str, 
    taxon_categories: List[str]
) -> Tuple[pl.DataFrame, pl.DataFrame]:

    print('[INFO] Importing taxonomic bipartite graphs and ICTV taxonomic hierarchy')
    species_annotation_file = os.path.join(result_dir, "Species_bipartite.graph")
    
    # Read and join the base species and VMR data
    species_df = pl.read_csv(species_annotation_file, separator='\t')
    vmr_df = pl.read_csv(vmr_mapping_file)
    
    merged_df = species_df.join(
        vmr_df.select(['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']),
        on='Species',
        how='left'
    )

    # Sequentially join all other taxon weight files
    # This is much cleaner than the pandas pre-allocation and merge loop
    for taxa in taxon_categories[1:]: # Skip 'Species', already processed
        unit_annotation_file = os.path.join(result_dir, f"{taxa}_bipartite.graph")
        
        # Check if file exists to avoid errors
        if not os.path.exists(unit_annotation_file):
            print(f"[WARN] File not found, skipping: {unit_annotation_file}")
            # Add a null column so downstream logic doesn't fail
            merged_df = merged_df.with_columns(pl.lit(None).alias(f'{taxa}_weight'))
            continue
            
        unit_df = pl.read_csv(
            unit_annotation_file, 
            separator='\t'
        ).select(
            ['Genome_ID', taxa, f'{taxa}_weight']
        )
        
        merged_df = merged_df.join(
            unit_df, 
            on=['Genome_ID', taxa], 
            how='left'
        )

    weights = ['Species_weight', 'Genus_weight', 'Family_weight', 'Order_weight', 'Class_weight', 'Phylum_weight', 'Kingdom_weight', 'Realm_weight']
    lineages = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']

    # Ensure all weight/lineage columns exist, filling with null if
    # any were missing (e.g., from skipped files)
    all_cols = weights + lineages + ["Genome_ID"] # Keep essential cols
    existing_cols = [col for col in all_cols if col in merged_df.columns]
    
    # Select existing columns and add nulls for any missing
    merged_df = merged_df.select(existing_cols)
    for col in weights + lineages:
        if col not in merged_df.columns:
            merged_df = merged_df.with_columns(pl.lit(None).alias(col))

    print('[INFO] Calculating average taxonomic score... (using fast Polars apply)')
    
    # Define the structure of the dictionary that _process_taxonomy_lists returns
    return_dtype = pl.Struct([
        pl.Field("average", pl.Float64),
        pl.Field("lineage", pl.String)
    ])
    
    # Pack columns into lists
    w_list_col = pl.list([pl.col(w) for w in weights])
    l_list_col = pl.list([pl.col(l) for l in lineages])

    # Apply the fast helper function using map_elements
    merged_df = merged_df.with_columns(
        pl.struct(w_list=w_list_col, l_list=l_list_col)
        .map_elements(
            lambda row: _process_taxonomy_lists(row["w_list"], row["l_list"]),
            return_dtype=return_dtype,
        )
        .alias("taxonomy_result")
    ).unnest("taxonomy_result").rename({"average": "lineage_score", "lineage": "lineage"})

    # Filter for valid lineages
    merged_df = merged_df.filter(
        pl.col('lineage').is_not_null() & (pl.col('lineage').str.count_matches(';') >= 3)
    )

    # Pad lineages to a standard length of 8
    # 'a;b;c' -> '-;-;-;-;-;a;b;c'
    merged_df = merged_df.with_columns(
        pl.col('lineage').str.split(';')
        .list.reverse()
        .list.pad_end(8, '-')
        .list.reverse()
        .str.join(';')
        .alias('lineage')
    )

    # Get unique lineages
    all_lineage_df = merged_df.select(
        ['Genome_ID', 'lineage', 'lineage_score']
    ).unique()

    print('[INFO] Determine best-fit taxonomic hierarchy')
    
    all_lineage_df = all_lineage_df.with_columns(
        pl.col('lineage').str.count_matches('-').alias('dash_count')
    )

    # --- High-score logic: min 'dash_count' (fewest dashes = most specific) ---
    high_score_df = all_lineage_df.filter(pl.col('lineage_score') >= 0.9)
    best_high_score_df = high_score_df.sort(
        "dash_count"
    ).group_by(
        'Genome_ID', maintain_order=False
    ).first()

    # --- Low-score logic: max 'lineage_score' ---
    low_score_df = all_lineage_df.filter(pl.col('lineage_score') < 0.9)
    
    # Filter out genomes that already have a high-confidence assignment
    low_score_df = low_score_df.join(
        best_high_score_df.select('Genome_ID'),
        on='Genome_ID',
        how='anti'
    ).filter(
        pl.col('lineage_score').is_not_null() # Drop nulls
    )

    best_low_score_df = low_score_df.sort(
        'lineage_score', descending=True
    ).group_by(
        'Genome_ID', maintain_order=False
    ).first()

    # Combine results
    best_lineage_df = pl.concat(
        [best_high_score_df, best_low_score_df]
    ).select(
        ['Genome_ID', 'lineage', 'lineage_score'] # Re-select to match pandas output
    )

    print('[INFO] Assigning confidence level to taxonomic hierarchy')
    
    # Assign confidence levels
    best_lineage_df = best_lineage_df.with_columns(
        pl.when(pl.col('lineage_score') >= 0.9).then(pl.lit('High-confidence'))
          .when(pl.col('lineage_score') >= 0.1).then(pl.lit('Medium-confidence'))
          .otherwise(pl.lit('Low-confidence'))
          .alias('Confidence_level')
    )

    # Replace first two elements for low/medium confidence
    best_lineage_df = best_lineage_df.with_columns(
        pl.when(pl.col('Confidence_level').is_in(['Low-confidence', 'Medium-confidence']))
        .then(
            # Replaces first two elements with '-': 'a;b;c;d' -> '-;-;c;d'
            pl.col('lineage').str.split(';').list.eval(
                pl.concat_list([
                    pl.lit(['-','-']), 
                    pl.element().list.slice(2)
                ])
            ).str.join(';')
        )
        .otherwise(pl.col('lineage'))
        .alias('lineage')
    )
    
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
			subprocess.run(["prodigal", "-p", "meta", "-i", merge_fasta, "-a", input_faa, "-f", "gff", "-o", merge_gff], stdout=log_file, stderr=log_file)
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
			taxonomy_assigning(target_blast_fp, VMR_csv_file, genome_cutoff_file, taxa, length_file, merge_gff, taxa_annot_file)
		else:
			print(f"[INFO] Using exited file: {taxa_annot_file} for {taxa} graph")
			
	# ===== complete the taxonomic classification and correction ===== 
	print("===== Determining viral lineages based on multi-partite graph ===== ")
	best_lineage_df, all_lineage_df = multipartite_taxonomic_graph_generator(result_folder, VMR_csv_file, taxon_categories)

	if args_parser.low_conf:
		def replace_first_four(lineage):
			if not isinstance(lineage, str):
				lineage = str(lineage) if not pd.isna(lineage) else ''
			elements = lineage.split(';')
			replaced_elements = ['-'] * 4
			replaced_elements += [str(e) if isinstance(e, str) else '-' for e in elements[4:8]]
			if len(replaced_elements) < 8:
				replaced_elements += ['-'] * (8 - len(replaced_elements))
			return ';'.join(replaced_elements[:8])
	
		best_lineage_df.loc[
			best_lineage_df['Confidence_level'] == 'Low-confidence', 
			'lineage'
		] = best_lineage_df.loc[
			best_lineage_df['Confidence_level'] == 'Low-confidence', 
			'lineage'
		].apply(replace_first_four)
	else:
		best_lineage_df = best_lineage_df[best_lineage_df['Confidence_level'] != 'Low-confidence']

	best_lineage_df['lineage'] = best_lineage_df['lineage'].astype(str)

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

