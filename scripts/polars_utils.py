import os

import polars as pl
from typing import List, Tuple, Dict, Any, Optional

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