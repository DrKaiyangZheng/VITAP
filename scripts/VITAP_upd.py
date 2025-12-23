#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import math
import argparse
import re
import os
import glob
import zipfile
import shutil
import subprocess
import random
import polars as pl
import pandas as pd
from sys import argv, exit
from Bio import SeqIO
import pyrodigal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from datetime import datetime
from collections import defaultdict
from uniref90_accession2taxid import uniref90_accession2taxid

# =========================================================
# Utility functions (UNMODIFIED)
# =========================================================

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

# =========================================================
# Download genome
# =========================================================

def download_and_process_genome(row):
    global completed_downloads
    virus_id = row[0]

    if virus_id in downloaded_ids:
        with counter_lock:
            progress_bar.update(1)
        return

    output_file = os.path.join(output_folder, f"{virus_id}.fasta")

    for attempt in range(10):
        with open("VITAP_VMR_update.log", "a") as log_file:
            subprocess.run(
                ["efetch", "-id", virus_id, "-format", "fasta", "-db", "nuccore"],
                stdout=open(output_file, "w"),
                stderr=log_file
            )
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            break

    with counter_lock:
        completed_downloads += 1
        progress_bar.update(1)

# =========================================================
# FASTA helpers (UNMODIFIED)
# =========================================================

def remove_invalid_lines(file_path):
    valid_bases = set('ATCGRYKMSWBDHVN')
    with open(file_path, 'r') as file:
        lines = file.readlines()

    valid_lines = [
        line for line in lines
        if line.startswith('>') or all(c.upper() in valid_bases for c in line.strip())
    ]

    with open(file_path, 'w') as file:
        file.writelines(valid_lines)

def extract_short_sequences(fasta_file, protein_file):
    fasta_ids = {record.id for record in SeqIO.parse(fasta_file, "fasta")}
    protein_ids = {
        record.id.rsplit('_', 1)[0]
        for record in SeqIO.parse(protein_file, "fasta")
    }
    short_ids = fasta_ids - protein_ids
    if not short_ids:
        return []
    return [
        record for record in SeqIO.parse(fasta_file, "fasta")
        if record.id in short_ids
    ]

def generate_short_gff(short_sequences, short_gff_file):
    with open(short_gff_file, "w") as gff:
        for record in short_sequences:
            gff.write(
                f"{record.id}\tSeqkit_translate\tCDS\t1\t{len(record)}\tNaN\tNaN\tNaN\tNaN\n"
            )

def run_pyrodigal(fasta_in, faa_out, gff_out):
    """
    Run gene prediction using pyrodigal in metagenomic mode.
    Output:
      - Protein FASTA (same ID style as prodigal CLI)
      - GFF file compatible with VITAP downstream steps
    """

    finder = pyrodigal.GeneFinder(meta=True)

    protein_records = []
    gff_lines = ["##gff-version 3"]

    for record in SeqIO.parse(fasta_in, "fasta"):
        genome_id = record.id
        genes = finder.find_genes(str(record.seq))

        for idx, gene in enumerate(genes, start=1):
            # ---- Protein FASTA ----
            protein_id = f"{genome_id}_{idx}"
            protein_seq = gene.translate()

            protein_records.append(
                SeqRecord(
                    Seq(protein_seq),
                    id=protein_id,
                    description=""
                )
            )

            # ---- GFF ----
            start = gene.begin + 1        # pyrodigal is 0-based
            end = gene.end
            strand = "+" if gene.strand == 1 else "-"

            gff_lines.append(
                "\t".join([
                    genome_id,
                    "pyrodigal",
                    "CDS",
                    str(start),
                    str(end),
                    ".",
                    strand,
                    "0",
                    f"ID={protein_id}"
                ])
            )

    # Write protein FASTA
    SeqIO.write(protein_records, faa_out, "fasta")

    # Write GFF
    with open(gff_out, "w") as f:
        f.write("\n".join(gff_lines) + "\n")

# =========================================================
# ORF count (Polars)
# =========================================================

def orf_count(gff_file):
    """
    Robust ORF counting for GFF without using pl.read_csv (avoids segfault on malformed lines).
    Logic unchanged:
      Genome_id = id.split('.')[0]
      ORF_number = number of CDS/feature rows per Genome_id (here: number of non-comment, tabbed records)
    Returns: Polars DataFrame with columns ["Genome_id", "ORF_number"]
    """
    counts = {}

    with open(gff_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            # GFF is tab-separated with 9 columns; we only need column 1 (seqid/id)
            # Use split with maxsplit to avoid heavy work on long attributes
            parts = line.rstrip("\n").split("\t", 1)
            if not parts or not parts[0]:
                continue

            genome_id = parts[0].split(".", 1)[0]
            counts[genome_id] = counts.get(genome_id, 0) + 1

    # Return as Polars DF (keeps downstream join style consistent)
    return pl.DataFrame(
        {
            "Genome_id": list(counts.keys()),
            "ORF_number": list(counts.values()),
        }
    )

# =========================================================
# Taxon cutoff (POLARS failed "Segment Fault", possibly due to macOS x86_64 Rosetta/Intel )
# =========================================================
def _read_diamond_align_as_polars(blast_results_file: str) -> pl.DataFrame:
    """
    Robust reader for Diamond tabular output (qseqid sseqid bitscore).
    Avoids polars CSV parser segfaults by using pure-Python line parsing.
    Skips malformed lines safely.
    """
    qseqid = []
    sseqid = []
    bitscore = []

    with open(blast_results_file, "rb") as f:
        for raw in f:
            # Skip comments/empty
            if not raw or raw.startswith(b"#") or raw in (b"\n", b"\r\n"):
                continue
            # Drop NUL-containing lines (often poison for parsers)
            if b"\x00" in raw:
                continue
            try:
                line = raw.decode("utf-8", errors="replace").rstrip("\n\r")
            except Exception:
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            qs = parts[0].strip()
            ss = parts[1].strip()
            bs = parts[2].strip()
            if not qs or not ss or not bs:
                continue

            try:
                bsv = float(bs)
            except Exception:
                continue

            qseqid.append(qs)
            sseqid.append(ss)
            bitscore.append(bsv)

    return pl.DataFrame(
        {"qseqid": qseqid, "sseqid": sseqid, "bitscore": bitscore},
        schema={"qseqid": pl.Utf8, "sseqid": pl.Utf8, "bitscore": pl.Float64},
    )

def _read_ictv_map_as_polars(ictv_file: str, taxon_level: str) -> pl.DataFrame:
    """
    Robust reader for ICTV VMR csv. Only keeps:
      "Virus GENBANK accession" and taxon_level column
    Uses Python csv module-like parsing via polars? -> We'll do pure python for max safety.
    """
    accessions = []
    taxa = []

    with open(ictv_file, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        # Column names must match your file header
        acc_key = "Virus GENBANK accession"
        if acc_key not in reader.fieldnames:
            raise KeyError(f"ICTV file missing column: {acc_key}")
        if taxon_level not in reader.fieldnames:
            raise KeyError(f"ICTV file missing column: {taxon_level}")

        for row in reader:
            acc = (row.get(acc_key) or "").strip()
            tx = (row.get(taxon_level) or "").strip()
            if not acc:
                continue
            accessions.append(acc)
            taxa.append(tx if tx != "" else None)

    return pl.DataFrame(
        {"genome_id": accessions, "taxa": taxa},
        schema={"genome_id": pl.Utf8, "taxa": pl.Utf8},
    )

def _read_genome_length_as_polars(genome_length_file: str) -> pl.DataFrame:
    """
    Robust reader for seqkit fx2tab length TSV.
    We only need '#id' -> genome id (before '.') and keep all other columns for faithful merge.
    """

    rows = []
    with open(genome_length_file, "r", encoding="utf-8", errors="replace", newline="") as f:
        # It's TSV with header, produced by seqkit fx2tab -H
        reader = csv.DictReader(f, delimiter="\t")
        if "#id" not in reader.fieldnames:
            raise KeyError("Genome length file missing column: #id")

        for r in reader:
            # Skip empty
            if not r or not (r.get("#id") or "").strip():
                continue
            rid = r["#id"].strip()
            gid = rid.split(".", 1)[0]
            r["id"] = gid
            rows.append(r)

    if not rows:
        # empty df with at least id col
        return pl.DataFrame({"id": []}, schema={"id": pl.Utf8})

    # Let polars infer other columns; id forced to Utf8
    df = pl.DataFrame(rows)
    if "id" in df.columns:
        df = df.with_columns(pl.col("id").cast(pl.Utf8))
    return df

def taxon_cutoff(blast_results_file, ictv_file, taxon_level,
                 genome_length_file, orf_count_df, taxon_threshold_output):

    # ---- blast ----
    blast_results = pd.read_csv(
        blast_results_file,
        sep="\t",
        names=["qseqid", "sseqid", "bitscore"]
    )

    ictv_data = pd.read_csv(ictv_file)

    blast_results["qseqid_genome_id"] = blast_results["qseqid"].apply(lambda x: x.split(".")[0])
    blast_results["sseqid_genome_id"] = blast_results["sseqid"].apply(lambda x: x.split(".")[0])

    ictv_data.set_index("Virus GENBANK accession", inplace=True)

    blast_results = blast_results.join(
        ictv_data[taxon_level], on="qseqid_genome_id"
    ).rename(columns={taxon_level: "qseqid_taxa"})

    blast_results = blast_results.join(
        ictv_data[taxon_level], on="sseqid_genome_id"
    ).rename(columns={taxon_level: "sseqid_taxa"})

    # ---- Step 1: self_hit ----
    blast_results["hit_type"] = blast_results.apply(
        lambda x: "self_hit" if x["qseqid"] == x["sseqid"] else "other",
        axis=1
    )

    # ---- Step 2: top hit taxa ----
    top_hits_taxa = (
        blast_results[blast_results["hit_type"] == "other"]
        .groupby("qseqid")["sseqid_taxa"]
        .first()
    )

    # ---- Step 3: mark top_hit ----
    blast_results["hit_type"] = blast_results.apply(
        lambda x: "top_hit"
        if x["hit_type"] != "self_hit"
        and x["sseqid_taxa"] == top_hits_taxa.get(x["qseqid"], None)
        else x["hit_type"],
        axis=1,
    )

    # ---- Step 4: weight1 ----
    blast_results["bitscore_calculation_weight1"] = blast_results.apply(
        lambda x: 1.2
        if x["hit_type"] == "top_hit"
        else (1 if x["hit_type"] == "self_hit" else 0.8),
        axis=1,
    )

    # ---- Step 5: weight2 ----
    qseqid_sseqid_taxa_counts = (
        blast_results.groupby(["qseqid", "sseqid_taxa"])
        .size()
        .reset_index(name="qseqid_sseqid_taxa_count")
    )

    qseqid_total_counts = blast_results["qseqid"].value_counts().to_dict()

    qseqid_sseqid_taxa_counts["sseqid_taxa_percentage"] = (
        qseqid_sseqid_taxa_counts.apply(
            lambda x: x["qseqid_sseqid_taxa_count"] / qseqid_total_counts[x["qseqid"]],
            axis=1,
        )
    )

    qseqid_sseqid_taxa_counts["bitscore_calculation_weight2"] = (
        qseqid_sseqid_taxa_counts["sseqid_taxa_percentage"]
        .apply(lambda x: 1.2 if x > 0.5 else 1)
    )

    blast_results = blast_results.merge(
        qseqid_sseqid_taxa_counts[
            ["qseqid", "sseqid_taxa", "bitscore_calculation_weight2"]
        ],
        on=["qseqid", "sseqid_taxa"],
    )

    blast_results.loc[
        blast_results["hit_type"] == "self_hit",
        "bitscore_calculation_weight2",
    ] = 1

    # ---- taxon bitscore ----
    blast_results[f"{taxon_level}_bitscore"] = (
        blast_results["bitscore"]
        * blast_results["bitscore_calculation_weight1"]
        * blast_results["bitscore_calculation_weight2"]
    )

    # ---- ORF occurrence ----
    qseq_genome2taxon_df = blast_results[
        ["qseqid", "qseqid_genome_id", "sseqid_taxa"]
    ].drop_duplicates()

    genome2taxon_df = (
        qseq_genome2taxon_df
        .groupby(["qseqid_genome_id", "sseqid_taxa"])
        .size()
        .reset_index(name="orf_occurance_taxon_count")
    )

    genome_taxon2taxon_df = blast_results[
        ["qseqid_genome_id", "qseqid_taxa", "sseqid_taxa", f"{taxon_level}_bitscore"]
    ].merge(genome2taxon_df, on=["qseqid_genome_id", "sseqid_taxa"])

    grouped = genome_taxon2taxon_df.groupby(
        ["qseqid_genome_id", "qseqid_taxa", "sseqid_taxa", "orf_occurance_taxon_count"]
    )

    sum_df = grouped[f"{taxon_level}_bitscore"].sum().reset_index(name="sum_bitscore")
    count_df = grouped.size().reset_index(name="count")

    grouped = pd.merge(
        sum_df,
        count_df,
        on=["qseqid_genome_id", "qseqid_taxa", "sseqid_taxa", "orf_occurance_taxon_count"],
    )

    grouped["total_count_per_qseqid_taxa"] = grouped.groupby(
        "qseqid_genome_id"
    )["count"].transform("sum")

    grouped["perc_sseqid_taxa"] = 10 * (
        grouped["count"] / grouped["total_count_per_qseqid_taxa"]
    )

    genome_length_df = pd.read_csv(genome_length_file, sep="\t")
    genome_length_df["id"] = genome_length_df["#id"].apply(lambda x: x.split(".")[0])
    genome_length_df.set_index("id", inplace=True)

    grouped = grouped.merge(
        genome_length_df, left_on="qseqid_genome_id", right_index=True
    )

    grouped = grouped.merge(
        orf_count_df.to_pandas(),
        left_on="qseqid_genome_id",
        right_on="Genome_id",
    )

    grouped["taxon_score"] = (
        (grouped["sum_bitscore"] / grouped["count"])
        * ((grouped["orf_occurance_taxon_count"] / grouped["ORF_number"]) ** 2)
        * grouped["perc_sseqid_taxa"]
    )

    taxon2taxon_df = grouped[["qseqid_taxa", "sseqid_taxa", "taxon_score"]]

    taxon_score_thresholds = []
    for taxa in taxon2taxon_df["qseqid_taxa"].unique():
        same = taxon2taxon_df.loc[
            (taxon2taxon_df["qseqid_taxa"] == taxa)
            & (taxon2taxon_df["sseqid_taxa"] == taxa),
            "taxon_score",
        ]
        diff = taxon2taxon_df.loc[
            (taxon2taxon_df["qseqid_taxa"] == taxa)
            & (taxon2taxon_df["sseqid_taxa"] != taxa),
            "taxon_score",
        ]

        if not diff.empty:
            cutoff = (same.min() + diff.max()) / 2
        else:
            cutoff = same.min() * 3 / 4

        taxon_score_thresholds.append(
            {taxon_level: taxa, f"{taxon_level}_score_cut-off": cutoff}
        )

    pd.DataFrame(taxon_score_thresholds).to_csv(
        taxon_threshold_output, sep="\t", index=False
    )

def delete_temp_files(path):
    temp_files = Path(path).glob('*.temp')

    for temp_file in temp_files:
        temp_file.unlink()

#=======================================================================
def upd(args):
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

    # ===== Prepare UniRef90 database =====
    print("[INFO] Preparing UniRef90 database")

    uniref90_fasta_gz = os.path.join(updated_DB_folder, "uniref90.fasta.gz")
    uniref90_fasta = os.path.join(updated_DB_folder, "uniref90.fasta")
    taxdmp_zip = os.path.join(updated_DB_folder, "taxdmp.zip")
    taxdmp_dir = os.path.join(updated_DB_folder, "taxdmp")
    accession2taxid_file = os.path.join(updated_DB_folder, "uniref90.accession2taxid")
    uniref90_dmnd = os.path.join(updated_DB_folder, "uniref90.dmnd")

    # -- download UniRef90 fasta.gz (supports resume)
    subprocess.run(
        ["wget", "-c",
         "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz",
         "-O", uniref90_fasta_gz],
        check=True
    )

    # -- gunzip if fasta not exists
    if not Path(uniref90_fasta).is_file():
        subprocess.run(["gunzip", "-c", uniref90_fasta_gz],
                       stdout=open(uniref90_fasta, "w"),
                       check=True)

    # -- download taxonomy dump
    subprocess.run(
        ["wget", "-c",
         "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",
         "-O", taxdmp_zip],
        check=True
    )

    # -- unzip taxdmp
    if not Path(taxdmp_dir).is_dir():
        subprocess.run(["unzip", "-o", taxdmp_zip, "-d", taxdmp_dir], check=True)

    # -- generate accession2taxid
    if not Path(accession2taxid_file).is_file():
        uniref90_accession2taxid(uniref90_fasta, accession2taxid_file)
    else:
        print("[INFO] UniRef90 accession2taxid exists, skipping.")

    # -- build diamond database with taxonomy
    if not Path(uniref90_dmnd).is_file():
        subprocess.run(
            [
                "diamond", "makedb",
                "--in", uniref90_fasta,
                "--db", os.path.join(updated_DB_folder, "uniref90"),
                "--taxonmap", accession2taxid_file,
                "--taxonnodes", os.path.join(taxdmp_dir, "nodes.dmp"),
            ],
            check=True
        )
    else:
        print("[INFO] UniRef90 DIAMOND database exists, skipping.")

    # ===== Prodigal =====
    db_prot_file = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.faa")
    db_gff_file = os.path.join(updated_DB_folder, f"VMR_genome_{db_name}.gff")
    if not Path(db_prot_file).is_file() or not Path(db_gff_file).is_file():
        print(f"[INFO] ORF calling for VMR_genome_{db_name}.fasta")
        with open("VITAP_VMR_update.log", "w") as log_file:
            run_pyrodigal(db_genome_file, db_prot_file, db_gff_file)
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

