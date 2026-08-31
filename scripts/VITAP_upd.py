#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import math
import argparse
import re
import os
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


ACCESSION_COLUMN = "Virus GENBANK accession"
SEQUENCE_ID_COLUMN = "VITAP sequence ID"
START_END_COLUMN = "Start/End site"
FULL_LENGTH = "full_length"

# =========================================================
# Utility functions (UNMODIFIED)
# =========================================================

def clean_virus_id(virus_id):
    if ':' in virus_id:
        virus_id = virus_id.split(':', 1)[1].strip()
    return virus_id.strip()


def normalize_reference_id(reference_id):
    """Normalize a reference identifier for VITAP taxonomy-map lookups."""
    reference_id = str(reference_id).strip()
    if "__" in reference_id:
        return reference_id
    return reference_id.split(".", 1)[0]


def mapping_id_column(columns):
    """Prefer the range-aware ID while remaining compatible with older databases."""
    if SEQUENCE_ID_COLUMN in columns:
        return SEQUENCE_ID_COLUMN
    if ACCESSION_COLUMN in columns:
        return ACCESSION_COLUMN
    raise KeyError(
        f"ICTV mapping file requires {SEQUENCE_ID_COLUMN!r} or {ACCESSION_COLUMN!r}."
    )


def validate_file_identifier(identifier, field_name):
    """Reject identifiers that are unsafe or unsuitable as FASTA filenames."""
    identifier = str(identifier).strip()
    if not identifier or not re.fullmatch(r"[A-Za-z0-9_.-]+", identifier):
        raise ValueError(
            f"Invalid {field_name}: {identifier!r}. "
            "Only letters, numbers, underscores, dots, and hyphens are supported."
        )
    return identifier

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
    virus_id = virus_id.strip()
    match = re.fullmatch(r'([^()\s]+)\s*\((\d+)\.(\d+)\)', virus_id)
    if match:
        start = int(match.group(2))
        end = int(match.group(3))
        if start < 1 or end < start:
            raise ValueError(
                f"Invalid 1-based inclusive coordinates in {virus_id!r}: "
                f"start={start}, end={end}"
            )
        start_end = f"{start}~{end}"
        virus_id = match.group(1)
    else:
        if "(" in virus_id or ")" in virus_id:
            raise ValueError(f"Could not parse accession coordinates: {virus_id!r}")
        start_end = FULL_LENGTH
    validate_file_identifier(virus_id, ACCESSION_COLUMN)
    return virus_id, start_end


def make_sequence_id(accession, start_end_sites):
    """Create a stable, unique reference ID for a full record or coordinate slice."""
    accession = validate_file_identifier(accession, ACCESSION_COLUMN)
    base_accession = normalize_reference_id(accession)
    if start_end_sites == FULL_LENGTH:
        return base_accession
    start, end = map(int, start_end_sites.split("~"))
    if start < 1 or end < start:
        raise ValueError(
            f"Invalid 1-based inclusive coordinates: {start_end_sites!r}"
        )
    return f"{base_accession}__{start}_{end}"

# =========================================================
# Download genome
# =========================================================

def download_and_process_genome(
    virus_id,
    output_folder,
    downloaded_ids,
    progress_bar,
    counter_lock,
):
    virus_id = validate_file_identifier(virus_id, ACCESSION_COLUMN)
    try:
        if virus_id in downloaded_ids:
            return

        output_file = os.path.join(output_folder, f"{virus_id}.fasta")
        partial_file = f"{output_file}.part"
        success = False

        for _ in range(10):
            if os.path.exists(partial_file):
                os.remove(partial_file)
            with open("VITAP_VMR_update.log", "a") as log_file, open(partial_file, "w") as out:
                result = subprocess.run(
                    ["efetch", "-id", virus_id, "-format", "fasta", "-db", "nuccore"],
                    stdout=out,
                    stderr=log_file,
                    check=False,
                )

            if result.returncode != 0 or not os.path.exists(partial_file):
                continue
            if os.path.getsize(partial_file) == 0:
                continue

            try:
                next(SeqIO.parse(partial_file, "fasta"))
            except (StopIteration, ValueError):
                continue

            os.replace(partial_file, output_file)
            success = True
            break

        if os.path.exists(partial_file):
            os.remove(partial_file)
        if not success:
            raise RuntimeError(
                f"Failed to download a valid FASTA record for {virus_id} after 10 attempts. "
                "See VITAP_VMR_update.log for details."
            )
    finally:
        with counter_lock:
            progress_bar.update(1)


def load_accession_record(fasta_file, accession):
    """Load the FASTA record matching an accession, ignoring a version suffix."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise RuntimeError(f"No FASTA record found in {fasta_file}")

    accession_base = normalize_reference_id(accession)
    matches = [
        record for record in records
        if normalize_reference_id(record.id) == accession_base
    ]
    if len(matches) == 1:
        return matches[0]
    if len(records) == 1:
        return records[0]
    raise RuntimeError(
        f"Could not uniquely match accession {accession!r} in {fasta_file}; "
        f"found {len(records)} FASTA records."
    )


def write_segment_fasta(record, start, end, sequence_id, output_file):
    """Write one 1-based inclusive interval atomically with a unique FASTA ID."""
    sequence_id = validate_file_identifier(sequence_id, SEQUENCE_ID_COLUMN)
    if start < 1 or end < start or end > len(record):
        raise ValueError(
            f"Coordinates {start}~{end} are outside record {record.id!r} "
            f"(length {len(record)})."
        )

    expected_length = end - start + 1
    segment = SeqRecord(
        record.seq[start - 1:end],
        id=sequence_id,
        name=sequence_id,
        description="",
    )
    partial_file = f"{output_file}.part"
    try:
        SeqIO.write([segment], partial_file, "fasta")
        written_records = list(SeqIO.parse(partial_file, "fasta"))
        if (
            len(written_records) != 1
            or written_records[0].id != sequence_id
            or len(written_records[0]) != expected_length
        ):
            raise RuntimeError(
                f"Validation failed while writing {sequence_id}: "
                f"expected one {expected_length}-bp FASTA record."
            )
        os.replace(partial_file, output_file)
    except Exception:
        if os.path.exists(partial_file):
            os.remove(partial_file)
        raise


def process_reference_sequences(
    rows,
    output_folder,
    accession_index,
    sequence_id_index,
    coordinates_index,
):
    """Extract all coordinate ranges safely and return FASTA files to merge."""
    rows_by_accession = defaultdict(list)
    seen_sequence_ids = {}

    for row in rows:
        accession = validate_file_identifier(row[accession_index], ACCESSION_COLUMN)
        sequence_id = validate_file_identifier(
            row[sequence_id_index], SEQUENCE_ID_COLUMN
        )
        coordinates = row[coordinates_index]
        expected_sequence_id = make_sequence_id(accession, coordinates)
        if sequence_id != expected_sequence_id:
            raise ValueError(
                f"{SEQUENCE_ID_COLUMN} {sequence_id!r} does not match "
                f"accession/range {accession!r}, {coordinates!r}; "
                f"expected {expected_sequence_id!r}."
            )
        definition = (accession, coordinates)
        previous = seen_sequence_ids.get(sequence_id)
        if previous is not None and previous != definition:
            raise ValueError(
                f"Duplicate {SEQUENCE_ID_COLUMN} {sequence_id!r} represents both "
                f"{previous} and {definition}."
            )
        seen_sequence_ids[sequence_id] = definition
        rows_by_accession[accession].append(row)

    for accession, accession_rows in rows_by_accession.items():
        partial_rows = [
            row for row in accession_rows
            if row[coordinates_index] != FULL_LENGTH
        ]
        if not partial_rows:
            continue

        input_fasta = os.path.join(output_folder, f"{accession}.fasta")
        record = load_accession_record(input_fasta, accession)
        for row in partial_rows:
            start, end = map(int, row[coordinates_index].split("~"))
            sequence_id = row[sequence_id_index]
            output_fasta = os.path.join(output_folder, f"{sequence_id}.fasta")
            write_segment_fasta(record, start, end, sequence_id, output_fasta)

        has_full_length = any(
            row[coordinates_index] == FULL_LENGTH for row in accession_rows
        )
        if not has_full_length:
            os.remove(input_fasta)

    expected_files = []
    seen_files = set()
    for row in rows:
        if row[coordinates_index] == FULL_LENGTH:
            fasta_file = os.path.join(
                output_folder, f"{row[accession_index]}.fasta"
            )
        else:
            fasta_file = os.path.join(
                output_folder, f"{row[sequence_id_index]}.fasta"
            )
        if fasta_file not in seen_files:
            if not os.path.isfile(fasta_file) or os.path.getsize(fasta_file) == 0:
                raise RuntimeError(f"Expected FASTA file is missing or empty: {fasta_file}")
            expected_files.append(fasta_file)
            seen_files.add(fasta_file)

    return expected_files

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
    Robust reader for ICTV VMR csv. Only keeps the preferred reference ID
    ("VITAP sequence ID" when present, otherwise "Virus GENBANK accession")
    and the requested taxon-level column.
    Uses Python csv module-like parsing via polars? -> We'll do pure python for max safety.
    """
    accessions = []
    taxa = []

    with open(ictv_file, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        acc_key = mapping_id_column(reader.fieldnames or [])
        if taxon_level not in reader.fieldnames:
            raise KeyError(f"ICTV file missing column: {taxon_level}")

        for row in reader:
            acc = normalize_reference_id(row.get(acc_key) or "")
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

    # Extract genome accession from protein ID (format: {genome_id}_{orf_num}).
    # Strip only the trailing ORF number, then normalize a legacy version suffix.
    blast_results["qseqid_genome_id"] = blast_results["qseqid"].apply(
        lambda x: normalize_reference_id(x.rsplit("_", 1)[0])
    )
    blast_results["sseqid_genome_id"] = blast_results["sseqid"].apply(
        lambda x: normalize_reference_id(x.rsplit("_", 1)[0])
    )

    reference_key = mapping_id_column(ictv_data.columns)
    ictv_data[reference_key] = ictv_data[reference_key].apply(normalize_reference_id)

    conflicting = (
        ictv_data.groupby(reference_key, dropna=False)[taxon_level]
        .nunique(dropna=False)
    )
    conflicting = conflicting[conflicting > 1]
    if not conflicting.empty:
        examples = ", ".join(map(str, conflicting.index[:5]))
        raise ValueError(
            f"Conflicting {taxon_level} assignments for reference IDs: {examples}"
        )

    ictv_data = ictv_data.drop_duplicates(subset=[reference_key], keep="first")
    ictv_data.set_index(reference_key, inplace=True)

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
    input_file = args.vmr
    output_file = args.out
    db_name = args.db
    output_folder = "VMR_Genome"
    os.makedirs(output_folder, exist_ok=True)
    VMR_csv_file = output_file

    # Delete empty files and interrupted temporary downloads first.
    for root, dirs, files in os.walk(output_folder):
        for name in files:
            file_path = os.path.join(root, name)
            if name.endswith(".part") or os.path.getsize(file_path) == 0:
                os.remove(file_path)

    # ===== Reformatting VMR table =====
    with open(input_file, "r", encoding="utf-8") as infile, open(output_file, "w", encoding="utf-8", newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        input_header = next(reader)
        if not input_header or input_header[0].strip() != ACCESSION_COLUMN:
            raise ValueError(
                f"The first input column must be {ACCESSION_COLUMN!r}."
            )
        reserved_columns = {SEQUENCE_ID_COLUMN, START_END_COLUMN}
        if reserved_columns.intersection(input_header):
            raise ValueError(
                "The input appears to be already reformatted; please provide the original VMR CSV."
            )
        output_header = input_header + [SEQUENCE_ID_COLUMN, START_END_COLUMN]
        writer.writerow(output_header)

        # Processing data rows
        for row_number, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != len(input_header):
                raise ValueError(
                    f"CSV row {row_number} has {len(row)} columns; "
                    f"expected {len(input_header)}."
                )
            virus_ids = row[0].split(";")
            for virus_id in virus_ids:
                cleaned_id = clean_virus_id(virus_id)
                if cleaned_id:
                    cleaned_id, start_end_sites = extract_start_end_sites(cleaned_id)
                    sequence_id = make_sequence_id(cleaned_id, start_end_sites)
                    taxonomy_row = fill_empty_cells(
                        [cleaned_id] + row[1:], input_header
                    )
                    writer.writerow(
                        taxonomy_row + [sequence_id, start_end_sites]
                    )

    # Y/N confirmation
    while True:
        user_input = input(f"=====\n[WARNING] You need to double check the VMR dataframe generated in last step [{output_file}], and ensure all start/end information is correct. Please pay special attention to these categories: Peduoviridae (eg. AE006468), Belpaoviridae (eg. LK928904), all GTA-viriform (Bartogtaviriformidae and Rhodogtaviriformidae). The related metedata of these categories is incorrect in ICTV_VMR-MSL38_210426. Please make sure that all mistakes have been fixed and saved, and enter 'N' to abort, then rerun the program. If no mistakes existed yet, enter 'Y' to continue (Y/N): \n=====\n")
        if user_input.lower() == 'y':
            break
        elif user_input.lower() == 'n':
            print("[INFO] Exiting the program.")
            exit()

    # ===== Downloading genome FASTA =====
    # ===== Load VMR rows =====
    with open(VMR_csv_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        vmr_header = next(reader)
        rows = [row for row in reader if row and any(cell.strip() for cell in row)]

    try:
        accession_index = vmr_header.index(ACCESSION_COLUMN)
        sequence_id_index = vmr_header.index(SEQUENCE_ID_COLUMN)
        coordinates_index = vmr_header.index(START_END_COLUMN)
    except ValueError as error:
        raise ValueError(
            f"Reformatted VMR CSV must contain {ACCESSION_COLUMN!r}, "
            f"{SEQUENCE_ID_COLUMN!r}, and {START_END_COLUMN!r}."
        ) from error

    malformed_rows = [
        index + 2 for index, row in enumerate(rows)
        if len(row) != len(vmr_header)
    ]
    if malformed_rows:
        raise ValueError(
            f"Malformed rows in reformatted VMR CSV: {malformed_rows[:10]}"
        )

    # Coordinates may have been corrected during the confirmation pause.
    # Regenerate the derived sequence IDs and persist the synchronized table.
    for row in rows:
        row[sequence_id_index] = make_sequence_id(
            row[accession_index], row[coordinates_index]
        )
    with open(VMR_csv_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(vmr_header)
        writer.writerows(rows)

    unique_accessions = list(dict.fromkeys(
        validate_file_identifier(row[accession_index], ACCESSION_COLUMN)
        for row in rows
    ))

    VMR_csv_file = output_file

    counter_lock = Lock()

    downloaded_ids = {
        f[:-6] for f in os.listdir(output_folder)
        if f.endswith(".fasta") and os.path.getsize(os.path.join(output_folder, f)) > 0
    }

    with ThreadPoolExecutor(max_workers=3) as executor:
        progress_bar = tqdm(total=len(unique_accessions), desc="Downloading genomes")
        futures = [
            executor.submit(
                download_and_process_genome,
                virus_id,
                output_folder,
                downloaded_ids,
                progress_bar,
                counter_lock,
            )
            for virus_id in unique_accessions
        ]

        try:
            for future in as_completed(futures):
                future.result()
        finally:
            progress_bar.close()

    reference_fasta_files = process_reference_sequences(
        rows,
        output_folder,
        accession_index,
        sequence_id_index,
        coordinates_index,
    )

    print(
        f"[INFO] All files successfully downloaded and processed "
        f"({len(reference_fasta_files)} unique reference sequences)"
    )

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
            for fasta_file in reference_fasta_files:
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
