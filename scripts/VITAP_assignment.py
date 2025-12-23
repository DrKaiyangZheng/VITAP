import csv
import argparse
from collections import Counter, defaultdict
import os
import glob
import subprocess
import time
import math
import pyrodigal
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import polars as pl

# Your fallback function is provided by another script/module
from uniref_taxa_fallback import uniref_taxa_fallback


# =========================================================
# Helpers
# =========================================================

TAXON_CATEGORIES = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
TAXON_INDEX = {t: i for i, t in enumerate(TAXON_CATEGORIES)}


def deplication_check(fasta_file: str):
    """Check duplicated Genome_IDs in the input FASTA file."""
    genome_ids = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                genome_ids.append(line[1:].split()[0])

    duplicated_ids = [k for k, v in Counter(genome_ids).items() if v > 1]
    if duplicated_ids:
        print("[WARNING] Duplicated Genome_IDs found in the input FASTA file:")
        for d in duplicated_ids:
            print(f"-> {d}")
        raise SystemExit(1)


def extract_short_sequences(fasta_file: str, protein_file: str):
    """Find genomes ignored by prodigal and return those genome sequences."""
    fasta_ids = {r.id for r in SeqIO.parse(fasta_file, "fasta")}

    protein_ids = set()
    for record in SeqIO.parse(protein_file, "fasta"):
        genome_id = record.id.rsplit("_", 1)[0]
        protein_ids.add(genome_id)

    short_ids = fasta_ids - protein_ids
    if not short_ids:
        return []

    return [r for r in SeqIO.parse(fasta_file, "fasta") if r.id in short_ids]

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

def orf_count(gff_file: str) -> pl.DataFrame:
    """
    Safe ORF counter for GFF (plain parsing).
    Returns: polars DataFrame with columns: id, ORF_number
    """
    counter = Counter()
    with open(gff_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            gid = line.split("\t", 1)[0]
            counter[gid] += 1

    return pl.DataFrame({"id": list(counter.keys()), "ORF_number": list(counter.values())})


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def is_blank_or_bracketed(label: str) -> bool:
    """Return True if label is None/empty/nan or contains [] which must be replaced by '-'."""
    if label is None:
        return True
    s = str(label).strip()
    if s == "" or s.lower() in {"nan", "none"}:
        return True
    if "[" in s or "]" in s:
        return True
    return False


def normalize_label(label: str) -> str:
    """Convert invalid labels to '-' and keep others as string."""
    return "-" if is_blank_or_bracketed(label) else str(label)


def find_vmr_base_name(db_dir: str) -> str:
    """
    Resolve the correct DIAMOND db for ICTV/VMR blast by using the *.gff base name.
    Example: VMR_genome_VMR-MSL_v40.gff -> base 'VMR_genome_VMR-MSL_v40'
             then DIAMOND db should be 'VMR_genome_VMR-MSL_v40.dmnd'
    """
    gff_files = glob.glob(os.path.join(db_dir, "*.gff"))
    if not gff_files:
        raise FileNotFoundError(f"No *.gff found under db_dir: {db_dir}")

    # If multiple gff exist, choose the first one deterministically (sorted)
    gff_files.sort()
    gff_path = gff_files[0]
    base = os.path.splitext(os.path.basename(gff_path))[0]
    return base


def resolve_ictv_dmnd(db_dir: str) -> str:
    base = find_vmr_base_name(db_dir)
    dmnd_path = os.path.join(db_dir, f"{base}.dmnd")
    if not Path(dmnd_path).is_file():
        raise FileNotFoundError(
            f"Expected ICTV/VMR DIAMOND db not found: {dmnd_path}\n"
            f"(Derived from gff base name '{base}')"
        )
    return dmnd_path


def resolve_vmr_csv(db_dir: str) -> str:
    csv_files = glob.glob(os.path.join(db_dir, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No *.csv found under db_dir: {db_dir}")
    csv_files.sort()
    return csv_files[0]


def resolve_vmr_fasta(db_dir: str) -> str:
    fasta_files = glob.glob(os.path.join(db_dir, "*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(f"No *.fasta found under db_dir: {db_dir}")
    fasta_files.sort()
    return fasta_files[0]

def normalize_lineage_field(taxon: str) -> str:
    if taxon is None:
        return "-"
    taxon = taxon.strip()
    if taxon == "":
        return "-"
    if "[" in taxon:
        return "-"
    return taxon

def normalize_lineage_string(lineage: str) -> str:
    """
    Normalize a semicolon-separated lineage string.
    Any field containing '[' is replaced by '-'.
    """
    parts = lineage.split(";")
    parts = [normalize_lineage_field(p) for p in parts]
    return ";".join(parts)

# =========================================================
# taxonomy_assigning (same logic as your stable version)
# =========================================================

def taxonomy_assigning(
    blast_results_file,
    ictv_file,
    genome_cutoff_file,
    taxon_level,
    genome_length_file,   # kept for interface compatibility; not used here
    gff_file,             # polars DF with id, ORF_number
    output_file,
    graph_output_file
):
    # ---- ORF_number map ----
    orf_number = {}
    for r in gff_file.to_dicts():
        orf_number[str(r["id"])] = int(r["ORF_number"])

    # ---- ICTV mapping: accession -> taxon ----
    accession_to_taxon = {}
    ictv_key = "Virus GENBANK accession"
    with open(ictv_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row.get(ictv_key, "")
            if acc is None:
                continue
            accession_to_taxon[str(acc)] = row.get(taxon_level)

    # ---- Read BLAST: qseqid sseqid bitscore ----
    blast_rows = []
    with open(blast_results_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            qseqid, sseqid, bitscore_s = parts[0], parts[1], parts[2]
            bitscore = safe_float(bitscore_s)

            qgenome = qseqid.rsplit("_", 1)[0]
            s_acc = sseqid.split(".", 1)[0]
            staxon = accession_to_taxon.get(s_acc)

            hit_type = "self_hit" if qseqid == sseqid else "other"

            blast_rows.append({
                "qseqid": qseqid,
                "qgenome": qgenome,
                "staxon": staxon,
                "bitscore": bitscore,
                "hit_type": hit_type
            })

    # ---- Step2: top hit taxa per qseqid (exclude self_hit) ----
    top_hit_taxa = {}
    for r in blast_rows:
        if r["hit_type"] == "other":
            q = r["qseqid"]
            if q not in top_hit_taxa:
                top_hit_taxa[q] = r["staxon"]

    # ---- Step3: mark top_hit ----
    for r in blast_rows:
        if r["hit_type"] != "self_hit":
            if r["staxon"] == top_hit_taxa.get(r["qseqid"]):
                r["hit_type"] = "top_hit"

    # ---- Step4: w1 ----
    for r in blast_rows:
        if r["hit_type"] == "top_hit":
            r["w1"] = 1.2
        elif r["hit_type"] == "self_hit":
            r["w1"] = 1.0
        else:
            r["w1"] = 0.8

    # ---- Step5: w2 by taxa frequency ----
    total_count = Counter()
    taxon_count = Counter()
    for r in blast_rows:
        total_count[r["qseqid"]] += 1
        taxon_count[(r["qseqid"], r["staxon"])] += 1

    w2_map = {}
    for (q, taxa), c in taxon_count.items():
        perc = c / total_count[q] if total_count[q] else 0.0
        w2_map[(q, taxa)] = 1.2 if perc > 0.5 else 1.0

    for r in blast_rows:
        if r["hit_type"] == "self_hit":
            r["w2"] = 1.0
        else:
            r["w2"] = w2_map.get((r["qseqid"], r["staxon"]), 1.0)

    # ---- taxon_bitscore ----
    bitscore_col = f"{taxon_level}_bitscore"
    for r in blast_rows:
        r[bitscore_col] = r["bitscore"] * r["w1"] * r["w2"]

    # ---- ORF occurrence count per (genome, taxa) ----
    seen = set()
    genome_taxon_orfcount = Counter()
    for r in blast_rows:
        key = (r["qseqid"], r["qgenome"], r["staxon"])
        if key in seen:
            continue
        seen.add(key)
        genome_taxon_orfcount[(r["qgenome"], r["staxon"])] += 1

    # ---- Aggregate sum and count by (genome, taxa, occ) ----
    sum_bitscore = defaultdict(float)
    count_hits = Counter()
    for r in blast_rows:
        g = r["qgenome"]
        t = r["staxon"]
        occ = genome_taxon_orfcount[(g, t)]
        k = (g, t, occ)
        sum_bitscore[k] += float(r[bitscore_col])
        count_hits[k] += 1

    total_count_per_genome = Counter()
    for (g, t, occ), c in count_hits.items():
        total_count_per_genome[g] += c

    # ---- taxon_score ----
    taxon_score_rows = []
    for (g, t, occ), c in count_hits.items():
        if c == 0:
            continue
        total_c = total_count_per_genome[g]
        perc_sseqid_taxa = 10 * (c / total_c) if total_c else 0.0
        orf_num = orf_number.get(g, 0)
        if orf_num == 0:
            score = 0.0
        else:
            score = (sum_bitscore[(g, t, occ)] / c) * ((occ / orf_num) ** 2) * perc_sseqid_taxa

        taxon_score_rows.append((g, t, score))

    # ---- thresholds ----
    cutoff_map = {}
    cut_col = f"{taxon_level}_score_cut-off"
    with open(genome_cutoff_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxa = row.get(taxon_level)
            cs = row.get(cut_col)
            if taxa is None:
                continue
            cutoff_map[str(taxa)] = safe_float(cs)

    # ---- output ----
    with open(output_file, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Genome_ID", taxon_level, f"{taxon_level}_weight"])
        for g, t, score in taxon_score_rows:
            cutoff = cutoff_map.get(str(t), float("nan"))
            if cutoff == cutoff and cutoff != 0:
                weight = score / cutoff
            else:
                weight = float("nan")
            w.writerow([g, t, weight])


# =========================================================
# Multipartite lineage generator + UniRef filter + UniRef-based fill
# =========================================================

def multipartite_taxonomic_graph_generator(result_dir, vmr_mapping_file, taxon_categories):
    """
    Generate best lineages and all candidate lineages.
    Applies UniRef90 fallback filtering as a pre-filter; and if all candidates are filtered out,
    a UniRef90-based lineage will be produced (lineage_score empty, confidence 'UniRef90-based').

    NOTE: Final hard "guarantee inclusion" + high-rank mapping from VMR CSV is done in assignment() as a last-resort.
    """

    def is_na(x):
        if x is None:
            return True
        if isinstance(x, str) and x.strip().lower() in {"nan", "none", ""}:
            return True
        try:
            return math.isnan(float(x))
        except Exception:
            return False

    def to_float(x):
        try:
            return float(x)
        except Exception:
            return float("nan")

    def mean_skip_na(vals):
        nums = [to_float(v) for v in vals if not is_na(v)]
        return float(sum(nums) / len(nums)) if nums else float("nan")

    def process_row(row, weights, lineages):
        # Step A: if any weight >= 0.6, return average from that level to the end and the lineage string
        for i, wcol in enumerate(weights):
            w = row.get(wcol)
            if (not is_na(w)) and to_float(w) >= 0.6:
                avg = mean_skip_na([row.get(c) for c in weights[i:]])
                lin = ";".join([str(row.get(c)) for c in lineages[i:]])
                return avg, lin

        # Step B: all < 0.6, compute total_max_average from Species end, keep at least 3 levels
        total_max_average = -1
        for i in range(len(weights), 2, -1):
            curr = [row.get(c) for c in weights[-i:]]
            if all(not is_na(v) for v in curr):
                curr_avg = mean_skip_na(curr)
                if curr_avg > total_max_average:
                    total_max_average = curr_avg

        # Step C: choose the deepest combination that meets >= 0.95*total_max_average and count>=3
        max_average = -1
        max_count = 0
        max_lineage = None

        for i in range(len(weights), 2, -1):
            curr = [row.get(c) for c in weights[-i:]]
            if all(not is_na(v) for v in curr):
                curr_avg = mean_skip_na(curr)
                if curr_avg >= 0.95 * total_max_average and i >= 3:
                    if i > max_count:
                        max_average = curr_avg
                        max_count = i
                        max_lineage = ";".join([str(row.get(c)) for c in lineages[-i:]])

        if max_average != -1:
            return max_average, max_lineage
        return "NaN", "NaN"

    def pad_lineage(lineage):
        elements = str(lineage).split(";")
        return ";".join(["-"] * (8 - len(elements)) + elements)

    def count_dashes(lineage):
        return str(lineage).count("-")

    def replace_first_two_elements(lineage):
        elements = str(lineage).split(";")
        if len(elements) < 2:
            elements = (["-"] * (2 - len(elements))) + elements
        if elements[0] != "-" or elements[1] != "-":
            elements[0] = "-"
            elements[1] = "-"
        return ";".join(elements)

    print("[INFO] Importing taxonomic bipartite graphs and ICTV taxonomic hierarchy")

    # ---- UniRef90 fallback map: gid -> level -> set(taxa_name) ----
    fallback_map = defaultdict(lambda: defaultdict(set))
    fallback_file = os.path.join(result_dir, "target_uniref90_taxa_fallback.out")
    if Path(fallback_file).is_file():
        with open(fallback_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gid = row.get("genome_id")
                taxa = row.get("taxa_name")
                level = row.get("taxon_level")
                if gid and taxa and level:
                    fallback_map[str(gid)][str(level)].add(str(taxa))

    # ---- 1) Read Species bipartite graph ----
    species_annotation_file = os.path.join(result_dir, "Species_bipartite.graph")
    species_rows = []
    with open(species_annotation_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species_rows.append(row)

    # ---- 2) Read VMR mapping CSV: Species -> higher lineage levels ----
    vmr_map = {}
    with open(vmr_mapping_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sp = row.get("Species")
            if sp is None:
                continue
            vmr_map[str(sp)] = {k: row.get(k) for k in TAXON_CATEGORIES}

    # ---- 3) Merge species rows with VMR lineage ----
    merged = []
    for r in species_rows:
        sp = r.get("Species")
        lineage = vmr_map.get(str(sp), {})
        rec = {
            "Genome_ID": r.get("Genome_ID"),
            "Species": sp,
            "Species_weight": to_float(r.get("Species_weight")),
        }
        for t in taxon_categories[1:]:
            rec[t] = lineage.get(t)
            rec[f"{t}_weight"] = None
        merged.append(rec)

    # ---- 4) Fill weights for each rank ----
    for taxa in taxon_categories[1:]:
        unit_file = os.path.join(result_dir, f"{taxa}_bipartite.graph")
        key2w = {}
        with open(unit_file, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gid = row.get("Genome_ID")
                tval = row.get(taxa)
                w = row.get(f"{taxa}_weight")
                key2w[(str(gid), str(tval))] = to_float(w)

        for rec in merged:
            gid = str(rec.get("Genome_ID"))
            tval = rec.get(taxa)
            rec[f"{taxa}_weight"] = key2w.get((gid, str(tval)))

    weights = [f"{t}_weight" for t in taxon_categories]
    lineages = list(taxon_categories)

    print(
        "[INFO] Calculating average taxonomic score for target genomes: \n"
        "  -> finding the lowest effective taxonomic level\n"
        "  -> determining the optimal taxonomic hierarchy termination point"
    )

    # ---- 5) Compute lineage_score and lineage ----
    for rec in merged:
        score, lineage = process_row(rec, weights, lineages)
        rec["lineage_score"] = score
        rec["lineage"] = lineage

    # ---- 6) Filter: at least 4 levels ----
    merged = [r for r in merged if str(r.get("lineage")).count(";") >= 3]

    # ---- 7) pad lineage and deduplicate ----
    all_lineage = {}
    for r in merged:
        gid = r["Genome_ID"]
        lin = pad_lineage(r["lineage"])
        score = r["lineage_score"]
        all_lineage[(gid, lin, score)] = True

    # ---- 7.5) UniRef90 pre-filter on candidate lineages ----
    all_lineage_rows = []
    filtered_out_by_uniref = set()  # genomes that had fallback constraints but all candidates were removed

    # group candidates by genome to later detect "all removed"
    candidates_by_gid = defaultdict(list)
    for (gid, lin, score) in all_lineage.keys():
        candidates_by_gid[gid].append((lin, score))

    for gid, cand_list in candidates_by_gid.items():
        kept_any = False
        for (lin, score) in cand_list:
            if gid in fallback_map:
                elements = lin.split(";")
                keep = True
                for level, allowed_taxa in fallback_map[gid].items():
                    if level not in TAXON_INDEX:
                        continue
                    idx = TAXON_INDEX[level]
                    if idx < len(elements) and elements[idx] not in allowed_taxa:
                        keep = False
                        break
                if not keep:
                    continue

            kept_any = True
            all_lineage_rows.append({
                "Genome_ID": gid,
                "lineage": lin,
                "lineage_score": score,
                "dash_count": count_dashes(lin),
            })

        if gid in fallback_map and (not kept_any):
            filtered_out_by_uniref.add(gid)

    print("[INFO] Determine best-fit taxonomic hierarchy")

    # ---- 8) best_high: score>=0.9 choose minimal dash_count ----
    best_high = {}
    for r in all_lineage_rows:
        sc = r["lineage_score"]
        if is_na(sc) or to_float(sc) < 0.9:
            continue
        g = r["Genome_ID"]
        if g not in best_high or r["dash_count"] < best_high[g]["dash_count"]:
            best_high[g] = r

    # ---- 9) best_low: score<0.9 and not in best_high, choose max score ----
    best_low = {}
    for r in all_lineage_rows:
        sc = r["lineage_score"]
        if is_na(sc) or to_float(sc) >= 0.9:
            continue
        g = r["Genome_ID"]
        if g in best_high:
            continue
        if g not in best_low or to_float(sc) > to_float(best_low[g]["lineage_score"]):
            best_low[g] = r

    best = list(best_high.values()) + list(best_low.values())

    print("[INFO] Assigning confidence level to taxonomic hierarchy")

    # ---- 10) confidence + replace_first_two_elements for Low/Medium ----
    best_lineage_rows = []
    for r in best:
        sc = to_float(r["lineage_score"])
        if sc >= 0.9:
            conf = "High-confidence"
        elif sc >= 0.1:
            conf = "Medium-confidence"
        else:
            conf = "Low-confidence"

        lin = r["lineage"]
        if conf in ("Low-confidence", "Medium-confidence"):
            lin = replace_first_two_elements(lin)

        best_lineage_rows.append({
            "Genome_ID": r["Genome_ID"],
            "lineage": lin,
            "lineage_score": r["lineage_score"],
            "Confidence_level": conf
        })

    # NOTE:
    # UniRef90-based "fill" is not done here, because the full mapping (higher ranks) is performed
    # more safely in assignment() using the VMR CSV and fallback rows (participation_index, taxon_level).

    return best_lineage_rows, all_lineage_rows


# =========================================================
# UniRef90-based fallback: guarantee inclusion of fallback genome IDs
# =========================================================

def normalize_taxon_level(level: str):
    """
    Normalize taxon level name to standard TAXON_CATEGORIES entry.
    Returns normalized level name or None if invalid.
    """
    if level is None:
        return None

    s = str(level).strip()
    if not s:
        return None

    # case-insensitive match
    s_fold = s.casefold()
    for t in TAXON_CATEGORIES:
        if s_fold == t.casefold():
            return t

    return None

def load_vmr_csv_rows(vmr_csv_file: str):
    with open(vmr_csv_file, "r", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)

def build_high_rank_lineage_from_vmr(vmr_rows, taxon_level: str, taxa_name: str):
    taxon_level_norm = normalize_taxon_level(taxon_level)
    if taxon_level_norm is None:
        return ["-"] * 8

    idx = TAXON_INDEX[taxon_level_norm]

    # case-insensitive taxa_name match
    taxa_name_norm = str(taxa_name).strip().casefold()

    match = None
    for r in vmr_rows:
        val = r.get(taxon_level_norm)
        if val is None:
            continue
        if str(val).strip().casefold() == taxa_name_norm:
            match = r
            break

    lineage = ["-"] * 8
    lineage[idx] = normalize_label(taxa_name)

    if match:
        for j in range(idx + 1, 8):
            lvl = TAXON_CATEGORIES[j]
            lineage[j] = normalize_label(match.get(lvl))

    return lineage

def parse_uniref_fallback_file(fallback_file: str):
    """
    Return dict gid -> list of fallback records
    record: {taxa_name, participation_index, taxon_level}
    """
    gid2recs = defaultdict(list)
    if not Path(fallback_file).is_file():
        return gid2recs

    with open(fallback_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row.get("genome_id")
            if not gid:
                continue
            gid2recs[str(gid)].append({
                "taxa_name": row.get("taxa_name"),
                "participation_index": row.get("participation_index"),
                "taxon_level": row.get("taxon_level"),
            })
    return gid2recs


def ensure_fallback_genomes_in_outputs(
    best_rows: list[dict],
    all_rows: list[dict],
    vmr_csv_file: str,
    fallback_file: str
):
    """
    If a Genome_ID appears in target_uniref90_taxa_fallback.out but is missing
    from best_rows and all_rows, append a UniRef90-based record:

    - lineage: built using VMR CSV (higher ranks above taxon_level)
    - lineage_score/participation_index: participation_index
    - Confidence_level: UniRef90-based
    - Any label containing [] is replaced by '-'
    """

    vmr_rows = load_vmr_csv_rows(vmr_csv_file)
    gid2recs = parse_uniref_fallback_file(fallback_file)

    best_gids = {r["Genome_ID"] for r in best_rows}
    all_gids = {r["Genome_ID"] for r in all_rows}

    appended_best = 0
    appended_all = 0

    for gid, recs in gid2recs.items():
        # Only act if missing from outputs
        need_best = gid not in best_gids
        need_all = gid not in all_gids
        if not (need_best or need_all):
            continue

        # Choose the record with max participation_index as the representative fallback
        # (This is a reasonable deterministic policy.)
        best_rec = None
        best_pi = -1.0
        for r in recs:
            pi = safe_float(r.get("participation_index"))
            if pi > best_pi:
                best_pi = pi
                best_rec = r

        if best_rec is None:
            continue

        taxon_level = str(best_rec.get("taxon_level") or "").strip()
        taxa_name = str(best_rec.get("taxa_name") or "").strip()
        pi_val = best_rec.get("participation_index")

        lineage_list = build_high_rank_lineage_from_vmr(vmr_rows, taxon_level, taxa_name)
        lineage_str = ";".join(lineage_list)

        # Append to all_lineages (no confidence column there)
        if need_all:
            all_rows.append({
                "Genome_ID": gid,
                "lineage": lineage_str,
                "lineage_score": pi_val,  # will be renamed in output header
                "dash_count": lineage_str.count("-"),  # keep internal consistency
            })
            all_gids.add(gid)
            appended_all += 1

        # Append to best_determined_lineages
        if need_best:
            best_rows.append({
                "Genome_ID": gid,
                "lineage": lineage_str,
                "lineage_score": pi_val,  # will be renamed in output header
                "Confidence_level": "UniRef90-based"
            })
            best_gids.add(gid)
            appended_best += 1

    return appended_best, appended_all


# =========================================================
# assignment()
# =========================================================

def assignment(args):
    input_fasta = args.fasta
    db_dir = args.db.rstrip("/") + "/"
    threads = str(args.cpu)
    result_dir = args.out.rstrip("/") + "/"

    os.makedirs(result_dir, exist_ok=True)

    start_wall = time.time()
    start_cpu = time.process_time()

    print("===== The VITAP (Viral Taxonomic Assignment Pipeline) v.1.8 =====")
    print("[INFO] Temporary files will be generated in the result directory.")

    # ---- Check duplicate ID ----
    deplication_check(input_fasta)

    # ---- ORF calling ----
    print(f"===== ORF calling for {input_fasta} =====")
    merge_fasta = os.path.join(result_dir, "merge_genome.fasta")
    merge_gff = os.path.join(result_dir, "merge_genome.gff")
    merge_faa = os.path.join(result_dir, "merge_genome.faa")
    selected_ref = os.path.join(result_dir, "ICTV_selected_genomes.fasta")

    vmr_fasta = resolve_vmr_fasta(db_dir)
    vmr_csv = resolve_vmr_csv(db_dir)

    # Merge target genome + random refs, then run prodigal if merge_faa not exists
    if not Path(merge_faa).is_file():
        print("[INFO] Ten ICTV reference genomes were randomly selected and merged to your genome set!")
        subprocess.run(["seqkit", "sample", "-n", "5", "-o", selected_ref, vmr_fasta], check=True)

        with open(merge_fasta, "w") as out:
            out.write(Path(input_fasta).read_text())
            out.write(Path(selected_ref).read_text())

        run_pyrodigal(
            fasta_in=merge_fasta,
            faa_out=merge_faa,
            gff_out=merge_gff,
        )
    else:
        print(f"[INFO] Using existing file: {merge_faa}")

    # Handle short sequences ignored by prodigal
    short_sequences = extract_short_sequences(merge_fasta, merge_faa)
    if short_sequences:
        print(f"[INFO] Processing short sequences in {merge_fasta} ignored by prodigal.")
        short_genome_file = os.path.join(result_dir, "target_short_genome.fasta")
        SeqIO.write(short_sequences, short_genome_file, "fasta")
        short_faa_file = os.path.join(result_dir, "target_short_genome.faa")
        subprocess.run(["seqkit", "translate", "-f", "6", "-F", "--clean", "-o", short_faa_file, short_genome_file], check=True)

        with open(merge_faa, "a") as final_out, open(short_faa_file, "r") as sf:
            final_out.write(sf.read())

        os.remove(short_genome_file)
        os.remove(short_faa_file)
    else:
        print("[INFO] Genome and ORF files have consistent non-redundant FASTA IDs. √")

    # ---- Genome length (kept for legacy; not used in the stable formula here) ----
    length_file = os.path.join(result_dir, "merged_genome_length.tsv")
    if not Path(length_file).is_file():
        print("[INFO] Length statistic for merged_genome.fasta")
        subprocess.run(["seqkit", "fx2tab", "-l", "-n", "-i", "-H", "-o", length_file, merge_fasta], check=True)
    else:
        print("[INFO] The merged_genome_length.tsv exists, skipping.")

    # ---- ORF count ----
    print("[INFO] Statistic of the number of ORF per genome")
    orf_df = orf_count(merge_gff)

    # ---- ICTV DIAMOND blastp ----
    print("===== Aligning to ICTV reference database based on BLAST-algorithm =====")
    ictv_dmnd = resolve_ictv_dmnd(db_dir)
    target_blast_fp = os.path.join(result_dir, "target_ICTV_blastp.align")

    if not Path(target_blast_fp).is_file():
        print(f"[INFO] Please delete the output file ({target_blast_fp}) if this step aborted!")
        print(f"[INFO] Please keep the output file ({target_blast_fp}) if this step finished!")
        subprocess.run(
            ["diamond", "blastp", "--quiet", "-p", threads, "-q", merge_faa, "-d", ictv_dmnd,
             "--sensitive", "-f", "6", "qseqid", "sseqid", "bitscore",
             "-o", target_blast_fp, "-k", "1000", "--max-hsps", "1", "-e", "1e-3"],
            check=True
        )
    else:
        print(f"[INFO] Using existing file: {target_blast_fp}")

    # ---- Taxonomy graphs ----
    print("===== Locating contig on taxonomic framework and generating genome community network =====")
    for taxa in TAXON_CATEGORIES:
        genome_cutoff_file = os.path.join(db_dir, f"{taxa}_genome.threshold")
        taxa_annot_file = os.path.join(result_dir, f"{taxa}_bipartite.graph")
        if not Path(taxa_annot_file).is_file():
            print(f"[INFO] Generating {taxa} graph for target genomes")
            taxonomy_assigning(
                target_blast_fp,
                vmr_csv,
                genome_cutoff_file,
                taxa,
                length_file,
                orf_df,
                taxa_annot_file,
                taxa_annot_file
            )
        else:
            print(f"[INFO] Using existing file: {taxa_annot_file} for {taxa} graph")

    # ---- UniRef90 filtering step ----
    print("===== Filtering lineage information based on UniRef90 =====")
    uniref_blast_out = os.path.join(result_dir, "target_uniref90_blastp.align")
    uniref_fallback_out = os.path.join(result_dir, "target_uniref90_taxa_fallback.out")
    uniref_dmnd = os.path.join(db_dir, "uniref90.dmnd")
    taxdmp_dir = os.path.join(db_dir, "taxdmp")

    if not Path(uniref_fallback_out).is_file():
        if not Path(uniref_blast_out).is_file():
            subprocess.run([
                "diamond", "blastp",
                "-q", merge_faa,
                "-d", uniref_dmnd,
                "-e", "1e-3",
                "-o", uniref_blast_out,
                "--max-hsps", "1",
                "-k", "100",
                "--outfmt", "6", "qseqid", "sseqid", "bitscore", "staxids",
                "--ultra-sensitive",
                "--taxonlist", "10239",
                "--threads", "10",
                "--quiet"
            ], check=True)

        uniref_taxa_fallback(
            uniref_blast_out,
            taxdmp_dir,
            merge_faa,
            uniref_fallback_out
        )
    else:
        print(f"[INFO] Using existing UniRef90 fallback file: {uniref_fallback_out}")

    # ---- Multipartite inference ----
    print("===== Determining viral lineages based on multi-partite graph ===== ")
    best_rows, all_rows = multipartite_taxonomic_graph_generator(result_dir, vmr_csv, TAXON_CATEGORIES)

    # ---- Low-confidence handling (only affects best_rows) ----
    if not args.low_conf:
        best_rows = [r for r in best_rows if r.get("Confidence_level") != "Low-confidence"]

    # ---- HARD fallback guarantee inclusion (your new requirement) ----
    appended_best, appended_all = ensure_fallback_genomes_in_outputs(
        best_rows=best_rows,
        all_rows=all_rows,
        vmr_csv_file=vmr_csv,
        fallback_file=uniref_fallback_out
    )
    if appended_best or appended_all:
        print(f"[INFO] UniRef90-based fallback appended: best={appended_best}, all={appended_all}")

    # ---- Write outputs ----
    print("===== Exporting results ===== ")
    best_out = os.path.join(result_dir, "best_determined_lineages.tsv")
    all_out = os.path.join(result_dir, "all_lineages.tsv")

    # Rename score column in output headers
    score_col = "lineage_score/participation_index"

    # best: Genome_ID, lineage, score, Confidence
    with open(best_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Genome_ID", "lineage", score_col, "Confidence_level"])
        for r in best_rows:
            clean_lineage = normalize_lineage_string(str(r.get("lineage")))

            w.writerow([
                r.get("Genome_ID"),
                clean_lineage,
                r.get("lineage_score"),
                r.get("Confidence_level")
            ])

    # all: Genome_ID, lineage, score
    with open(all_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Genome_ID", "lineage", score_col])
        for r in all_rows:
            clean_lineage = normalize_lineage_string(str(r.get("lineage")))

            w.writerow([
                r.get("Genome_ID"),
                clean_lineage,
                r.get("lineage_score")
            ])

    # ---- Timing ----
    end_wall = time.time()
    end_cpu = time.process_time()
    print(f"[INFO] Time-consuming (wall clock): {(end_wall - start_wall) / 3600:.1f} hours")
    print(f"[INFO] Time-consuming (CPU): {(end_cpu - start_cpu) / 3600:.1f} hours")
    print("===== All done =====")


if __name__ == "__main__":
    assignment()