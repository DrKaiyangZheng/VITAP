#!/usr/bin/env python3
import sys
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
import pandas as pd
import os

THRESHOLD = 0.25


# --------------------------------------------------
# taxonkit: taxid -> phylum, class, order
# --------------------------------------------------
def run_taxonkit(taxids, taxdump):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tf:
        for t in taxids:
            tf.write(t + "\n")
        tf.flush()
        taxid_file = tf.name

    cmd = (
        f"taxonkit lineage --data-dir {taxdump} --show-lineage-taxids {taxid_file} | "
        f"taxonkit reformat --data-dir {taxdump} -f '{{p}},{{c}},{{o}}'"
    )

    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, check=True
    )

    taxid2tax = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        taxid = parts[0]
        taxstr = parts[-1]
        ranks = taxstr.split(",")
        while len(ranks) < 3:
            ranks.append("NA")
        ranks = [r if r else "NA" for r in ranks[:3]]
        taxid2tax[taxid] = tuple(ranks)

    Path(taxid_file).unlink()
    return taxid2tax


# --------------------------------------------------
# prodigal FAA: genome_id -> num_orf
# --------------------------------------------------
def count_orfs_from_prodigal(faa_file):
    genome2orf = defaultdict(int)
    with open(faa_file) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                genome_id = seq_id.rsplit("_", 1)[0]
                genome2orf[genome_id] += 1
    return genome2orf


# --------------------------------------------------
# contextual NA filling
# --------------------------------------------------
def fill_na_with_parent(df):
    df["class"] = df.apply(
        lambda r: r["class"]
        if r["class"] != "NA"
        else (f"c.{r['phylum']}" if r["phylum"] != "NA" else "NA"),
        axis=1,
    )
    df["order"] = df.apply(
        lambda r: r["order"]
        if r["order"] != "NA"
        else (f"o.{r['class']}" if r["class"] != "NA" else "NA"),
        axis=1,
    )
    return df


# --------------------------------------------------
# rank aggregation
# --------------------------------------------------
def build_rank_raw(df, rank):
    g = (
        df.groupby(["genome_id", rank], as_index=False)
        .agg(
            participation_num=("protein_id", "nunique"),
            sum_bitscore=("bitscore", "sum"),
        )
    )

    genome_orf = (
        df[["genome_id", "num_orf"]]
        .drop_duplicates(subset=["genome_id"])
    )

    g = g.merge(genome_orf, on="genome_id", how="left")
    g["participation_freq"] = g["participation_num"] / g["num_orf"]

    return g


# --------------------------------------------------
# normalize + participation_index
# --------------------------------------------------
def finalize_rank_df(raw_df):
    if raw_df.empty:
        return raw_df

    df = raw_df.copy()
    df["norm_sum_bitscore"] = (
        df["sum_bitscore"]
        / df.groupby("genome_id")["sum_bitscore"].transform("max")
    )
    df["participation_index"] = (
        df["participation_freq"] * df["norm_sum_bitscore"]
    )
    return df


# --------------------------------------------------
# semantic fallback
# --------------------------------------------------
def semantic_fallback(rank_df, rank, level_name, genome_pool):
    records = []

    for genome_id in genome_pool:
        sub = rank_df[
            (rank_df["genome_id"] == genome_id)
            & (rank_df["participation_index"] >= THRESHOLD)
        ]

        if len(sub) != 1:
            continue

        taxa = sub.iloc[0][rank]

        if taxa.startswith(("o.", "c.")):
            continue

        records.append(
            {
                "genome_id": genome_id,
                "taxa_name": taxa,
                "participation_index": sub.iloc[0]["participation_index"],
                "taxon_level": level_name,
            }
        )

    return pd.DataFrame(
        records,
        columns=["genome_id", "taxa_name", "participation_index", "taxon_level"],
    )


# ==================================================
# PUBLIC API
# ==================================================
def uniref_taxa_fallback(
    diamond_file,
    taxdump,
    prodigal_faa,
    out_final,
    debug_dir=None,
):
    """
    UniRef-based taxonomic semantic fallback (order → class → phylum)

    Parameters
    ----------
    diamond_file : str
        DIAMOND 输出 TSV（protein_id, subject_id, bitscore, taxid）
    taxdump : str
        NCBI taxonomy dump 目录（taxonkit 使用）
    prodigal_faa : str
        Prodigal 预测的 FAA 文件
    out_final : str
        最终输出 TSV
    debug_dir : str or None
        若提供，则输出所有中间 debug 文件
    """
    DEBUG = debug_dir is not None
    if DEBUG:
        os.makedirs(debug_dir, exist_ok=True)

    # ---------- read diamond ----------
    records = []
    taxids = set()

    with open(diamond_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            protein_id = parts[0]
            subject_id = parts[1]
            bitscore = float(parts[2])
            taxid = parts[3]
            genome_id = protein_id.rsplit("_", 1)[0]

            records.append(
                {
                    "protein_id": protein_id,
                    "subject_id": subject_id,
                    "bitscore": bitscore,
                    "taxid": taxid,
                    "genome_id": genome_id,
                }
            )
            taxids.add(taxid)

    df = pd.DataFrame(records)

    # ---------- taxonomy ----------
    taxid2tax = run_taxonkit(taxids, taxdump)
    tax_df = df["taxid"].map(taxid2tax).apply(pd.Series)
    tax_df.columns = ["phylum", "class", "order"]
    df = pd.concat([df, tax_df], axis=1)

    df = fill_na_with_parent(df)
    df = df[df["phylum"] != "NA"]

    if DEBUG:
        df[
            [
                "protein_id",
                "subject_id",
                "bitscore",
                "taxid",
                "phylum",
                "class",
                "order",
            ]
        ].to_csv(
            os.path.join(debug_dir, "diamond_with_lineage.tsv"),
            sep="\t",
            index=False,
        )

    # ---------- ORF ----------
    genome2orf = count_orfs_from_prodigal(prodigal_faa)
    df["num_orf"] = df["genome_id"].map(genome2orf).fillna(0).astype(int)

    # ---------- deduplicate ----------
    df = (
        df.sort_values("bitscore", ascending=False)
        .drop_duplicates(
            subset=[
                "protein_id",
                "phylum",
                "class",
                "order",
                "genome_id",
                "num_orf",
            ],
            keep="first",
        )
    )

    # ---------- rank dfs ----------
    rank_dfs = {}
    for rank in ["order", "class", "phylum"]:
        rdf = finalize_rank_df(build_rank_raw(df, rank))
        rank_dfs[rank] = rdf

        if DEBUG:
            rdf.rename(columns={rank: "taxa_name"}).to_csv(
                os.path.join(debug_dir, f"{rank}.tsv"),
                sep="\t",
                index=False,
            )

    # ---------- semantic fallback ----------
    remaining = set(genome2orf.keys())
    res_df = pd.DataFrame(
        columns=["genome_id", "taxa_name", "participation_index", "taxon_level"]
    )

    for rank in ["order", "class", "phylum"]:
        hit = semantic_fallback(rank_dfs[rank], rank, rank, remaining)
        res_df = pd.concat([res_df, hit], ignore_index=True)
        remaining -= set(hit["genome_id"])

    res_df.to_csv(out_final, sep="\t", index=False)


# ==================================================
# CLI
# ==================================================
def main():
    if len(sys.argv) not in (5, 6):
        sys.exit(
            "Usage:\n"
            "  python uniref_taxa_fallback.py "
            "<diamond.tsv> <taxdump/> <prodigal.faa> <final.tsv> [debug_dir/]"
        )

    diamond_file = sys.argv[1]
    taxdump = sys.argv[2]
    prodigal_faa = sys.argv[3]
    out_final = sys.argv[4]
    debug_dir = sys.argv[5] if len(sys.argv) == 6 else None

    uniref_taxa_fallback(
        diamond_file,
        taxdump,
        prodigal_faa,
        out_final,
        debug_dir=debug_dir,
    )


if __name__ == "__main__":
    main()
    