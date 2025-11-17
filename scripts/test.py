from VITAP_assignment_20241007 import taxonomy_assigning
from polars_utils import taxonomy_assigning as taxonomy_polars
import os
import tempfile
import shutil
import pandas as pd 
# compare two dataframes
from pandas.testing import assert_frame_equal

import time

# --------------- TEST # 1: Taxonomy assignment step --------------- #
print("[INFO] Testing taxonomy assignment step")
db_folder = "./test/DB_MSL-demo/"
target_blast_fp = "./test/toy_set_annotation_result/target_ICTV_blastp.align" 
VMR_csv_file = "./test/DB_MSL-demo/VMR_taxonomy_map_MSL-demo.csv"
length_file = "./test/toy_set_annotation_result/merged_genome_length.tsv"
orf_path = "./test/toy_set_annotation_result/merge_genome.gff"
# create tmp folder for result
result_folder = tempfile.mkdtemp(prefix="VITAP_test_")

taxa = "Species"
genome_cutoff_file = os.path.join(db_folder, f"{taxa}_genome.threshold")
taxa_annot_file = os.path.join(result_folder, f"{taxa}_bipartite.graph")
taxa_annot_file_polars = os.path.join(result_folder, f"{taxa}_bipartite_polars.graph")
blast_graph_file = os.path.join(result_folder, "genome2genome.nwk")
time1 = time.time()
taxonomy_assigning(target_blast_fp, VMR_csv_file, genome_cutoff_file, taxa, length_file, orf_path, taxa_annot_file, taxa_annot_file)
time2 = time.time()
print(f"[INFO] Time taken for pandas version: {time2 - time1} seconds")
time3 = time.time()
taxonomy_polars(target_blast_fp, VMR_csv_file, genome_cutoff_file, taxa, length_file, orf_path, taxa_annot_file_polars)
time4 = time.time()
print(f"[INFO] Time taken for polars version: {time4 - time3} seconds")
print(f"[INFO] Speedup: {(time2 - time1)/(time4 - time3)}x")

# compare the two output files
# load and compare two files
species_graph_pandas = pd.read_csv(taxa_annot_file, sep="\t").sort_values(["Genome_ID", "Species_weight"])
species_graph_polars = pd.read_csv(taxa_annot_file_polars, sep="\t").sort_values(["Genome_ID", "Species_weight"])

# compare two dataframe
assert_frame_equal(species_graph_pandas.reset_index(drop=True), 
                   species_graph_polars.reset_index(drop=True),
                   check_dtype=False)
print("[INFO] The outputs from pandas and polars versions are identical")


# cleanup
shutil.rmtree(result_folder)
