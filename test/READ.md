##Runing the test of VITAP to make sure it works on your system
###When you run the following command and enter "Y":
```
chmod +x ./scripts/*
./scripts/VITAP upd --vmr ICTV_VMR_MSL-demo.csv -o ICTV_VMR_MSL-demo_reformat.csv -d MSL-demo
```
###You shall see the following output (which should be directly displayed in your terminal console):
```
Downloading genomes: 100%|██████████████████████████████████████████████████████████████████████████| 73/73 [01:10<00:00,  1.03it/s]
[INFO] All files successfully downloaded and processed
[INFO] Merging file to VMR_genome_MSL-demo.fasta...
[INFO] Length statitic for VMR_genome_MSL-demo.fasta
[INFO] ORF calling for VMR_genome_MSL-demo.fasta
[INFO] Processing short sequences ignored by prodigal.
[INFO] Genome and ORF files have consistent on non-redundant FASTA IDs. √
[INFO] Cleaning the GFF file
[INFO] Statistic of the number of ORF per genome
[INFO] Moving ICTV_VMR_MSL-demo_reformat.csv to DB_MSL-demo as DB_MSL-demo/VMR_taxonomy_map_MSL-demo.csv
[INFO] Building ICTV reference protein databse.
diamond v0.9.14.115 | by Benjamin Buchfink <buchfink@gmail.com>
Licensed under the GNU AGPL <https://www.gnu.org/licenses/agpl.txt>
Check http://github.com/bbuchfink/diamond for updates.

#CPU threads: 12
Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)
Database file: DB_MSL-demo/VMR_genome_MSL-demo.faa
Opening the database file...  [9.9e-05s]
Loading sequences...  [0.009135s]
Masking sequences...  [0.047632s]
Writing sequences...  [0.004336s]
Loading sequences...  [4e-06s]
Writing trailer...  [2.3e-05s]
Closing the input file...  [5e-06s]
Closing the database file...  [0.00111s]
Processed 7868 sequences, 2409950 letters.
Total time = 0.062398s
[INFO] Self-aligning of ICTV reference proteins.
[INFO] Calculating best-fit taxonomic threshold for Species
[INFO] Calculating best-fit taxonomic threshold for Genus
[INFO] Calculating best-fit taxonomic threshold for Family
[INFO] Calculating best-fit taxonomic threshold for Order
[INFO] Calculating best-fit taxonomic threshold for Class
[INFO] Calculating best-fit taxonomic threshold for Phylum
[INFO] Calculating best-fit taxonomic threshold for Kingdom
[INFO] Calculating best-fit taxonomic threshold for Realm
[INFO] All updating steps finished.
```
###At this point, two folders have been generated: "VMR_Genome/" and "DB_MSL-demo/", the latter being the database automatically generated and used by VITAP.
###Then, when you run the following command:
```
./scripts/VITAP assignment -i toy_set.fasta -d DB_MSL-demo -o toy_set_annotation_result -p 6
```
###You shall see the following output:
```
===== The VITAP (Viral Taxonomic Assignment Pipeline, main program) v.1.5, Copyright 2024 Kaiyang Zheng, Viral/microbial diversity Lab. Ocean University of China ===== 
[INFO] A number of temporary files will be generated in your result directory, please keep them or delete them followed by the guidlines.
===== ORF calling for toy_set.fasta =====
[INFO] Ten ICTV reference genomes were randomly selected and merged to your genome set! 
[INFO] Genome and ORF files have consistent non-redundant FASTA IDs. √
[INFO] Length statitic for merged_genome.fasta
[INFO] Statistic of the number of ORF per genome
===== Aligning to ICTV reference database based on BLAST-algorithm =====
[INFO] Please delete the output file (toy_set_annotation_result/target_ICTV_blastp.align) if this step aborted!
[INFO] Please keep the output file (toy_set_annotation_result/target_ICTV_blastp.align) if this step finished!
===== Locating contig on taxonomic framework and generating genome community network =====
[INFO] Generating Species graph for target genomes
[INFO] Generating Genus graph for target genomes
[INFO] Generating Family graph for target genomes
[INFO] Generating Order graph for target genomes
[INFO] Generating Class graph for target genomes
[INFO] Generating Phylum graph for target genomes
[INFO] Generating Kingdom graph for target genomes
[INFO] Generating Realm graph for target genomes
===== Determining viral lineages based on multi-partite graph ===== 
[INFO] Importing taxonomic bipartite graphs and ICTV taxonomic hierarchy
[INFO] Calculating average taxonomic score for target genomes: 
  -> finding the lowest effective taxonomic level
  -> determining the optimal taxonomic hierarchy termination point
[INFO] Determine best-fit taxonomic hierarchy
[INFO] Assigning confidence level to taxonomic hierarchy
===== Exporting results ===== 
[INFO] Time-consuming (wall clock): 0.0 hours
[INFO] Time-consuming (CPU): 0.0 hours
===== All done ===== 
```
###One folders shall be generated: "toy_set_annotation_result/". At this point, you should have successfully reproduced all the result files in the test/ directory, and VITAP should be working properly. Have fun!
