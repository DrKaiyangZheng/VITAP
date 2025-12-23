![logo](/images/logo.png)
# **The VIral TAxonomic Assignment Pipeline (VITAP)**
VITAP (Viral Taxonomic Assignment Pipeline), a high-accuracy pipeline for classifying DNA and RNA viral fragments. VITAP significantly enhances the precision of viral fragment classification by combining alignment-based methods with graph algorithm and providing confidence scores for each taxonomic unit. It can automatically update its alignment database to align with the latest viral reference releases from the International Committee on Taxonomy of Viruses (ICTV), effectively classifying viral fragments as short as 1000 base pairs.
## The workflow of VITAP
![logo](/images/workflow.png)
The VITAP contains two main parts: the update part enables to generate the VITAP-specific database based on current viral taxonomic information from ICTV; the prediction part utilizes previously built database to perform taxonomic assignment. The VITAP utilized an algorithm to determine the best-fit taxonomic lineage, and provides result with various confidence level.
## Undate to VITAP v.1.10
We recognized that a major limitation of relying solely on the ICTV VMR for viral taxonomic assignments is the restricted diversity of protein collections, which negatively affects overall annotation rates and may introduce potential biases in the assignment of certain lineages. To address this issue, we incorporated the external protein database UniRef90 in this version to calibrate and supplement taxonomic assignments.

In the UniRef90-based taxonomy fallback, a single sequence often corresponds to multiple open reading frames (ORFs), and these ORFs may produce BLASTp hits to the same or to different taxa. To quantify the representativeness and support of a given taxon within a genome, we defined a `Participation Index` (PI). The `PI` is used to quantify the representativeness of a taxonomic unit inferred from UniRef90 annotations. For a given genome G with a total of `|P_G|` predicted ORFs and a taxonomic unit T, the `PI` is defined as:

`PI(G, T) = |P_{G,T}| / |P_G|`,

where |P_{G,T}| denotes the number of ORFs in genome G whose UniRef90 hits are assigned to taxon T. The `PI` ranges from 0 to 1 and reflects the proportion of ORFs supporting a given taxonomic assignment.

By default, `PI` is set to 0.25, corresponding to a scenario in which 50% of the ORFs in genome G (referring to geNomad) are assigned to taxon A, and the ratio of the summed bitscore for taxon A to the maximum summed bitscore across all taxa in genome G is 0.5. This harmonization parameter helps avoid misleading taxonomic assignments arising from two scenarios: (1) a single ORF exhibiting an exceptionally high bitscore to taxon A, or (2) multiple ORFs mapping to taxon A but with generally low sequence similarity.

In a study of viromes from Mariana Trench sediments (unpublished), the annotation rate increased by 11% compared with VITAP v1.7.

In addition, we replaced *pandas* with *polaris* in several steps and substituted *prodigal* with *pyrodigal*, resulting in a substantial improvement in computational efficiency. Many thanks to [Valentyn](https://github.com/valentynbez).

Because UniRef90 has been integrated, prebuilt databases are no longer provided.
## Installation
The VITAP is compiled with Python3 and tested on Mac OSX and Linux CentOS, and can run on the above two systems.
### Manually deployed environment
```
conda install -n base -c conda-forge mamba -y
git clone https://github.com/DrKaiyangZheng/VITAP.git
conda install -n base -c conda-forge mamba -y
mamba create -n vitap_v.1.10 -c conda-forge -c bioconda \
  python>=3.9,<3.12 \
  diamond=2.1.16 \
  entrez-direct=16.2 \
  seqkit \
  taxonkit \
  pandas \
  polars \
  pyarrow \
  tqdm \
  biopython \
  pyrodigal \
  -y
```
### Deploy the environment via provided images
```
conda install -n base -c conda-forge mamba -y
git clone https://github.com/DrKaiyangZheng/VITAP.git
cd VITAP
mamba env create -f env.yaml -n vitap_v.1.10
```
### Installation through conda/mamba (Recommanded)
```
mamba create -n vitap_v.1.10 -c conda-forge vitap=1.10
```
### Use the pre-built database
#### Outdated version, not recommended for traditional use! ####
The pre-built database (a hybrid database integrated VMR-MSL37, NCBI RefSeq209, and IMG/VR v.4) can be retrived from [Figshare repository](https://doi.org/10.6084/m9.figshare.25426159.v3). The database is compressed in *.zip* format, with the compressed package being approximately 0.62 GB and about 2.1 GB after decompression (sha256: `13ac7e4790cfdad2a7f0d6c3ebc8a66badfab98138693bb86f7b2a71a91cb951`). You need to place all the folders generated from decompressing the database into the VITAP directory. There should be 3 folders; ignore the MacOSX folder, as I compressed it on my own Mac.

## Usage
### Preparation of database utilized by VITAP
If you need to prepare or update the VITAP database yourself, the following are the standard steps provided:
#### 1. Explore the ICTV website to retrive and download the latest [VMR_MSL file](https://ictv.global/vmr);
#### 2. The downloaded VMR_MSL file is a table with many columns in Excel format. Only **9** columns are required, including **Virus GENBANK accession**,	**Realm**,	**Kingdom**,	**Phylum**,	**Class**,	**Order**,	**Family**,	**Genus**, and	**Species**. You need to extract these columns from original one, prepare another table, and save as *.csv* format.
#### 3. After you have prepared the *.csv* file, use **VITAP upd** to generate or update the database.
```
cd VITAP
./VITAP upd --vmr ./ICTV_VMR_MSL38/ICTV_VMR-MSL38.csv -o ./ICTV_VMR_MSL38/ICTV_VMR-MSL38_reformat.csv -d VMR-MSL
```
This command shall generate a reformat table from a 9-column table, which can be recognized and use by VITAP.
##### **[WARNING]** You need to carefully check if the locus information is correct. Since the original ICTV_VMR_MSL uses Excel for storage, and the start and end loci are separated by a ".", such as “AE006468 (844298.877981)”, Excel might treat it as a floating-point number (omitting the final 0), which could cause errors in the program. Please pay special attention to these categories: Peduoviridae (e.g., AE006468), Belpaoviridae (e.g., LK928904), all GTA-viriform (Bartogtaviriformidae and Rhodogtaviriformidae). The related metadata of these categories is incorrect in ICTV_VMR-MSL38_210426. Please ensure all errors have been corrected and saved.
Once you have verified the table, VITAP will automatically download genomes from GenBank based on genome accessions, and process those endogenous viral genomes that require the extraction of specific loci.
##### **[NOTE]** If you only need to follow the ICTV taxonomic framework, then you just need to wait for the program to finish running; however, if you wish to integrate certain other genomes and their classification information into VITAP, you will need to manually add this information to the *_reformat.csv* file during the confirmation step. Also, add the related genomes to the *VMR_Genome* folder, with each genome stored in its own FASTA file. Let the programe continue after you complete these two steps. You can prepare by referring to the folder *DB_hybrid_MSL37_RefSeq209_IMGVR* in the pre-built database, which integrates some viral genomes from NCBI RefSeq209 and IMG/VR (v.4).
#### 3. VITAP will create a folder prefixed with *DB_* and suffixed with the label you assigned, which serves as the VITAP database (e.g. *DB_NAME*).
##### **[NOTE]** You may choose to delete the *VMR_Genome* folder, as VITAP does not require this folder for performing taxonomic assignments. However, for the convenience of future updates or customization of the VITAP database, I recommend retaining this folder. Each time VITAP runs the update program, it checks against the *_reformat.csv* to see if the *VMR_Genome* folder exists and if any genomes are missing, and it will download those genomes that are present in *_reformat.csv* but not in *VMR_Genome*. This allows VITAP to avoid re-downloading all genomes with each update. This is also why I ask you to manually add taxonomic information to *_reformat.csv* and place those genomes that cannot be downloaded from GenBank into the VMR_Genome folder in advance when creating a custom VITAP database, to prevent any errors in the program.
### Perform taxonomic assignments
#### 1. Performing taxonomic assignment with VITAP is very easy:
```
cd VITAP
./VITAP assignment -i test.fasta -d ./DB_VMR-MSL/ -o test_result
```
#### 2. A folder called *test_result* shall be generated, containing 17 files. The *best_determined_lineages.tsv* should be the one you need.
### The result generated by VITAP
#### 1. The 3 open reading frames related files: *merge_genome.fasta*, *merge_genome.faa*, *merge_genome.gff*
#### 2. The length information of genomes (*merged_genome_length.tsv*).
#### 3. A alignment file generated by DIAMOND against viral reference protein sets (*target_ICTV_blastp.align*);
#### 4. The 8 bipartite network files describe the quantified relationship between genomes and taxonomic hierarchies (from Realm to Species)(*Genus_bipartite.graph*);
#### 5. A reference genome file in FASTA format (*ICTV_selected_genomes.fasta*). The VITAP will ramdomly select several reference genomes and add to your genome file. These reference genomes will be performed taxonomic assignments together with yours;
#### 6. A alignment file generated by DIAMOND against UniRef90 protein sets (*target_uniref90.align*);
#### 7. A UniRef90-based fallback file used to refine and supplement the assignment results;
#### 8. The 2 taxonomic assignment results: *all_lineages.tsv* includes all possible taxonomic lineages; *best_determined_lineages.tsv* includes best-fit taxonomic lineage of a genome;
#### 9. A log file (*VITAP_run.log*).
## Citation
If VITAP contributes to your work, I am grateful if you can cite the [Zheng, K. et al. VITAP: a high precision tool for DNA and RNA viral classification based on meta-omic data. Nat Commun 16, 2226 (2025)](https://www.nature.com/articles/s41467-025-57500-7).


