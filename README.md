![logo](/images/logo.png)
# **The VIral TAxonomic Assignment Pipeline (VITAP)**
VITAP (Viral Taxonomic Assignment Pipeline), a high-accuracy pipeline for classifying DNA and RNA viral fragments. VITAP significantly enhances the precision of viral fragment classification by combining alignment-based methods with graph algorithm and providing confidence scores for each taxonomic unit. It can automatically update its alignment database to align with the latest viral reference releases from the International Committee on Taxonomy of Viruses (ICTV), effectively classifying viral fragments as short as 1000 base pairs.
## The workflow of VITAP
![logo](/images/workflow.png)
The VITAP contains two main parts: the update part enables to generate the VITAP-specific database based on current viral taxonomic information from ICTV; the prediction part utilizes previously built database to perform taxonomic assignment. The VITAP utilized an algorithm to determine the best-fit taxonomic lineage, and provides result with various confidence level.
## Installation
The VITAP is compiled with Python3 and tested on Mac OSX and Linux CentOS, and can run on the above two systems.
### Required Dependencies
####   · biopython  ≥1.78
####   · diamond  ≥0.9
####   · entrez-direct  ≥16.2
####   · networkx  ≥3.1
####   · numpy  ≥1.25
####   · pandas  ≥1.5.3
####   · prodigal  ≥2.6.3
####   · python  ≥3.9
####   · scipy  ≥1.10
####   · seqkit  ≥2.5
####   · tqdm  ≥4.65
####   · wget  1.21
All these dependences can be installed using Anaconda
### Installation from source code
For MacOS users,
```
git clone https://github.com/DrKaiyangZheng/VITAP.git
cd VITAP
conda env create -f VITAP_conda_environment_OSX.yaml -n VITAP
```
For Linux users,
```
git clone https://github.com/DrKaiyangZheng/VITAP.git
cd VITAP
conda env create -f VITAP_conda_environment_Linux.yaml -n VITAP
```
### Installation through conda (Recommanded)
```
conda install bioconda::vitap
```
### Use the pre-built database
You may wish to directly download the database rather prepare by your own before using VITAP. The pre-built database (a hybrid database integrated VMR-MSL37, NCBI RefSeq209, and IMG/VR v.4) can be retrived from [Figshare repository](https://doi.org/10.6084/m9.figshare.25426159.v2]). The database is compressed in *.zip* format, with the compressed package being approximately 0.62 GB and about 2.1 GB after decompression (sha256: `13ac7e4790cfdad2a7f0d6c3ebc8a66badfab98138693bb86f7b2a71a91cb951`). You need to place all the folders generated from decompressing the database into the VITAP directory. There should be 3 folders; ignore the MacOSX folder, as I compressed it on my own Mac.
## Usage
### Preparation of database utilized by VITAP
If you need to prepare or update the VITAP database yourself, the following are the standard steps provided:
#### 1. Explore the ICTV website to retrive and download the latest [VMR_MSL file](https://ictv.global/msl);
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
#### 2. A alignment file generated by DIAMOND against viral reference protein sets (*target_ICTV_blastp.align*);
#### 3. The 8 bipartite network files describe the quantified relationship between genomes and taxonomic hierarchies (from Realm to Species)(*Genus_bipartite.graph*);
#### 4. A reference genome file in FASTA format (*ICTV_selected_genomes.fasta*). The VITAP will ramdomly select several reference genomes and add to your genome file. These reference genomes will be performed taxonomic assignments together with yours;
#### 5. The 2 taxonomic assignment results: *all_lineages.tsv* includes all possible taxonomic lineages; *best_determined_lineages.tsv* includes best-fit taxonomic lineage of a genome;
#### 6. A log file (*VITAP_run.log*).
## Citation
The manuscript of VITAP has not been published. If VITAP contributes to your work, I am grateful if you can cite the [preprint](https://doi.org/10.21203/rs.3.rs-4406120/v1) or this repository.


