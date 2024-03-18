![logo](/images/logo.png)
# **The VIral TAxonomic Assignment Pipeline (VITAP)**
VITAP (Viral Taxonomic Assignment Pipeline), a high-accuracy pipeline for classifying DNA and RNA viral fragments. VITAP significantly enhances the precision of viral fragment classification by combining alignment-based methods with graph algorithm and providing confidence scores for each taxonomic unit. It can automatically update its alignment database to align with the latest viral reference releases from the International Committee on Taxonomy of Viruses (ICTV), effectively classifying viral fragments as short as 1000 base pairs.
## The workflow of VITAP
![logo](/images/workflow.png)
The VITAP workflow includes two main sections: generation of taxonomic-specific database, and taxonomic assignments for target genomes. The first step enables users to generate a VITAP-specific database according to every release of ICTV proposal. The genomes included in viral metadata resource (VMR) are automatically retrieved and downloaded from GenBank, and used to generate viral reference protein database. To differentiate different taxonomic signal possessed by different viral protein, the self-alignment is performed to determine diverse taxonomic-specific proteins alignment cut-offs. The taxonomic-specific proteins alignment cut-offs were used to filter the proteins possessed various alignment scores. The filtered protein alignment scores were used to calculate the taxonomic units’ cut-offs. Hence, the VITAP-specific database includes viral reference protein database, taxonomic-specific proteins alignment cut-offs, taxonomic units’ cut-offs, and VMR information. The second step is taxonomic prediction of target genomes utilizing VITAP-specific database. The proteins of target genomes are firstly aligned to viral reference proteins. The proteins with taxonomic signals were filtered from the protein alignment between target genomes and viral references, and are used to calculate the taxonomic scores. The taxonomic scores are used in best-fit taxonomic unit determination. Based on their taxonomic scores compared to related cut-offs, this part of results is defined as high-/medium-confidence results. Then, based on previous taxonomic assignments and self-alignment of target genome on protein level, the taxonomic units are extended by utilizing weighted projection. This part of results is defined as low-confidence results.
## Installation
The VITAP is compiled with Python3 and tested on Mac OSX and Linux CentOS, and can run on the above two systems.
### Required Dependencies
####   · biopython  1.78
####   · diamond  0.9
####   · entrez-direct  16.2
####   · networkx  3.1
####   · numpy  1.25
####   · pandas  1.5.3
####   · prodigal  2.6.3
####   · python  3.10.13
####   · scipy  1.10.1
####   · seqkit  2.5.1
####   · tqdm  4.65.0
####   · wget  1.21.4
All these dependences can be installed using Anaconda
### Easy way to configuration environment
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
### Use the pre-built database
You may wish to directly download the database rather prepare by your own before using VITAP. The pre-built database can be retrived from [Figshare repository](https://doi.org/10.6084/m9.figshare.25426159.v1). The database is compressed in *.zip* format, with the compressed package being approximately 1 GB and about 4.5 GB after decompression. The hash value (MD5) is `2dd77e083b1d35569201cf90adb7373f`. You need to place all the folders generated from decompressing the database into the VITAP directory. There should be 3 folders; ignore the MacOSX folder, as I compressed it on my own Mac.
### Install using Anaconda (Recommend)
The VITAP now has been distributed by Anaconda. You can easily install VITAP using conda:
```
conda install -c bioconda VITAP
```
Please note, you need to download the pre-built database or prepare by you own.
## Usage
### Preparation of database utilized by VITAP
If you need to prepare or update the VITAP database yourself, the following are the standard steps provided:
#### 1. Explore the ICTV website to retrive and download the latest [VMR_MSL file](https://ictv.global/msl);
#### 2. The downloaded VMR_MSL file is a table with many columns in Excel format. Only **9** columns are required, including **Virus GENBANK accession**,	**Realm**,	**Kingdom**,	**Phylum**,	**Class**,	**Order**,	**Family**,	**Genus**, and	**Species**. You need to extract these columns from original one, prepare another table, and save as *.csv* format. Please see folder *ICTV_VMR_MSL38* in the pre-built database, which has been deposited at [Figshare repository](https://doi.org/10.6084/m9.figshare.25426159.v1). This folder contains three files, including a original ICTV_VMR-MSL38 file, a ICTV_VMR-MSL38 file in *.csv* format containing **9** columns, and a processed ICTV_VMR-MSL38 file automatically generated by VITAP update program. 
#### 3. After you have prepared the *.csv* file, use **VITAP upd** to generate or update the database.
```
cd VITAP
./scripts/VITAP.sh upd --vmr ./ICTV_VMR_MSL38/ICTV_VMR-MSL38_IMGVR_221222.csv -o ./ICTV_VMR_MSL38/ICTV_VMR-MSL38_IMGVR_221222_reformat.csv
```
This command shall generate a reformat table from a 9-column table, which can be recognized and use by VITAP.
##### **[WARNING]** You need to carefully check if the locus information is correct. Since the original ICTV_VMR_MSL uses Excel for storage, and the start and end loci are separated by a ".", such as “AE006468 (844298.877981)”, Excel might treat it as a floating-point number (omitting the final 0), which could cause errors in the program. Please pay special attention to these categories: Peduoviridae (e.g., AE006468), Belpaoviridae (e.g., LK928904), all GTA-viriform (Bartogtaviriformidae and Rhodogtaviriformidae). The related metadata of these categories is incorrect in ICTV_VMR-MSL38_210426. Please ensure all errors have been corrected and saved.
Once you have verified the table, VITAP will automatically download genomes from GenBank based on genome accessions, and process those endogenous viral genomes that require the extraction of specific loci.
##### **[NOTE]** If you only need to follow the ICTV taxonomic framework, then you just need to wait for the program to finish running; however, if you wish to integrate certain other genomes and their classification information into VITAP, you will need to manually add this information to the *_reformat.csv* file during the confirmation step. Also, add the related genomes to the *VMR_Genome* folder, with each genome stored in its own FASTA file. Let the programe continue after you complete these two steps. You can prepare by referring to the folder *ICTV_VMR_MSL38* in the pre-built database, which integrates some viral genomes from IMG/VR (v.4).
#### 3. VITAP will create a folder prefixed with *DB_* and suffixed with the update date, which serves as the VITAP database (e.g. *DB_20240115*).
##### **[NOTE]** You may choose to delete the *VMR_Genome* folder, as VITAP does not require this folder for performing taxonomic assignments. However, for the convenience of future updates or customization of the VITAP database, I recommend retaining this folder. Each time VITAP runs the update program, it checks against the *_reformat.csv* to see if the *VMR_Genome* folder exists and if any genomes are missing, and it will download those genomes that are present in *_reformat.csv* but not in *VMR_Genome*. This allows VITAP to avoid re-downloading all genomes with each update. This is also why I ask you to manually add taxonomic information to *_reformat.csv* and place those genomes that cannot be downloaded from GenBank into the VMR_Genome folder in advance when creating a custom VITAP database, to prevent any errors in the program.


