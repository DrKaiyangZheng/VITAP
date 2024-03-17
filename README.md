![logo](/images/logo.png)
# **The VIral TAxonomic Assignment Pipeline (VITAP)**
VITAP (Viral Taxonomic Assignment Pipeline), a high-accuracy pipeline for classifying DNA and RNA viral fragments. VITAP significantly enhances the precision of viral fragment classification by combining alignment-based methods with graph algorithm and providing confidence scores for each taxonomic unit. It can automatically update its alignment database to align with the latest viral reference releases from the International Committee on Taxonomy of Viruses (ICTV), effectively classifying viral fragments as short as 1000 base pairs.
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
```
conda env create -f VITAP_environment.yaml -n VITAP
```
You need to download the database before using VITAP
```
wget http://www.XXXXXX.org
```

