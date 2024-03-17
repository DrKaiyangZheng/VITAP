![logo](/images/logo.png)
# **The VIral TAxonomic Assignment Pipeline (VITAP)**
VITAP (Viral Taxonomic Assignment Pipeline), a high-accuracy pipeline for classifying DNA and RNA viral fragments. VITAP significantly enhances the precision of viral fragment classification by combining alignment-based methods with graph algorithm and providing confidence scores for each taxonomic unit. It can automatically update its alignment database to align with the latest viral reference releases from the International Committee on Taxonomy of Viruses (ICTV), effectively classifying viral fragments as short as 1000 base pairs.
## Installation
The VITAP is compiled with Python3 and tested on Mac OSX and Linux CentOS, and can run on the above two systems.
### Required Dependencies
Name               Version   
biopython          1.78
blas               1.0
bottleneck         1.3.5
brotli-python      1.0.9
bzip2              1.0.8
ca-certificates    2023.11.17
certifi            2023.11.17
charset-normalizer 3.3.2
diamond            0.9.14
entrez-direct      16.2
idna               3.6
intel-openmp       2023.1.0
libblas            3.9.0
libcblas           3.9.0
libcxx             14.0.6
libffi             3.4.4
libgfortran        5.0.0
libgfortran5       13.2.0
libidn2            2.3.4
liblapack          3.9.0
libunistring       0.9.10
llvm-openmp        17.0.6
mkl                2023.1.0
mkl-service        2.4.0
mkl_fft            1.3.8
mkl_random         1.2.4
ncurses            6.4
networkx           3.1
numexpr            2.8.4
numpy              1.25.2
numpy-base         1.25.2
openssl            3.2.0
packaging          23.2
pandas             1.5.3
pip                23.2.1
platformdirs       4.1.0
pooch              1.8.0
prodigal           2.6.3
pysocks            1.7.1
python             3.10.13
python-dateutil    2.8.2
python_abi         3.10
pytz               2023.3.post1
readline           8.2
requests           2.31.0
scipy              1.10.1
seqkit             2.5.1
setuptools         68.0.0
six                1.16.0
sqlite             3.41.2
tbb                2021.8.0
tk                 8.6.12
tqdm               4.65.0
tzdata             2023c
urllib3            2.1.0
wget               1.21.4
wheel              0.38.4
xz                 5.4.2
zlib               1.2.13

All these dependences can be installed using Anaconda
