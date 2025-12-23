from setuptools import setup, find_packages
from pathlib import Path

curr_dir = Path(__file__).parent
long_description = (curr_dir / "README.md").read_text()

# https://docs.python.org/3/distutils/sourcedist.html
# https://github.com/pypa/sampleproject/blob/master/setup.py
setup(
    name='VITAP',
    version='1.8',
    packages=find_packages(exclude=['images']),
    description='Viral Taxonomic Assignment Pipeline',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPL-3',
    url='https://github.com/DrKaiyangZheng/VITAP/',
    author='Kaiyang Zheng',
    author_email='zhengkaiyang@stu.ouc.edu.cn',
    keywords=['virus', 'metagenomics', 'metatranscriptomics', 'viral taxonomy'],
    python_requires='>=3.9',
    # https://packaging.python.org/discussions/install-requires-vs-requirements/
    install_requires=[
        'pandas>=1.5',
        'polars>=0.19',
        'pyarrow >=20',
        'tqdm>=4.65.0',
        'biopython>=1.78',
        'pyrodigal>=3.6',
        'setuptools>=65.6',
    ],

    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GPL3 License",
        "Programming Language :: Python :: 3",
    ],
    scripts=['scripts/VITAP',
             'scripts/VITAP_assignment.py',
             'scripts/VITAP_upd.py',
             'scripts/uniref_taxa_fallback.py',
             'scripts/uniref90_accession2taxid.py'
    ]
)
