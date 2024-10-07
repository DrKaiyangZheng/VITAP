from setuptools import setup, find_packages
from pathlib import Path

curr_dir = Path(__file__).parent
long_description = (curr_dir / "README.md").read_text()

# https://docs.python.org/3/distutils/sourcedist.html
# https://github.com/pypa/sampleproject/blob/master/setup.py
setup(
    name='VITAP',
    version='1.7',
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
        'numpy >=1.25',
        'pandas >=1.5',
        'entrez-direct >=16.2',
        'seqkit >=2.5',
        'prodigal >=2.6',
        'setuptools >=65.6',
        'scipy >=1.10.1',
        'networkx >=3.1',
        'tqdm >=4.65.0',
        'biopython >=1.78',
        'diamond >=0.9',
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GPL3 License",
        "Programming Language :: Python :: 3",
    ],
    scripts=['scripts/VITAP',
             'scripts/VITAP_assignment_20241007.py',
             'scripts/VITAP_upd_20241007.py'
    ]
)
