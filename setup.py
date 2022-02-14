import setuptools
from distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

long_description='ProtView is designed to present statistics of in silico digestions and provide useful information, ' \
                 'such as the protein sequence coverage, peptides covering exon-exon junctions, and the percentage of ' \
                 'junctions or residues in the data that are covered by peptides of a digest. README and instructions ' \
                 'can be found on the [ProtView github page](https://github.com/SSPuliasis/ProtView)'

setup(
    name='protview',
    version='1.0.4',
    packages=['protview'],
    url='https://github.com/SSPuliasis/ProtView',
    license='GPL v3',
    author='Sophia Puliasis',
    author_email='SSPuliasis@dundee.ac.uk',
    description='Digest scheme evaluation in proteomic & proteogenomic context',
    long_description = long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'numpy',
        'pandas',
        'gffpandas==1.2.0',
        'BioPython',
        'rpg==1.1.0'
    ],
    python_requires='>=3.8',
    setup_requires=['wheel'],
    entry_points={
        'console_scripts':[
            'protview=protview.ProtView:main'
        ]
    }
)
