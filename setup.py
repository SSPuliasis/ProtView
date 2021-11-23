import setuptools
from distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

long_description='README and instructions can be found on the [ProtView github page](https://github.com/SSPuliasis/ProtView)'

setup(
    name='protview',
    version='0.0.1',
    packages=['modules'],
    url='https://github.com/SSPuliasis/ProtView',
    license='GPL v3',
    author='Sophia Puliasis',
    author_email='SSPuliasis@dundee.ac.uk',
    description='Digest scheme evaluation in proteomic & proteogenomic context',
    long_description = long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'numpy==1.16.5',
        'pandas==0.25.1',
        'gffpandas==1.2.0',
        'BioPython==1.76',
        'rpg==1.1.0'
    ],
    python_requires='>3.7',
    setup_requires=['wheel']
)
