[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
# Table of Contents

* [Overview](#overview)
* [Installation](#installation)
* [Input files](#input-files)
* [Walk-through](#walk-through)
    * [RPG digest](#rpg-digest)
    * [Peptide Processing](#peptide-processing)
    * [Coding Sequence Extraction](#coding-sequence-extraction)
    * [General Summary Statistics](#general-summary-statistics)
    * [Genomic co-ordinate conversion](#genomic-co-ordinate-conversion)
    * [Exon-exon Junction Covering Peptide Identification & Summary Statistics](#junction-covering-peptides--statistics)
    * [Unique Peptide Counter](#unique-peptide-counter)
* [Contact](#contact)

# Overview
ProtView is designed to present statistics of in silico digestions and provide useful information, 
such as the protein sequence coverage, peptides covering exon-exon junctions, and the percentage 
of junctions or residues in the data that are covered by peptides of a digest (Figure 1). It allows 
the users to see the portions of the gene/transcript on the genome that are covered by proteomic 
data. It takes the protein sequence (fasta) and the coding sequence annotations (gff3 format) on 
the genome as inputs. It incorporates Rapid Peptides Generator (RPG) 
([Maillet, 2019](https://academic.oup.com/nargab/article/2/1/lqz004/5581718)), which carries out 
the in-silico digestion. It then maps the digested proteins back to transcripts/genes on the genome, 
which allows the comparisons of the transcript/ gene sequences visible to proteomics experiments 
between different digestions.

![alt text](./README_images/protview_context.png)

**Figure 1:** outline of the ProtView workflow

# Installation
Supported Python version: 3.7.4

ProtView requires the following packages:
* [rpg](https://rapid-peptide-generator.readthedocs.io/en/latest/userguide.html#installation) version 1.1.0
* [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html?highlight=install) version 0.25.1
* [gffpandas](https://gffpandas.readthedocs.io/en/latest/installation.html) version 1.2.0
* [numpy](https://numpy.org/install/) version 1.16.5
* [BioPython](https://biopython.org/wiki/Download) version 1.76

To install the dependencies please follow the instructions in their documentation.

The ProtView repository can be downloaded as a zip file or cloned. To clone, use the command:
`git clone https://github.com/SSPuliasis/ProtView.git`

# Input Files
Required input for an analysis with ProtView:
* Protein sequence(s) in fasta format
* Corresponding gff3 file(s)

Any additional input required for downstream functions is generated by ProtView during the analysis.

# Walk-through
To view the available ProtView modules use `python path\to\ProtView\protview.py -h`. To get help for an
individual module use `python path\to\ProtView\protview.py <module name> -h`

The *Arabidopsis thaliana* genes AT1G666600 and AT1G66610 (Figure 2) are included as example data to 
walk through how ProtView works. Araport11 protein sequences were retrieved from the 
[tair database](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FSequences%2FAraport11_blastsets) 
in fasta format and the Araport11 GFF3 files were downloaded via [jbrowse](https://www.araport.org/data/araport11).
These files can be found in the example_data folder.

![alt text](./README_images/jbrowse_at1g66600_at1g66610.PNG)

**Figure 2:** Jbrowse depiction of the example

## RPG digest
ProtView incorporates Rapid Peptides Generator (RPG) ([Maillet, 2019](https://academic.oup.com/nargab/article/2/1/lqz004/5581718))
to carry out *in-silico* digests. RPG takes a fasta protein sequence as input, more instructions for 
running RPG can be found [here](https://rapid-peptide-generator.readthedocs.io/en/latest/).

RPG can carry out digests in sequential or concurrent mode, where a sequence is cleaved by multiple enzymes. 
ProtView can then combine peptides from sequential digests, creating parallel enzyme combinations.

To digest the example proteins with Trypsin, Asp-N, and Glu-C and save the output as 'at1g6600_at1g66610_rpg.fasta':
`rpg -i example_data\at1g66600_at1g66610.fasta -o example_data\at1g66600_at1g66610_rpg.fasta -e 2 23 42 -q`

## Peptide Processing
Processes RPG-generated peptides in fasta format, allowing the user to specify the number of missed cleavages 
allowed per peptide, create parallel enzyme digests, and filter for peptides containing a specific residue
or of a certain amino acid length.

Required Arguments:
* **input_file**: name of the RPG peptides file

Options:
* **-mc, --miscleavage**:  Mis-cleavage value, default = 0
* **-e, --enzymes**: Enzymes to create a parallel digest with, default = none
* **-r, --residue**: Residue to be filtered for, default = none
* **-min, --min_len**: Minimum peptide length to filter for, default = 7
* **-max, --max_len**: Maximum peptide length to filter for, default = 35

The options are carried out in the order in which they are given here and the output is saved after
each step.

To process the RPG peptides generated above, allowing for no missed cleavages, creating a parallel digest
of Trypsin and Asp-N, and then filter for the default 7-35 amino acid length range:
`python path\to\ProtView\protview.py rpg_output_processing at1g66600_at1g66610_rpg.fasta -e Trypsin Asp-N`

To filter the individual enzyme digests, it must be ran again without the parallel digest option.
`python path\to\ProtView\protview.py rpg_output_processing at1g66600_at1g66610_rpg.fasta`

**Note**:
* When creating a parallel digest, the enzymes need to entered exactly as they are in the RPG output (e.g. 'Asp-N'
 and not 'aspn', asp-n', 'Asp-n' etc.). ProtView will not raise an error if an incorrect enzyme name is given, however 
 the output will be inaccurate and not represent the parallel digest.

## Coding Sequence Extraction
Extracts coding sequence (CDS) information from gff3 files and prepares this information for the 
downstream analysis. The resulting tables retain relevant information from the original gff3 file, with 
the addition of relative protein sequence coordinates of the CDS start and end positions, unique CDS IDs, 
and the lengths and unique IDs of introns in between CDSs. Relative protein sequence coordinates are 
calculated and CDS and intron IDs consist of chromosome, start position, end position, and strand.

Required Arguments:
* **input_file**: input file name of gff3 file

To extract CDS information from the gff3 file for the example Arabidopsis proteins:

`python path\to\ProtView\protview.py cds_extraction at1g66600_at1g66610.gff3`

The resulting CDS information is saved separately for each strand in csv format 
('at1g66600_at1g66610_+\_cdsdf.csv, at1g66600_at1g66610\_-_cdsdf.csv,')

**Note**: 
- ProtView will generate a CDS information file for each DNA strand, regardless of whether the input only 
contains protein(s) on one of the strands. It is important to keep empty CDS information files generated, as they 
are required as input for steps in the downstream analysis.

## General Summary Statistics
This gives a table with columns for the total number of peptides generated, both before and after 
filtering criteria were applied, peptide length distributions, and protein sequence coverage. The 
calculation of residue coverage is optional and carried out if residues are provided by the user.

Required Arguments:
* **-fasta, --fasta_files**: original fasta sequence files used in the digests
* **-u, --unfiltered_rpg_files**: unfiltered rpg peptides in csv format
* **-f, --filtered_rpg_files**: filtered rpg peptides in csv format

Optional Arguments:
* **-o, --output_name**: name of the output file, ending in .csv. Default = summary_table.csv
* **-r, --residues**: residues to calculate coverage for

To generate summary statistics of the example data for the individual and parallel digests, including
Cysteine (C), Serine (S) and Threonine (T) coverage:
```
python path\to\ProtView\protview.py summary_stats -fasta at1g66600_at1g66610.fasta -u at1g66600_at1g66610_rpg.csv 
at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel.csv -f at1g66600_at1g66610_rpg_len_7_35.csv 
at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv -r C S T
```

## Genomic co-ordinate conversion
Relative peptide coordinates from the digest output can converted to the outer bounds of their 
corresponding coordinates on the genome. The resulting table contains the isoform, both genomic and 
relative protein start and end positions for each peptide, and the enzymes used to generate the peptide. 

This function works best on one protein at a time and the output allows for genomic coverage visualization using
tools such as [Gviz](https://link.springer.com/protocol/10.1007%2F978-1-4939-3578-9_16).

Arguments:
* **rpg_file**: csv file containing the peptides and their proteomic co-orindates
* **cds_files**: both csv files containing extracted coding sequences for each DNA strand, ending in 
+/-_cdsdf.csv

To calculate, the genomic coordinates of the filtered parallel Trypsin:Asp-N digest in our example:
```
python path\to\ProtView\protview.py gen_coords at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv 
at1g66600_at1g66610_+_cdsdf.csv at1g66600_at1g66610_-_cdsdf.csv
```
**Note**:
 - The purpose of this function is to enable peptide visualisation on the genome and it works best when ran on
 small sets of input

## Junction-covering peptides & statistics
Peptides are filtered for those that cover eon-exon junctions. The recommended input for identifying 
junction-covering peptides is the filtered digest results. Positive outcomes are saved in the same csv 
format as the digest results, with an additional column for junction location.

Arguments:
* **rpg_file**: csv file containing the peptides to be filtered
* **cds_files**: both csv files containing extracted coding sequences for each DNA strand, ending in 
+/-_cdsdf.csv
* **output_name**: desired name of output csv file, ProtView will automatically differentiate between 
DNA strands

To filter the individual digests from the example
```
python path\to\ProtView\protview.py junction_peptides at1g66600_at1g66610_rpg_len_7_35.csv 
at1g66600_at1g66610_+_cdsdf.csv at1g66600_at1g66610_-_cdsdf.csv single_digest_junction_peptides.csv
```
To filter the parallel digest, first argument changes to the csv containing peptides from this digest and another 
output name is provided
```
python path\to\ProtView\protview.py junction_peptides at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv 
at1g66600_at1g66610_+_cdsdf.csv at1g66600_at1g66610_-_cdsdf.csv parallel_digest_junction_peptides.csv
```
By giving 'single_digest_junction_peptides.csv' as the output argument, ProtView adds '\_+\_' or '\_-\_' to the 
end to distinguish between strands, resulting in files named 'single_digest_junction_peptides_+\_.csv' and 
'single_digest_junction_peptides_-_.csv'

A summary table can be generated from the junction-covering peptide results. Each table includes the number of 
junction-covering peptides generated by each enzyme, the number of unique junctions that an enzyme covers (to avoid 
double counting of splice junctions shared between transcripts), and a junction coverage percentage, which is the 
percentage of the total junctions available in the isoforms being examined that are covered by an enzyme. 

Arguments:
* **-pept, --junction_spanning_peptides**: csv files containing the junction covering peptides
* **-cds, --cds_files**: both csv files containing extracted coding sequences for each DNA strand, 
ending in +/-_cdsdf.csv
* **-out, --output_name**: desired name of output csv summary file, default: junction_summary.csv

To summarise the junction-covering peptides identified above:
```
python path\to\ProtView\protview.py junction_summary -pept single_digest_junction_peptides_+_.csv 
single_digest_junction_peptides_-_.csv parallel_digest_junction_peptides_+_.csv parallel_digest_junction_peptides_-_.csv 
-cds at1g66600_at1g66610_+_cdsdf.csv at1g66600_at1g66610_-_cdsdf.csv
```
## Unique Peptide Counter
ProtView can count the number of unique peptides per digest and append the calculation as a column to either the
general proteomic summary table or the junction-covering peptide summary table. Unique peptides are defined as those 
with sequences that can only be found in one isoform in a digest.

Arguments:
* **table_name**: the proteomic or junction summary statistics table for the column to be appended to
* **rpg_files**: the filtered digest results in csv format

To calculate the  number of unique peptides in the example digest and append it to the proteomic summary table:
```
python path\to\ProtView\protview.py unique_count summary_table.csv at1g66600_at1g66610_rpg_len_7_35.csv 
at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv
```
# Contact
For further information please contact <SSPuliasi@dundee.ac.uk>