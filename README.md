[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
# Table of Contents

* [Overview](#overview)
* [Installation](#installation)
* [Input files](#input-files)
* [Walk-through](#walk-through)
    * [Input Sanitation](#input-sanitation)
    * [RPG digest](#rpg-digest)
    * [Peptide Processing](#peptide-processing)
    * [Coding Sequence Extraction](#coding-sequence-extraction)
    * [General Summary Statistics](#general-summary-statistics)
    * [Genomic co-ordinate conversion](#genomic-co-ordinate-conversion)
    * [Exon-exon Junction Covering Peptide Identification & Summary Statistics](#junction-covering-peptides--statistics)
    * [Unique Peptide Counter](#unique-peptide-counter)
* [Contact](#contact)

# Overview
ProtView is designed to present statistics of *in silico* digestions and provide useful information, 
such as the protein sequence coverage, peptides covering exon-exon junctions, and the percentage 
of junctions or residues in the data that are covered by peptides of a digest (Figure 1). It allows 
the users to see the portions of the gene/transcript on the genome that are covered by proteomic 
data. It takes the protein sequence (fasta) and the coding sequence annotations (gff3 format) on 
the genome as inputs. It incorporates Rapid Peptides Generator (RPG) 
([Maillet, 2019](https://academic.oup.com/nargab/article/2/1/lqz004/5581718)), which carries out 
the in-silico digestion. It then maps the digested proteins back to transcripts/genes on the genome, 
which allows the comparisons of the transcript/ gene sequences visible to proteomics experiments 
between different digestions.

![diagram of workflow in proteomic and proteogenomic context](https://github.com/SSPuliasis/ProtView/blob/master/README_images/protview_context.png)

**Figure 1:** outline of the ProtView workflow

# Installation
Supported Python version: 3.7.4

ProtView requires the following packages:
* [rpg](https://rapid-peptide-generator.readthedocs.io/en/latest/userguide.html#installation) version 1.1.0
* [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html?highlight=install) version 0.25.1
* [gffpandas](https://gffpandas.readthedocs.io/en/latest/installation.html) version 1.2.0
* [numpy](https://numpy.org/install/) version 1.16.5
* [BioPython](https://biopython.org/wiki/Download) version 1.76

We recommend creating a separate environment to run ProtView in the command line, to avoid clashes between package versions in the main
environment. To create and activate a separate environment named 'protview_environment', and install ProtView and 
dependencies:
```
conda create -n protview_environment
conda activate protview_environment
pip install protview
```

# Input Files
Required input for an analysis with ProtView:
* Protein sequence(s) in fasta format
* Corresponding gff3 file(s)

Any additional input required for downstream functions is generated by ProtView during the analysis.

# Walk-through
ProtView is intended to be used in the command line. To view the available ProtView modules use `protview -h`. To get 
help for an individual module use `protview <module name> -h`

The *Arabidopsis thaliana* genes AT1G666600 and AT1G66610 (Figure 2) are included as example data to 
walk through how ProtView works. Araport11 protein sequences were retrieved from the 
[tair database](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FSequences%2FAraport11_blastsets) 
in fasta format and the Araport11 GFF3 files were downloaded via [jbrowse](https://www.araport.org/data/araport11).
These files can be found and downloaded from the protview_example_data folder on the github page. Once the working directory has 
been set to the directory containing the intended input files, example data or otherwise, the analysis can begin.

![Jbrowse depiction of the example](https://github.com/SSPuliasis/ProtView/blob/master/README_images/jbrowse_at1g66600_at1g66610.PNG)

**Figure 2:** Jbrowse depiction of the example

## Input Sanitation
This step replaces underscores and '>' characters in the fasta protein descriptions to avoid downstream errors in the 
analysis. To run the fasta input sanitation function of the fasta file in the example data:

`protview fasta_input_sanitation at1g66600_at1g66610.fasta`

## RPG digest
ProtView incorporates Rapid Peptides Generator (RPG) ([Maillet, 2019](https://academic.oup.com/nargab/article/2/1/lqz004/5581718))
to carry out *in-silico* digests. RPG takes a fasta protein sequence as input, more instructions for 
running RPG can be found [here](https://rapid-peptide-generator.readthedocs.io/en/latest/).

RPG can carry out digests in sequential (default) or concurrent mode (`-d c`), where a sequence is cleaved by multiple enzymes. 
ProtView can then combine peptides from sequential digests, creating parallel enzyme combinations.

To digest the example proteins with Trypsin, Asp-N, and Glu-C and save the output as 'at1g6600_at1g66610_rpg.fasta':
`rpg -i at1g66600_at1g66610.fasta -o at1g66600_at1g66610_rpg.fasta -e 2 23 42 -q`. The argument
`-d c` can be added if carrying out a concurrent digest with the enzymes, rather than sequential, and concurrent
digests are then treated the same as single enzyme digests throughout the workflow.

## Peptide Processing
Processes RPG-generated peptides in fasta format, allowing the user to specify the number of missed cleavages 
allowed per peptide, create parallel enzyme digests, and filter for peptides containing a specific residue
or of a certain amino acid length.

RPG treats mis-cleavage as a percentage of how frequently cleavage is missed at each theoretical cleavage site. The 
aim of ProtView is to give the theoretical upper limit of peptides that could be identified in an experiment and allow 
for *in silico* results to be compared to experimental data from software algorithms that treat mis-cleavage as the number 
of missed cleavages allowed per peptide . A function for generating mis-cleaved peptides is therefore included in 
ProtView, which allows the generation of mis-cleaved peptides by concatenating adjacent peptide sequences from the RPG 
digest up to n times, where n is the mis-cleavage number specified by the user.

Required Arguments:
* **input_file**: name of the RPG peptides file

Options:
* **-mc, --miscleavage**:  Mis-cleavage value, default = 0
* **-e, --enzymes**: Enzymes to create a parallel digest with, default = none
* **-r, --residue**: Residue to be filtered for, default = none
* **-min, --min_len**: Minimum peptide length to filter for, default = 7
* **-max, --max_len**: Maximum peptide length to filter for, default = 35

The options are carried out in the order in which they are given here and the output is saved after each step.

Output files are named \<input_file>_\<parameter_suffixes>.csv
Parameter suffixes:
* mc\<x>, where x = the number of missed cleavages, if mis-cleavage is used
* \<enzymes>_parallel, if a parallel digest is created
* len_\<min length>_\<max length>, if filtered for amino acid length
* \<residue> , if filtered for a specific residue

Output columns:
* **FASTA_description**: Description of the protein from th original FASTA file
* **cleavage_position**: Protein end coordinate of the peptide
* **peptide_size**: peptide length in amino acids
* **mol_weight**: molecular weight
* **isoelectric point**
* **sequence**
* **peptide_start**: Protein start coordinate of the peptide

To process the RPG peptides generated above, allowing for no missed cleavages, creating a parallel digest
of Trypsin and Asp-N, and then filter for the default 7-35 amino acid length range:
`protview rpg_output_processing at1g66600_at1g66610_rpg.fasta -e Trypsin Asp-N`

To filter the individual enzyme digests, it must be ran again without the parallel digest option.
`protview rpg_output_processing at1g66600_at1g66610_rpg.fasta`

**Notes**:
* When creating a parallel digest, the enzymes need to entered exactly as they are in the RPG output (e.g. 'Asp-N'
 and not 'aspn', asp-n', 'Asp-n' etc.)
* Any files containing peptides from parallel digests or containing missed cleavages will have 'parallel' or 'mc' in 
their file names. It is important not to rename these files, as some of the statistics are calculated differently for 
mis-cleaved peptides and combined digests.

## Coding Sequence Extraction
Extracts coding sequence (CDS) information from gff3 files and prepares this information for the 
downstream analysis. The resulting tables retain relevant information from the original gff3 file, with 
the addition of relative protein sequence coordinates of the CDS start and end positions, unique CDS IDs, 
and the lengths and unique IDs of introns in between CDSs. Relative protein sequence coordinates are 
calculated and CDS and intron IDs consist of chromosome, start position, end position, and strand.

Required Arguments:
* **input_file**: input file name of gff3 file

To extract CDS information from the gff3 file for the example Arabidopsis proteins:

`protview cds_extraction at1g66600_at1g66610.gff3`

The resulting CDS information is saved separately for each strand in csv format 
('at1g66600_at1g66610_+\_cdsdf.csv, at1g66600_at1g66610\_-_cdsdf.csv,')

Output columns:
Information is retained from the original gff3 files under the same columns names. Columns added by ProtView are:
* **cds_id**: unique ID of the CDS
* **intron_id**: unique ID of the adjacent intron before the CDS 
* **intron_length**: length of the adjacent intron
* **cum_intron**: cumulative length of introns before the CDS
* **cds_first_start**: start position of the first CDS of an isoform
* **protein_start**: protein start coordinate of the CDS
* **protein_end**: protein end coordinate of the CDS

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

Output columns:
* **enzyme**: the protease or protease combination used in the digest
* **total_peptides**: the total number of peptides generated by that enzyme
* **mean_length**: mean length of all peptides generated by this enzyme
* **median_length**: median length of all peptides generated by this enzyme
* **filtered_peptides**: number of peptides that meet filtering criteria
* **coverage**: coverage % of the original digested protein sequence

To generate summary statistics of the example data for the individual and parallel digests, including
Cysteine (C), Serine (S) and Threonine (T) coverage:
```
protview summary_stats -fasta at1g66600_at1g66610.fasta -u at1g66600_at1g66610_rpg.csv 
at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel.csv -f at1g66600_at1g66610_rpg_len_7_35.csv 
at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv -r C S T
```

**Note:**
- Please ensure that any files containing peptides with missed cleavages or parallel digests have 'mc' or 'parallel'
in their file names (standard ProtView output unless the names have been manually changed) before carrying out this step

## Genomic co-ordinate conversion
Relative peptide coordinates from the digest output can converted to the outer bounds of their 
corresponding coordinates on the genome. The resulting table contains the isoform, both genomic and 
relative protein start and end positions for each peptide, and the enzymes used to generate the peptide. 

This function works best on one protein at a time and the output allows for genomic coverage visualization using
tools such as [Gviz](https://link.springer.com/protocol/10.1007%2F978-1-4939-3578-9_16).

Required Arguments:
* **rpg_file**: csv file containing the peptides and their proteomic co-orindates
* **cds_files**: both csv files containing extracted coding sequences for each DNA strand, ending in 
+/-_cdsdf.csv

To calculate, the genomic coordinates of the filtered parallel Trypsin:Asp-N digest in our example:
```
protview gen_coords at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv 
at1g66600_at1g66610_+_cdsdf.csv at1g66600_at1g66610_-_cdsdf.csv
```

Output columns:
* **isoform**: the isoform that the peptide originates from
* **Protein start coordinate**: relative start coordinate on the protein sequence
* **Genomic start coordinate**: calculated start coordinate on the genome
* **Protein end coordinate**: relative end coordinate on the protein sequence
* **Genomic end coordinate**: calculated end coordinate on the genome
* **enzyme**: the protease or protease combination used in the digest

**Note**:
 - The purpose of this function is to enable peptide visualisation on the genome and it works best when ran on
 small sets of input

## Junction-covering peptides & statistics
Peptides are filtered for those that cover exon-exon junctions. The recommended input for identifying 
junction-covering peptides is the filtered digest results. Positive outcomes are saved in the same csv 
format as the digest results, with an additional column for junction location.

Required Arguments:
* **rpg_file**: csv file containing the peptides to be filtered
* **cds_files**: both csv files containing extracted coding sequences for each DNA strand, ending in 
+/-_cdsdf.csv
* **output_name**: desired name of output csv file, ProtView will automatically differentiate between 
DNA strands

To filter the individual digests from the example
```
protview junction_peptides at1g66600_at1g66610_rpg_len_7_35.csv at1g66600_at1g66610_+_cdsdf.csv 
at1g66600_at1g66610_-_cdsdf.csv single_digest_junction_peptides.csv
```
To filter the parallel digest, first argument changes to the csv containing peptides from this digest and another 
output name is provided
```
protview junction_peptides at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv 
at1g66600_at1g66610_+_cdsdf.csv at1g66600_at1g66610_-_cdsdf.csv parallel_digest_junction_peptides.csv
```
By giving 'single_digest_junction_peptides.csv' as the output argument, ProtView adds '\_+\_' or '\_-\_' to the 
end to distinguish between strands, resulting in files named 'single_digest_junction_peptides_+\_.csv' and 
'single_digest_junction_peptides_-_.csv'

Output columns from this step are the same as the peptide processing step, with the addition of:
* **junction**: the protein coordinate of the junction that is covered

A summary table can be generated from the junction-covering peptide results. Each table includes the number of 
junction-covering peptides generated by each enzyme, the number of unique junctions that an enzyme covers (to avoid 
double counting of splice junctions shared between transcripts), and a junction coverage percentage, which is the 
percentage of the total junctions available in the isoforms being examined that are covered by an enzyme. 

Required Arguments:
* **-pept, --junction_spanning_peptides**: csv files containing the junction covering peptides
* **-cds, --cds_files**: both csv files containing extracted coding sequences for each DNA strand, 
ending in +/-_cdsdf.csv
Optional Arguments
* **-out, --output_name**: desired name of output csv summary file, default: junction_summary.csv

To summarise the junction-covering peptides identified above:
```
protview junction_summary -pept single_digest_junction_peptides_+_.csv single_digest_junction_peptides_-_.csv 
parallel_digest_junction_peptides_+_.csv parallel_digest_junction_peptides_-_.csv -cds at1g66600_at1g66610_+_cdsdf.csv 
at1g66600_at1g66610_-_cdsdf.csv
```

Output columns:
* **enzyme**: the protease or protease combination used in the digest
* **junction_spanning_peptides**: total number of peptides that cover splice junctions
* **unique_junctions_covered**: number of unique junctions covered (to avoid double counting of those shared between isoforms)
* **total_junctions_covered**: the total number of junctions covered by peptides
* **total_junction_coverage**: % of available junctions that are covered by peptides

## Unique Peptide Counter
ProtView can count the number of unique peptides per digest and append the calculation as a column to either the
general proteomic summary table or the junction-covering peptide summary table. Unique peptides are defined as those 
with sequences that can only be found in one isoform in a digest.

Required Arguments:
* **table_name**: the proteomic or junction summary statistics table for the column to be appended to
* **rpg_files**: the filtered digest results in csv format

To calculate the  number of unique peptides in the example digest and append it to the proteomic summary table:
```
protview unique_count summary_table.csv at1g66600_at1g66610_rpg_len_7_35.csv 
at1g66600_at1g66610_rpg_Trypsin_Asp-N_parallel_len_7_35.csv
```

Output column:
This column can be appended to either summary table generated
* **isoform unique peptides**

# Contact
For further information please contact <SSPuliasi@dundee.ac.uk>