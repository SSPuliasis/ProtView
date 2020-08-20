import os
import pandas as pd

# takes the raw rpg results as input, withoug '.fasta' extension in the name
def process_rpg_output(input_name):
    og_rpg_file = input_name + '.fasta'

    # step 1: add '_' to end of every other line before converting to one line
    with open(og_rpg_file, 'r') as istr:
        with open(og_rpg_file + '.new.fasta', 'w') as ostr:
            for i, line in enumerate(istr):
                # Get rid of the trailing newline (if any).
                line = line.rstrip('\n')
                if i % 2 == 0:
                    line += '_'
                print(line, file=ostr)
    # step 2: converting to a single line
    filename = (og_rpg_file + '.new.fasta')
    open(filename + '.one.txt', 'w').write(''.join(open(filename).read().split('\n')[0:]))
    # step 3: replacing ">" with a new line and ">"
    filein = (filename + '.one.txt')
    f = open(filein, 'r')
    filedata = f.read()
    f.close()
    newdata = filedata.replace(">", "\n>")
    f = open(filein, 'w')
    f.write(newdata)
    f.close()
    # step 4: get rid of the intermediate files
    os.remove(filename)
    os.remove(og_rpg_file)
    # step 5: change the name for the output results file
    outputfile = input_name + '.fasta'
    os.rename(filein, outputfile)
    # Renaming the column headers and setting the delimiter
    headers = ["FASTA_description", "sequential_number", "enzyme", "cleavage_position",
               "peptide_size", "mol_weight", "isoelectric_point", "sequence"]
    og_rpg = pd.read_csv(input_name + '.fasta', delimiter='_', header=None, names=headers)
    og_rpg['peptide_start'] = abs(og_rpg['cleavage_position'] - og_rpg['peptide_size'])
    og_rpg['peptide_start'] = og_rpg['peptide_start'] + 1
    parent_isoform = []
    for item in og_rpg['FASTA_description']:
        split_description = item.split()
        parent_isoform.append(split_description[0])
    og_rpg['parent'] = parent_isoform
    # save new table
    og_rpg.to_csv(input_name + '_unfiltered.csv')
    # Can now remove the rpg output that is in fasta format (not a table)
    os.remove(input_name + '.fasta')


# FILTERING BY LENGTH
def filter_by_length(min_len, max_len, input_file, output_file):
    unfiltered_rpg = pd.read_csv(input_file)
    filtered_rpg = unfiltered_rpg[(unfiltered_rpg.peptide_size <= max_len) & (unfiltered_rpg.peptide_size >= min_len)]
    filtered_rpg.to_csv(output_file)


# Filter peptides in the processed RPG output format (CSV) for peptides containing
# certain amino acids
def filter_for_residue(residue, input_file, output_file):
    unfiltered_rpg = pd.read_csv(input_file)
    filtered_rpg = pd.DataFrame()
    for index, row in unfiltered_rpg.iterrows():
        if residue in row['sequence']:
            filtered_rpg = filtered_rpg.append(row)
    filtered_rpg.to_csv(output_file)


# EXAMPLE
# name of the rpg output file to be processed in '' without .fasta at the end
process_rpg_output('at1g6600_at1g6610_fwd_rpg')
filter_by_length(7, 35, 'at1g6600_a1g6610_fwd_rpg_unfiltered.csv', 'at1g6600_at1g6610_fwd_rpg_filtered_len.csv')
filter_for_residue('C', 'at1g6600_at1g6610_fwd_rpg_filtered_len.csv',
                   'at1g6600_at1g6610_fwd_rpg_filtered_len_cysteine.csv')