import os
import pandas as pd

# defining the function
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
    og_rpg['peptide_start'] = og_rpg['peptide_start']+1

    parent_isoform = []
    for item in og_rpg['FASTA_description']:
        split_description = item.split()
        parent_isoform.append(split_description[0])
    og_rpg['parent'] = parent_isoform

    # save new table
    og_rpg.to_csv(input_name + '_unfiltered.csv')

    # Can now remove the rpg output that is in fasta format (not a table)
    os.remove(input_name + '.fasta')

    # OPTIONAL STEP: FILTERING BY LENGTH
    # filter rows for length between 7 and 35
    filtered_rpg = og_rpg[(og_rpg.peptide_size <= 35) & (og_rpg.peptide_size >= 7)]
    # can add in more filtering options
    filtered_rpg.to_csv(input_name + '_filtered.csv')


# to run the function
# replace input_name with the name of the rog output file to be processed in '' without .fasta at the end
process_rpg_output('at1g6600_at1g6610_rpg')