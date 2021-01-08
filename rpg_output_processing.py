import os
import pandas as pd
import shutil
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# takes the raw rpg results as input, withoug '.fasta' extension in the name
# and replace isoelectric point by rpg with the biopython iso point
# remove '>' from start of parent and correct in all subsequent scripts
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
    og_rpg = pd.read_csv(input_name + '.fasta', delimiter='_', header=None, names=headers, keep_default_na=False)
    og_rpg = og_rpg.drop('sequential_number', axis =1)
    og_rpg['peptide_start'] = abs(og_rpg['cleavage_position'] - og_rpg['peptide_size'])
    og_rpg['peptide_start'] = og_rpg['peptide_start'] + 1
    parent_isoform = []
    for item in og_rpg['FASTA_description']:
        split_description = item.split()
        parent_isoform.append(split_description[0])
    og_rpg['parent'] = parent_isoform
    # removing rows with 'na' and sequences containing 'x' and '*'
    og_rpg = og_rpg.dropna()
    og_rpg = og_rpg[~og_rpg.sequence.str.contains('X')]
    og_rpg = og_rpg[~og_rpg.sequence.str.contains('\*')]
    # save new table
    og_rpg.to_csv(input_name + '_unfiltered.csv')
    # Can now remove the rpg output that is in fasta format (not a table)
    os.remove(input_name + '.fasta')

#miscleavage needs to be done BEFORE any filtering/analysis
def create_miscleavage(unfiltered_rpg_filename, mc, output_name):
    input_file = pd.read_csv(unfiltered_rpg_filename, keep_default_na=False)
    master_df = pd.DataFrame()
    for enzyme in sorted(set(input_file.enzyme)):
        by_enzyme = input_file.loc[(input_file.enzyme == enzyme)]
        for isoform in sorted(set(by_enzyme.parent)):
            baby_df = pd.DataFrame()
            by_isoform_and_enzyme = by_enzyme.loc[(by_enzyme.parent == isoform)]
            sequencelist = []
            FASTA_description = []
            start_positions = []
            cleavage_positions = []
            for index, pept in zip(by_isoform_and_enzyme.index, by_isoform_and_enzyme.sequence):
                sequencelist.append(pept)
                n = 1
                FASTA_description.append(by_isoform_and_enzyme.FASTA_description[index])
                start_positions.append(by_isoform_and_enzyme.peptide_start[index])
                cleavage_positions.append(by_isoform_and_enzyme.cleavage_position[index])
                while n <= mc:
                    if index + n in by_isoform_and_enzyme.index:
                        pept += input_file.sequence[index + n]
                        sequencelist.append(pept)
                        FASTA_description.append(by_isoform_and_enzyme.FASTA_description[index])
                        start_positions.append(by_isoform_and_enzyme.peptide_start[index])
                        cleavage_positions.append(input_file.cleavage_position[index + n])
                        n += 1
                    else:
                        break
            # start appending columns to a new df
            baby_df['sequence'] = sequencelist
            baby_df['FASTA_description'] = FASTA_description
            baby_df['enzyme'] = enzyme
            baby_df['parent'] = isoform
            baby_df['peptide_size'] = baby_df['sequence'].apply(len)
            baby_df['peptide_start'] = start_positions
            baby_df['cleavage_position'] = cleavage_positions
            master_df = master_df.append(baby_df)

    mol_weights = []
    iso_points = []
    for sequence in master_df['sequence']:
        analysed_seq = ProteinAnalysis(sequence)
        mol_w = analysed_seq.molecular_weight()
        iso_p = analysed_seq.isoelectric_point()
        mol_weights.append(mol_w)
        iso_points.append(iso_p)
    master_df['isoelectric_point'] = iso_points
    master_df['molecular_weight'] = mol_weights
    master_df.to_csv(output_name)

# combining digests in parallel
def create_parallel_digest(input_file, output_file, *enzymes):
    filein = pd.read_csv(input_file)
    filein = filein.drop('Unnamed: 0', axis = 1) # does this need to be removed earlier in rpg output processing?
    parallel_df = pd.DataFrame()
    newname= ""
    for protease_name in list(enzymes):
        parallel_df = parallel_df.append(filein.loc[(filein.enzyme == protease_name)])
        newname +=(':'+protease_name)
        parallel_enzymes = newname[1:]
    parallel_df['enzyme'] = parallel_enzymes
    parallel_df = pd.DataFrame.drop_duplicates(parallel_df)
    parallel_df.to_csv(output_file)
    print('The parallel enzyme combination for '+parallel_enzymes+' has been saved as '+output_file)

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
    filtered_rpg = unfiltered_rpg[unfiltered_rpg['sequence'].str.contains(residue)]
    filtered_rpg.to_csv(output_file)

# merge rpg outputs into one file
def merge_files(output_file_name, *input_file_names):
    with open(output_file_name, "wb") as ouput:
        with open(input_file_names[0], "rb") as file:
            shutil.copyfileobj(file, ouput)
        if len(input_file_names[1:]) > 1:
           for name in input_file_names[1:]:
               with open(name, "rb") as file:
                   file.readline()
                   shutil.copyfileobj(file, ouput)
        elif len(input_file_names[1:]) == 1:
            name = ''.join(input_file_names[1:])
            print(name)
            print(type(name))
            with open(name, "rb") as file:
                file.readline()
                shutil.copyfileobj(file, ouput)
    print("The files have been merged")

# EXAMPLE
# name of the rpg output file to be processed in '' without .fasta at the end
process_rpg_output('at1g6600_at1g6610_fwd_rpg')
filter_by_length(7, 35, 'at1g6600_a1g6610_fwd_rpg_unfiltered.csv', 'at1g6600_at1g6610_fwd_rpg_filtered_len.csv')
filter_for_residue('C', 'at1g6600_at1g6610_fwd_rpg_filtered_len.csv',
                   'at1g6600_at1g6610_fwd_rpg_filtered_len_cysteine.csv')
create_parallel_digest('at1g6600_at1g6610_rev_rpg_filtered.csv', 'tryp:lysc_rev_rpg_filtered.csv', 'Trypsin', 'Lys-C')