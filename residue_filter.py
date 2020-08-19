# Filter peptides in the processed RPG output format (CSV) for peptides containing
# certain amino acids

import pandas as pd


def filter_for_residue(residue, input_file, output_file):
    unfiltered_file = pd.read_csv(input_file)
    filtered_file = pd.DataFrame()
    for index, row in unfiltered_file.iterrows():
        if residue in row['sequence']:
            filtered_file = filtered_file.append(row)
    filtered_file.to_csv(output_file)


# EXAMPLE
# input and output file names need to have csv extension and residue should be capitalised
filter_for_residue('C', 'at1g6600_at1g6610_fwd_rpg_filtered.csv', 'test1.csv')