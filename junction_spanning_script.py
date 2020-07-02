import pandas as pd

# SET NAME TO THE RPG FILE, DEFAULT IS THE FILTERED VERSION
input_name = 'chr5'

posstrand_cdsdf = pd.read_csv(input_name + '_+_cdsdf.csv')
negstrand_cdsdf = pd.read_csv(input_name + '_-_cdsdf.csv')
# default rpg file here is the filtered one. Can change 'filtered_rpg' to 'og_rpg

# POSITIVE STRAND
filtered_rpg = pd.read_csv(input_name + '_fwd_rpg_filtered.csv')  # read in the rpg output file

dict_df = pd.DataFrame()
dict_df['Isoform'] = posstrand_cdsdf['Parent']
dict_df['intron_id'] = posstrand_cdsdf['intron_id']
dict_df['junction'] = posstrand_cdsdf['protein_start'] - 0.5

dict_df = dict_df[dict_df.intron_id != 'na']
dict_df = dict_df.set_index(['Isoform', 'junction'])

junctions_dict = dict_df.to_dict('index')

parent_isoform = []
for item in filtered_rpg['FASTA_description']:
    split_description = item.split()
    parent_isoform.append(split_description[0])
filtered_rpg['parent'] = parent_isoform

junction_spanning = pd.DataFrame()
for isoform, junction in junctions_dict.keys():
    rpg_by_isoform = filtered_rpg[filtered_rpg['parent'] == '>' + isoform]
    for index, row in rpg_by_isoform.iterrows():
        if row['cleavage_position'] > junction:
            if row['peptide_start'] < junction:
                row['junction'] = junction
                row['intron_id'] = junctions_dict[isoform, junction]['intron_id']
                junction_spanning = junction_spanning.append(row)

junction_spanning = junction_spanning.drop(['Unnamed: 0'], axis=1)
junction_spanning = junction_spanning.drop(['sequential_number'], axis=1)

junction_spanning.to_csv(input_name + '_+_junction_spanning.csv')

# NEGATIVE STRAND
filtered_rpg = pd.read_csv(input_name + '_rev_rpg_filtered.csv')  # read in the rpg output file

dict_df = pd.DataFrame()
dict_df['Isoform'] = negstrand_cdsdf['Parent']
dict_df['intron_id'] = negstrand_cdsdf['intron_id']
dict_df['junction'] = negstrand_cdsdf['protein_end'] - 0.5

dict_df = dict_df[dict_df.intron_id != 'na']
dict_df = dict_df.set_index(['Isoform', 'junction'])

junctions_dict = dict_df.to_dict('index')

parent_isoform = []
for item in filtered_rpg['FASTA_description']:
    split_description = item.split()
    parent_isoform.append(split_description[0])
filtered_rpg['parent'] = parent_isoform

junction_spanning = pd.DataFrame()
for isoform, junction in junctions_dict.keys():
    rpg_by_isoform = filtered_rpg[filtered_rpg['parent'] == '>' + isoform]
    for index, row in rpg_by_isoform.iterrows():
        if row['cleavage_position'] > junction:
            if row['peptide_start'] < junction:
                row['junction'] = junction
                row['intron_id'] = junctions_dict[isoform, junction]['intron_id']
                junction_spanning = junction_spanning.append(row)

junction_spanning = junction_spanning.drop(['Unnamed: 0'], axis=1)
junction_spanning = junction_spanning.drop(['sequential_number'], axis=1)

junction_spanning.to_csv(input_name + '_-_junction_spanning.csv')