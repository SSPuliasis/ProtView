# import the libraries needed
import pandas as pd

# define the functions:
def convert_pos_lhs_to_gen(prot_position, firststart, cumintron, cdsid):
    global merged_gen
    gen_position = firststart + cumintron + 3 * (prot_position - 1)
    merged_gen[cds_parent, gen_position] = {}
    merged_gen[cds_parent, gen_position] = cdsid
    # print('start position', prot_position, 'corresponds to', gen_position, 'in', cdsid)

def convert_pos_rhs_to_gen(prot_position, firststart, cumintron, cdsid):
    global merged_gen
    gen_position = firststart + cumintron + 3 * (prot_position - 1) + 2
    merged_gen[cds_parent, gen_position] = {}
    merged_gen[cds_parent, gen_position] = cdsid
    # print('end positions', prot_position, 'corresponds to', gen_position, 'in', cdsid)

def convert_neg_lhs_to_gen(prot_position, firststart, cumintron, cdsid):
    global merged_gen
    gen_position = firststart - cumintron - 3 * (prot_position - 1) - 2
    merged_gen[cds_parent, gen_position] = {}
    merged_gen[cds_parent, gen_position] = cdsid
    # print('start position', prot_position, 'corresponds to', gen_position, 'in', cdsid)

def convert_neg_rhs_to_gen(prot_position, firststart, cumintron, cdsid):
    global merged_gen
    gen_position = firststart - cumintron - 3 * (prot_position - 1)
    merged_gen[cds_parent, gen_position] = {}
    merged_gen[cds_parent, gen_position] = cdsid
    # print('end positions', prot_position, 'corresponds to', gen_position, 'in', cdsid)

############################
# POSITIVE STRAND
# importing and tidying up the datasets
input_test = pd.read_csv('at1g6600_at1g6610_+_junction_spanning.csv')
processed_cds = pd.read_csv('posstrand_cdsdf.csv')

# make an alphabetical list of the enzymes from the rpg output
enzyme_set = set(input_test['enzyme'].tolist())
enzyme_list = []
for protease in enzyme_set:
    enzyme_list.append(protease)
enzyme_list = sorted(enzyme_list)

# define this df outside the loop
overall_df = pd.DataFrame()

# filter the input dataframe by enzyme
for protease in enzyme_list:
    filtered_input_test = input_test[input_test['enzyme'] == protease]  # change thie to filtered_input_test if it works

    merged_prot = {}

    for pos, rpg_parent in zip(filtered_input_test.peptide_start, filtered_input_test.parent):
        for cds_start, cds_end, cds_id, cds_parent in zip(processed_cds.protein_start, processed_cds.protein_end,
                                                          processed_cds.cds_id, processed_cds.Parent):
            if rpg_parent.endswith(cds_parent):
                if pos in range(cds_start, cds_end + 1):
                    merged_prot[cds_parent, pos] = {}
                    merged_prot[cds_parent, pos] = cds_id
                else:
                    continue

    merged_gen = {}

    for dict_key, dict_cds_id in merged_prot.items():
        for df_cds_id, first_start, cum_intron, cds_parent in zip(processed_cds.cds_id, processed_cds.cds_first_start,
                                                                  processed_cds.cum_intron, processed_cds.Parent):
            if dict_cds_id == df_cds_id:
                if cds_parent == dict_key[0]:
                    pos_number = dict_key[1]
                    convert_pos_rhs_to_gen(pos_number, first_start, cum_intron, df_cds_id)

    isoform_start_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_start_pos', 'gen_start_pos'])
    for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
        if protkeys[0][0] == genkeys[0][0]:  # isoform
            row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                        'pept_start_pos': protkeys[0][1], 'gen_start_pos': int(genkeys[0][1])}
            isoform_start_positions_df = isoform_start_positions_df.append(row_dict, ignore_index=True)

    ###################################

    merged_prot = {}

    for pos, rpg_parent in zip(filtered_input_test.cleavage_position, filtered_input_test.parent):
        for cds_start, cds_end, cds_id, cds_parent in zip(processed_cds.protein_start, processed_cds.protein_end,
                                                          processed_cds.cds_id, processed_cds.Parent):
            if rpg_parent.endswith(cds_parent):
                if pos in range(cds_start, cds_end + 1):
                    # print('The position ', pos, ' is in the cds ', cds_id, 'with parent', rpg_parent)
                    merged_prot[cds_parent, pos] = {}
                    merged_prot[cds_parent, pos] = cds_id
                else:
                    continue

    merged_gen = {}

    for dict_key, dict_cds_id in merged_prot.items():
        for df_cds_id, first_start, cum_intron, cds_parent in zip(processed_cds.cds_id, processed_cds.cds_first_start,
                                                                  processed_cds.cum_intron, processed_cds.Parent):
            if dict_cds_id == df_cds_id:
                if cds_parent == dict_key[0]:
                    pos_number = dict_key[1]
                    convert_pos_lhs_to_gen(pos_number, first_start, cum_intron, df_cds_id)

    isoform_end_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_end_pos', 'gen_end_pos'])
    for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
        if protkeys[0][0] == genkeys[0][0]:  # isoform
            row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                        'pept_end_pos': protkeys[0][1], 'gen_end_pos': int(genkeys[0][1])}
            isoform_end_positions_df = isoform_end_positions_df.append(row_dict, ignore_index=True)

    # creating a new dataframe that combines the two
    isoform_positions_df = pd.DataFrame()
    isoform_positions_df['isoform'] = isoform_start_positions_df['isoform']  # WHY???
    isoform_positions_df['pept_start_pos'] = isoform_start_positions_df['pept_start_pos']
    isoform_positions_df['gen_start_pos'] = isoform_start_positions_df['gen_start_pos']
    isoform_positions_df['pept_end_pos'] = isoform_end_positions_df['pept_end_pos']
    isoform_positions_df['gen_end_pos'] = isoform_end_positions_df['gen_end_pos']
    isoform_positions_df['enzyme'] = protease

    overall_df = overall_df.append(isoform_positions_df)

overall_df.to_csv('at1g6600_at1g6610_+_positions_df.csv')

#############################
# NEGATIVE STRAND
# importing and tidying up the datasets
input_test = pd.read_csv('at1g6600_at1g6610_-_junction_spanning.csv')
processed_cds = pd.read_csv('negstrand_cdsdf.csv')

# make an alphabetical list of the enzymes from the rpg output
enzyme_set = set(input_test['enzyme'].tolist())
enzyme_list = []
for protease in enzyme_set:
    enzyme_list.append(protease)
enzyme_list = sorted(enzyme_list)

# define this df outside the loop
overall_df = pd.DataFrame()

# filter the input dataframe by enzyme
for protease in enzyme_list:
    filtered_input_test = input_test[input_test['enzyme'] == protease]  # change thie to filtered_input_test if it works

    merged_prot = {}

    for pos, rpg_parent in zip(filtered_input_test.cleavage_position, filtered_input_test.parent):
        for cds_start, cds_end, cds_id, cds_parent in zip(processed_cds.protein_start, processed_cds.protein_end,
                                                          processed_cds.cds_id, processed_cds.Parent):
            if rpg_parent.endswith(cds_parent):
                if pos in range(cds_end, cds_start + 1):
                    merged_prot[cds_parent, pos] = {}
                    merged_prot[cds_parent, pos] = cds_id
                else:
                    continue

    merged_gen = {}

    for dict_key, dict_cds_id in merged_prot.items():
        for df_cds_id, first_start, cum_intron, cds_parent in zip(processed_cds.cds_id, processed_cds.cds_first_end,
                                                                  processed_cds.cum_intron, processed_cds.Parent):
            if dict_cds_id == df_cds_id:
                if cds_parent == dict_key[0]:
                    pos_number = dict_key[1]
                    convert_neg_lhs_to_gen(pos_number, first_start, cum_intron, df_cds_id)

    isoform_start_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_start_pos', 'gen_start_pos'])
    for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
        if protkeys[0][0] == genkeys[0][0]:  # isoform
            row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                        'pept_start_pos': protkeys[0][1], 'gen_start_pos': int(genkeys[0][1])}
            isoform_start_positions_df = isoform_start_positions_df.append(row_dict, ignore_index=True)

    ##############################################
    merged_prot = {}

    for pos, rpg_parent in zip(filtered_input_test.peptide_start, filtered_input_test.parent):
        for cds_start, cds_end, cds_id, cds_parent in zip(processed_cds.protein_start, processed_cds.protein_end,
                                                          processed_cds.cds_id, processed_cds.Parent):
            if rpg_parent.endswith(cds_parent):
                if pos in range(cds_end, cds_start + 1):
                    # print('The position ', pos, ' is in the cds ', cds_id, 'with parent', rpg_parent)
                    merged_prot[cds_parent, pos] = {}
                    merged_prot[cds_parent, pos] = cds_id
                else:
                    continue
    merged_gen = {}

    for dict_key, dict_cds_id in merged_prot.items():
        for df_cds_id, first_start, cum_intron, cds_parent in zip(processed_cds.cds_id, processed_cds.cds_first_end,
                                                                  processed_cds.cum_intron, processed_cds.Parent):
            if dict_cds_id == df_cds_id:
                if cds_parent == dict_key[0]:
                    pos_number = dict_key[1]
                    convert_neg_rhs_to_gen(pos_number, first_start, cum_intron, df_cds_id)

    isoform_end_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_end_pos', 'gen_end_pos'])
    for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
        if protkeys[0][0] == genkeys[0][0]:  # isoform
            row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                        'pept_end_pos': protkeys[0][1], 'gen_end_pos': int(genkeys[0][1])}
            isoform_end_positions_df = isoform_end_positions_df.append(row_dict, ignore_index=True)

    # creating a new dataframe that combines the two
    isoform_positions_df = pd.DataFrame()
    isoform_positions_df['isoform'] = isoform_start_positions_df['isoform']  # WHY???
    isoform_positions_df['pept_start_pos'] = isoform_start_positions_df['pept_start_pos']
    isoform_positions_df['gen_start_pos'] = isoform_start_positions_df['gen_start_pos']
    isoform_positions_df['pept_end_pos'] = isoform_end_positions_df['pept_end_pos']
    isoform_positions_df['gen_end_pos'] = isoform_end_positions_df['gen_end_pos']
    isoform_positions_df['enzyme'] = protease

    overall_df = overall_df.append(isoform_positions_df)

overall_df.to_csv('at1g6600_at1g6610_-_positions_df.csv')