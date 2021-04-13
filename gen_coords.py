# import the libraries needed
import pandas as pd

def calculate_gen_coords(peptides_file, cds_files):
    input_test = pd.read_csv(peptides_file)
    file_name = peptides_file.replace('.csv', '')
    for file in cds_files:
        if file.endswith('_+_cdsdf.csv'):
            posstrand_cdsdf = file
        elif file.endswith('_-_cdsdf.csv'):
            negstrand_cdsdf = file
        else:
            print('error in cds file names, they should end with _+/-_cdsdf')

    ############################
    # POSITIVE STRAND
    # importing and tidying up the datasets
    processed_cds = pd.read_csv(posstrand_cdsdf)

    isoformlist = set(processed_cds['Parent'].tolist())
    strand_peptides = input_test[input_test['parent'].isin(isoformlist)]

    # make an alphabetical list of the enzymes from the rpg output
    enzyme_set = set(strand_peptides['enzyme'].tolist())
    enzyme_list = []
    for protease in enzyme_set:
        enzyme_list.append(protease)
    enzyme_list = sorted(enzyme_list)

    # define this df outside the loop
    overall_df = pd.DataFrame()

    # filter the input dataframe by enzyme
    for protease in enzyme_list:
        peptides_by_enzyme = strand_peptides[strand_peptides['enzyme'] == protease]  # change thie to peptides_by_enzyme if it works

        merged_prot = {}

        for pos, rpg_parent in zip(peptides_by_enzyme.peptide_start, peptides_by_enzyme.parent):
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
                        gen_position = first_start + cum_intron + 3 * (pos_number - 1) + 2
                        merged_gen[cds_parent, gen_position] = {}
                        merged_gen[cds_parent, gen_position] = df_cds_id
                        # print('end positions', prot_position, 'corresponds to', gen_position, 'in', cdsid)

        isoform_start_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_start_pos', 'gen_start_pos'])
        for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
            if protkeys[0][0] == genkeys[0][0]:  # isoform
                row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                            'pept_start_pos': protkeys[0][1], 'gen_start_pos': int(genkeys[0][1])}
                isoform_start_positions_df = isoform_start_positions_df.append(row_dict, ignore_index=True)

        ###################################

        merged_prot = {}

        for pos, rpg_parent in zip(peptides_by_enzyme.cleavage_position, peptides_by_enzyme.parent):
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
                        gen_position = first_start + cum_intron + 3 * (pos_number - 1)
                        merged_gen[cds_parent, gen_position] = {}
                        merged_gen[cds_parent, gen_position] = df_cds_id
                        # print('start position', prot_position, 'corresponds to', gen_position, 'in', cdsid)

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

    overall_df.to_csv(file_name+'_+_positions_df.csv', index=False)

    #############################
    # NEGATIVE STRAND
    # importing and tidying up the datasets
    processed_cds = pd.read_csv(negstrand_cdsdf)

    isoformlist = set(processed_cds['Parent'].tolist())
    strand_peptides = input_test[input_test['parent'].isin(isoformlist)]

    # make an alphabetical list of the enzymes from the rpg output
    enzyme_set = set(strand_peptides['enzyme'].tolist())
    enzyme_list = []
    for protease in enzyme_set:
        enzyme_list.append(protease)
    enzyme_list = sorted(enzyme_list)

    # define this df outside the loop
    overall_df = pd.DataFrame()

    # filter the input dataframe by enzyme
    for protease in enzyme_list:
        peptides_by_enzyme = input_test[input_test['enzyme'] == protease]  # change thie to peptides_by_enzyme if it works

        merged_prot = {}

        for pos, rpg_parent in zip(peptides_by_enzyme.cleavage_position, peptides_by_enzyme.parent):
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
                        gen_position = first_start - cum_intron - 3 * (pos_number - 1) - 2
                        merged_gen[cds_parent, gen_position] = {}
                        merged_gen[cds_parent, gen_position] = df_cds_id
                        # print('start position', prot_position, 'corresponds to', gen_position, 'in', cdsid)

        isoform_start_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_start_pos', 'gen_start_pos'])
        for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
            if protkeys[0][0] == genkeys[0][0]:  # isoform
                row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                            'pept_start_pos': protkeys[0][1], 'gen_start_pos': int(genkeys[0][1])}
                isoform_start_positions_df = isoform_start_positions_df.append(row_dict, ignore_index=True)

        ##############################################
        merged_prot = {}

        for pos, rpg_parent in zip(peptides_by_enzyme.peptide_start, peptides_by_enzyme.parent):
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
                        gen_position = first_start - cum_intron - 3 * (pos_number - 1)
                        merged_gen[cds_parent, gen_position] = {}
                        merged_gen[cds_parent, gen_position] = df_cds_id
                        # print('end positions', prot_position, 'corresponds to', gen_position, 'in', cdsid)

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

    overall_df.to_csv(file_name+'_-_positions_df.csv', index=False)