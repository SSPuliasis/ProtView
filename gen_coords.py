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

    print(posstrand_cdsdf, negstrand_cdsdf)
    ############################
    # POSITIVE STRAND
    # importing and tidying up the datasets
    processed_cds = pd.read_csv(posstrand_cdsdf)

    isoformlist = set(processed_cds['Parent'].tolist())
    strand_peptides = input_test[input_test['parent'].isin(isoformlist)]

    # make an alphabetical list of the enzymes from the rpg output
    enzyme_list = sorted(set(strand_peptides['enzyme'].tolist()))

    # define this df outside the loop
    overall_df = pd.DataFrame()

    # filter the input dataframe by enzyme
    for protease in enzyme_list:
        peptides_by_enzyme = strand_peptides[strand_peptides['enzyme'] == protease]

        merged_prot = {}
        merged_gen = {}

        for index, cds_row in processed_cds.iterrows():
            peptides_by_enzyme_parent = peptides_by_enzyme.loc[(peptides_by_enzyme['parent'].str.endswith(cds_row.Parent))
                                                               & (peptides_by_enzyme.peptide_start <= cds_row.protein_end)
                                                               & (peptides_by_enzyme.peptide_start >= cds_row.protein_start)]
            for pos in peptides_by_enzyme_parent.peptide_start:
                merged_prot[cds_row.Parent, pos] = cds_row.cds_id

            for dict_key, dict_cds_id in merged_prot.items():
                if dict_cds_id == cds_row.cds_id and cds_row.Parent == dict_key[0]:
                    pos_number = dict_key[1]
                    gen_position = cds_row.cds_first_start + cds_row.cum_intron + 3 * (pos_number - 1)
                    merged_gen[cds_row.Parent, gen_position] = cds_row.cds_id
                    # print('end positions', prot_position, 'corresponds to', gen_position, 'in', cdsid)

        isoform_start_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_start_pos', 'gen_start_pos'])
        for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
            if protkeys[0][0] == genkeys[0][0]:  # isoform
                row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                            'pept_start_pos': protkeys[0][1], 'gen_start_pos': int(genkeys[0][1])}
                isoform_start_positions_df = isoform_start_positions_df.append(row_dict, ignore_index=True)

        ###################################

        merged_prot = {}
        merged_gen = {}

        for index, cds_row in processed_cds.iterrows():
            peptides_by_enzyme_parent = peptides_by_enzyme.loc[(peptides_by_enzyme['parent'].str.endswith(cds_row.Parent))
                                                               & (peptides_by_enzyme.cleavage_position <= cds_row.protein_end)
                                                               & (peptides_by_enzyme.cleavage_position >= cds_row.protein_start)]
            for pos in peptides_by_enzyme_parent.cleavage_position:
                merged_prot[cds_row.Parent, pos] = cds_row.cds_id

            for dict_key, dict_cds_id in merged_prot.items():
                if dict_cds_id == cds_row.cds_id and cds_row.Parent == dict_key[0]:
                    pos_number = dict_key[1]
                    gen_position = cds_row.cds_first_start + cds_row.cum_intron + 3 * (pos_number - 1) + 2
                    merged_gen[cds_row.Parent, gen_position] = cds_row.cds_id
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

    overall_df.to_csv(file_name+'_+_positions_df.csv', index=False)

    #############################
    # NEGATIVE STRAND
    # importing and tidying up the datasets
    processed_cds = pd.read_csv(negstrand_cdsdf)

    isoformlist = set(processed_cds['Parent'].tolist())
    strand_peptides = input_test[input_test['parent'].isin(isoformlist)]

    # make an alphabetical list of the enzymes from the rpg output
    enzyme_list = sorted(set(strand_peptides['enzyme'].tolist()))

    # define this df outside the loop
    overall_df = pd.DataFrame()

    # filter the input dataframe by enzyme
    for protease in enzyme_list:
        peptides_by_enzyme = input_test[input_test['enzyme'] == protease]

        merged_prot = {}
        merged_gen = {}

        for index, cds_row in processed_cds.iterrows():
            peptides_by_enzyme_parent = peptides_by_enzyme.loc[(peptides_by_enzyme['parent'].str.endswith(cds_row.Parent))
                                                               & (peptides_by_enzyme.cleavage_position >= cds_row.protein_end)
                                                               & (peptides_by_enzyme.cleavage_position <= cds_row.protein_start)]
            for pos in peptides_by_enzyme_parent.cleavage_position:
                merged_prot[cds_row.Parent, pos] = cds_row.cds_id

            for dict_key, dict_cds_id in merged_prot.items():
                if dict_cds_id == cds_row.cds_id and cds_row.Parent == dict_key[0]:
                    pos_number = dict_key[1]
                    gen_position = cds_row.cds_first_end - cds_row.cum_intron - 3 * (pos_number - 1) - 2
                    merged_gen[cds_row.Parent, gen_position] = cds_row.cds_id
                    # print('end positions', prot_position, 'corresponds to', gen_position, 'in', cdsid)

        isoform_start_positions_df = pd.DataFrame(columns=['isoform', 'cds_id', 'pept_start_pos', 'gen_start_pos'])
        for protkeys, genkeys in zip(merged_prot.items(), merged_gen.items()):
            if protkeys[0][0] == genkeys[0][0]:  # isoform
                row_dict = {'isoform': protkeys[0][0], 'cds_id': protkeys[1],
                            'pept_start_pos': protkeys[0][1], 'gen_start_pos': int(genkeys[0][1])}
                isoform_start_positions_df = isoform_start_positions_df.append(row_dict, ignore_index=True)

        ##############################################
        merged_prot = {}
        merged_gen = {}

        for index, cds_row in processed_cds.iterrows():
            peptides_by_enzyme_parent = peptides_by_enzyme.loc[(peptides_by_enzyme['parent'].str.endswith(cds_row.Parent))
                                                               & (peptides_by_enzyme.peptide_start >= cds_row.protein_end)
                                                               & (peptides_by_enzyme.peptide_start <= cds_row.protein_start)]
            for pos in peptides_by_enzyme_parent.peptide_start:
                merged_prot[cds_row.Parent, pos] = cds_row.cds_id

            for dict_key, dict_cds_id in merged_prot.items():
                if dict_cds_id == cds_row.cds_id and cds_row.Parent == dict_key[0]:
                    pos_number = dict_key[1]
                    gen_position = cds_row.cds_first_end - cds_row.cum_intron - 3 * (pos_number - 1)
                    merged_gen[cds_row.Parent, gen_position] = cds_row.cds_id
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