import pandas as pd

def junction_spanning(cds_files, rpg_file, output_file_name):
    output_name_without_extension = output_file_name.replace('.csv', '')
    for filename in cds_files:
        # POSITIVE STRAND
        if '+' in filename:
            posstrand_cdsdf = pd.read_csv(filename)
            posstrand_cdsdf = posstrand_cdsdf[posstrand_cdsdf['protein_start'] < posstrand_cdsdf[
                'protein_end']]  # remove small CDSs (<3AA) that end up with start >=end
            processed_rpg = pd.read_csv(rpg_file)  # read in the rpg output file

            dict_df = pd.DataFrame()
            dict_df['Isoform'] = posstrand_cdsdf['Parent']
            dict_df['intron_id'] = posstrand_cdsdf['intron_id']
            dict_df = dict_df[dict_df.intron_id != 'na']

            dict_df['junction'] = posstrand_cdsdf['protein_start'] - 0.5

            dict_df = dict_df.set_index(['Isoform', 'junction'])
            junctions_dict = dict_df.to_dict('index')

            junction_spanning = pd.DataFrame()
            for isoform, junction in junctions_dict.keys():
                rpg_by_isoform = processed_rpg[processed_rpg['parent'] == isoform]
                for index, row in rpg_by_isoform.iterrows():
                    if row['cleavage_position'] > junction:
                        if row['peptide_start'] < junction:
                            row['junction'] = junction
                            row['intron_id'] = junctions_dict[isoform, junction]['intron_id']
                            row = row.to_frame().transpose()
                            junction_spanning= pd.concat([junction_spanning, row])


            junction_spanning.to_csv(output_name_without_extension+'_+_.csv', index=False)

        # NEGATIVE STRAND
        elif '-' in filename:
            negstrand_cdsdf = pd.read_csv(filename)
            negstrand_cdsdf = negstrand_cdsdf[negstrand_cdsdf['protein_start']>negstrand_cdsdf[
                'protein_end']] #remove small CDSs (<3AA) that end up with start >=end (opposite
            # way round to positive strand)
            processed_rpg = pd.read_csv(rpg_file)  # read in the rpg output file

            dict_df = pd.DataFrame()
            dict_df['Isoform'] = negstrand_cdsdf['Parent']
            dict_df['intron_id'] = negstrand_cdsdf['intron_id']
            dict_df['junction'] = negstrand_cdsdf['protein_end'] - 0.5

            dict_df = dict_df[dict_df.intron_id != 'na']
            dict_df = dict_df.set_index(['Isoform', 'junction'])
            junctions_dict = dict_df.to_dict('index')

            junction_spanning = pd.DataFrame()
            for isoform, junction in junctions_dict.keys():
                rpg_by_isoform = processed_rpg[processed_rpg['parent'] == isoform]
                for index, row in rpg_by_isoform.iterrows():
                    if row['cleavage_position'] > junction:
                        if row['peptide_start'] < junction:
                            row['junction'] = junction
                            row['intron_id'] = junctions_dict[isoform, junction]['intron_id']
                            row = row.to_frame().transpose()
                            junction_spanning = pd.concat([junction_spanning, row])

            junction_spanning.to_csv(output_name_without_extension+'_-_.csv', index=False)
