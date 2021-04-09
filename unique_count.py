import pandas as pd


def add_unique_pept_column(table_name, rpg_files):
    unique_count_df = pd.DataFrame(columns=['enzyme', 'isoform_unique_peptides'])
    index = 0
    summary_table = pd.read_csv(table_name)
    for enzyme in sorted(set(summary_table.enzyme)):
        merged_by_enzyme = pd.DataFrame()
        for name in rpg_files:
            print(name)
            seqs_to_check = pd.read_csv(name)
            by_enzyme = seqs_to_check.loc[seqs_to_check['enzyme'] == enzyme]
            merged_by_enzyme = merged_by_enzyme.append(by_enzyme)

        seq_iso_dict = merged_by_enzyme.groupby(['sequence'])['parent'].apply(list).to_dict()
        count = sum(len(v) == 1 for v in seq_iso_dict.values())
        unique_count_df.loc[index] = [enzyme, count]
        index += 1

    unique_count_df = unique_count_df.sort_values('enzyme')
    summary_table = summary_table.sort_values('enzyme').reset_index(drop=True)
    summary_table['isoform unique peptides'] = unique_count_df['isoform_unique_peptides']
    summary_table.to_csv(table_name)