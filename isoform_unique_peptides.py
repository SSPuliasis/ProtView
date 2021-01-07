import pandas as pd
from Bio import SeqIO

# note that this MUST to be ran on full proteome fasta at once
# because a seq that is unique in chr1 could come up in chr2 for example

def filter_isoform_unique_peptides(digest_results, full_proteome_fasta, output_name):
    output_df = pd.DataFrame()
    seqs_to_check = pd.read_csv(digest_results)
    for enzyme in set(seqs_to_check.enzyme):
        by_enzyme = seqs_to_check.loc[seqs_to_check['enzyme'] == enzyme]
        by_enzyme_nonredundant = by_enzyme.drop_duplicates(['sequence'], keep=False)

        seq_iso_dict = by_enzyme_nonredundant.groupby(['sequence'])['parent'].apply(list).to_dict()
        copy_seq_iso_dict = seq_iso_dict.copy()
        for sequence, isoform in seq_iso_dict.items():
            for single_isoform in isoform:
                for fasta_entry in SeqIO.parse(full_proteome_fasta, "fasta"):
                    if sequence in fasta_entry.seq:
                        if '>' + fasta_entry.id != single_isoform:
                            if sequence in copy_seq_iso_dict.keys():
                                del copy_seq_iso_dict[sequence]

        for peptide_sequence in copy_seq_iso_dict.keys():
            temp_df = by_enzyme.loc[by_enzyme['sequence'] == peptide_sequence]
            output_df = output_df.append(temp_df)
        print(enzyme, len(output_df))
    output_df.to_csv(output_name)
