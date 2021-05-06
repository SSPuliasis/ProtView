#Takes all sequences for an enzyme in the given input files, removes duplicates
# and then compares the sequences between enzymes, generating a Venn diagram that
# shows how many peptide sequences are generated by both of the given enzymes.

import pandas as pd
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
#% matplotlib inline

def compare_seqs(enzyme1, enzyme2, *inputfilenames):
    # obtain set lists of sequences for each enzyme
    enzyme1_seqs = []
    for name in inputfilenames:
        input_file = pd.read_csv(name)
        input_file = input_file.loc[(input_file.enzyme == enzyme1)]
        seqlist = input_file['sequence'].tolist()
        enzyme1_seqs.extend(seqlist)
    enzyme1_seqs = sorted(set(enzyme1_seqs))

    enzyme2_seqs = []
    for name in inputfilenames:
        input_file = pd.read_csv(name)
        input_file = input_file.loc[(input_file.enzyme == enzyme2)]
        seqlist = input_file['sequence'].tolist()
        enzyme2_seqs.extend(seqlist)
    enzyme2_seqs = sorted(set(enzyme2_seqs))

    # comparisons between set lists
    non_match_a = set(enzyme1_seqs) - set(enzyme2_seqs)
    non_match_b = set(enzyme2_seqs) - set(enzyme1_seqs)
    non_match = list(non_match_a) + list(non_match_b)

    match = list(set(enzyme1_seqs) & set(enzyme2_seqs))

    a = (len(enzyme1_seqs)) - len(match)
    c = (len(enzyme2_seqs)) - len(match)

    # sanity check. Will generate diagram if all ok
    if a + c == len(non_match):
        print(a, c, len(match))
        venn2(subsets=(a, c, len(match)), set_labels=(enzyme1, enzyme2))
    else:
        print('Mistake here. Check numbers manually')