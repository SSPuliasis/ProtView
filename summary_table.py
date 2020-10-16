import pandas as pd
from Bio import SeqIO

# enter the names of all of the input files as lists
unfiltered_rpg_files = ["chr1_fwd_rpg_unfiltered.csv", "chr2_fwd_rpg_unfiltered.csv",
                        "chr3_fwd_rpg_unfiltered.csv", "chr4_fwd_rpg_unfiltered.csv",
                        "chr5_fwd_rpg_unfiltered.csv",
                        "chr1_rev_rpg_unfiltered.csv", "chr2_rev_rpg_unfiltered.csv",
                        "chr3_rev_rpg_unfiltered.csv", "chr4_rev_rpg_unfiltered.csv",
                        "chr5_rev_rpg_unfiltered.csv"]
filtered_rpg_files = ['chr1_fwd_rpg_filtered_len_cysteine.csv', 'chr2_fwd_rpg_filtered_len_cysteine.csv',
                      'chr3_fwd_rpg_filtered_len_cysteine.csv', 'chr4_fwd_rpg_filtered_len_cysteine.csv',
                      'chr5_fwd_rpg_filtered_len_cysteine.csv',
                      'chr1_rev_rpg_filtered_len_cysteine.csv', 'chr2_rev_rpg_filtered_len_cysteine.csv',
                      'chr3_rev_rpg_filtered_len_cysteine.csv', 'chr4_rev_rpg_filtered_len_cysteine.csv',
                      'chr5_rev_rpg_filtered_len_cysteine.csv']
fasta_files = ['chr1_fwd.fasta', 'chr2_fwd.fasta', 'chr3_fwd.fasta', 'chr4_fwd.fasta', 'chr5_fwd.fasta',
               'chr1_rev.fasta', 'chr2_rev.fasta', 'chr3_rev.fasta', 'chr4_rev.fasta', 'chr5_rev.fasta']


# creating the summary table
def create_summary_table(unfiltered_rpg_files, filtered_rpg_files, fasta_files, output_name):
    summary_table = pd.DataFrame()

    # TOTAL PEPTIDES GENERATED
    testname = []
    testenzyme = []
    testtotal = []
    for name in unfiltered_rpg_files:
        inputfile = pd.read_csv(name)
        x = inputfile['enzyme'].tolist()
        print(name)
        for enzyme in set(x):
            testname.append(name)
            testenzyme.append(enzyme)
            testtotal.append(x.count(enzyme))

    temp_table = pd.DataFrame()
    temp_table["input"] = testname
    temp_table['enzyme'] = testenzyme
    temp_table['total peptides'] = testtotal
    # they need to be ordered by enzyme before the lists are compared
    temp_table = temp_table.sort_values(by=['enzyme'])

    # back to lists, but ordered by enzyme
    testenzyme = temp_table['enzyme'].tolist()
    testtotal = temp_table['total peptides'].tolist()

    for index, (a, c) in enumerate(zip(testenzyme, testtotal)):
        # print(index, a, c)
        # print(testenzyme[index], testenzyme[index - 1], testtotal[index])
        if testenzyme[index] == testenzyme[index - 1]:
            testtotal[index] += testtotal[index - 1]
        elif testenzyme[index] != testenzyme[index - 1]:
            continue

    temp_table['enzyme'] = testenzyme
    temp_table['cumulative_totals'] = testtotal

    enzymelist = []
    totallist = []
    temp_table['next_enzyme'] = temp_table['enzyme'].shift(periods=-1)
    for index, (a, b, c) in enumerate(zip(temp_table.enzyme, temp_table.next_enzyme, temp_table.cumulative_totals)):
        if a != b:
            enzymelist.append(a)
            totallist.append(c)
        else:
            continue

    summary_table['enzyme'] = enzymelist
    summary_table['total peptides'] = totallist

    #  TOTAL FILTERED PEPTIDES GENERATED
    testname = []
    testenzyme = []
    testtotal = []
    for name in filtered_rpg_files:
        inputfile = pd.read_csv(name)
        x = inputfile['enzyme'].tolist()
        print(name)
        for enzyme in set(x):
            testname.append(name)
            testenzyme.append(enzyme)
            testtotal.append(x.count(enzyme))

    temp_table = pd.DataFrame()
    temp_table["input"] = testname
    temp_table['enzyme'] = testenzyme
    temp_table['filtered peptides'] = testtotal
    # they need to be ordered by enzyme before the lists are compared
    temp_table = temp_table.sort_values(by=['enzyme'])
    # back to lists, but ordered by enzyme
    testenzyme = temp_table['enzyme'].tolist()
    testtotal = temp_table['filtered peptides'].tolist()

    for index, (a, c) in enumerate(zip(testenzyme, testtotal)):
        if testenzyme[index] == testenzyme[index - 1]:
            testtotal[index] += testtotal[index - 1]
        elif testenzyme[index] != testenzyme[index - 1]:
            continue

    temp_table['enzyme'] = testenzyme
    temp_table['cumulative_totals'] = testtotal

    totallist = []
    temp_table['next_enzyme'] = temp_table['enzyme'].shift(periods=-1)
    for index, (a, b, c) in enumerate(zip(temp_table.enzyme, temp_table.next_enzyme, temp_table.cumulative_totals)):
        if a != b:
            totallist.append(c)
        else:
            continue
    summary_table['filtered peptides'] = totallist
    ###

    # PROTEIN SEQUENCE COVERAGE
    coverage_table = pd.DataFrame()

    protein_len_sum = 0
    for name in fasta_files:
        sizes = [len(rec) for rec in SeqIO.parse(name, "fasta")]
        print(name, sum(sizes))
        protein_len_sum += sum(sizes)

    for name in filtered_rpg_files:
        rpg_file = pd.read_csv(name)
        rpg_file = rpg_file.sort_values(['parent', 'enzyme', 'peptide_start'])

        for protein in set(rpg_file.parent):
            protein_df = rpg_file.loc[(rpg_file.parent == protein)]
            protein_df = protein_df.rename(columns={'cleavage_position': 'end', 'peptide_start': 'start'})

            for enzyme in set(protein_df.enzyme):
                enzyme_df = protein_df.loc[(protein_df.enzyme == enzyme)]
                covered_regions_list = []
                minval = -10000
                maxval = -100000
                enzyme_df['covered_regions'] = enzyme_df[['start', 'end']].values.tolist()
                arraylist = enzyme_df['covered_regions'].tolist()
                for i in range(len(arraylist)):
                    a = arraylist[i]
                    if a[0] > maxval:
                        if i != 0:
                            covered_regions_list.append([minval, maxval])
                        maxval = a[1]
                        minval = a[0]
                    else:
                        if a[1] >= maxval:
                            maxval = a[1]
                if maxval != -100000 and [minval, maxval] not in covered_regions_list:
                    covered_regions_list.append([minval, maxval])
                temp_table = pd.DataFrame()
                temp_table['covered_regions'] = covered_regions_list
                temp_table['protein'] = protein
                temp_table['enzyme'] = enzyme
                temp_table['input_file'] = name
                temp_table[['start', 'end']] = pd.DataFrame(temp_table.covered_regions.tolist(), index=temp_table.index)
                temp_table['peptide_length'] = temp_table['end'] - temp_table['start'] + 1
                temp_table = temp_table.drop(['start', 'end'], axis=1)
                coverage_table = coverage_table.append(temp_table, sort=True)
    coverage_summary_table = pd.DataFrame(coverage_table.groupby(['enzyme'], as_index=False)['peptide_length'].sum())
    total_protein_length = protein_len_sum
    coverage_summary_table['protein_length'] = total_protein_length
    coverage_summary_table['coverage'] = (coverage_summary_table['peptide_length'] / coverage_summary_table[
        'protein_length']) * 100
    # write.csv(coverage_summary_table, 'coverage_summary_table.csv')

    summary_table['coverage'] = coverage_summary_table['coverage']

    # MEAN PEPTIDE LENGTHS
    digest_count = []
    for enzyme in coverage_summary_table['enzyme']:
        digests_combined = enzyme.count(':') + 1
        digest_count.append(digests_combined)

    coverage_summary_table['combined_digests'] = digest_count

    summary_table['mean length'] = (coverage_summary_table['protein_length'] * coverage_summary_table[
        'combined_digests']) / summary_table['total peptides']
    # summary_table

    # changing the column orders so that mena length is after total no. of unfiltered pepts
    summary_table = summary_table[['enzyme', 'total peptides', 'mean length', 'filtered peptides', 'coverage']]

    summary_table.to_csv(output_name)

##### WITH RESIDUE FILTERING
# FILTERED RPG FILES ABOVE NEED TO HAVE BEEN FILTERED FOR RESIDUE FOR THIS TO WORK
# COUNTING FREQ OF RESIDUES IN THE SEQUENCES
residue_freq_in_seq = 0
for file in fasta_files:
    for rec in SeqIO.parse(file, "fasta"):
        freq_in_fasta = rec.seq.count("C")
        residue_freq_in_seq += freq_in_fasta

# COUNTING FREQ OF RESIDUES IN FILTERED DATA
# NEED SEPARATE FOR EACH ENZYME, MAKE NEW TABLE

newlist = []
enzymelist = []
freqlist = []
temp_table = pd.DataFrame()
for name in filtered_rpg_files:
    rpg_file = pd.read_csv(name)
    x = rpg_file['enzyme'].tolist()
    for enzyme in sorted(set(x)):
        seqlist = []
        rpg_by_enzyme = rpg_file[rpg_file.enzyme == enzyme]
        seqlist = rpg_by_enzyme['sequence'].tolist()
        frequency = 0
        for sequence in seqlist:
            frequency += sequence.count("C")
        freqlist.append(frequency)
        enzymelist.append(enzyme)
temp_table['enzyme'] = enzymelist
temp_table['residue_freq'] = freqlist

x = temp_table['enzyme'].tolist()
enzymelist = []
freqlist = []
for enzyme in sorted(set(x)):
    temp_table_by_enzyme = temp_table[temp_table.enzyme == enzyme]
    enzymelist.append(enzyme)
    freqlist.append(sum(temp_table_by_enzyme['residue_freq']))  # append to datafrane

# WHERE IS TEMP TABLE BY ENZYME?!?!

temp_table = pd.DataFrame()
temp_table['enzyme'] = enzymelist
temp_table['residue_freq'] = freqlist

# ADDING COVERAGE COLUMN TO EXISTING TABLE
summary_table['residue_coverage_%'] = (temp_table['residue_freq'] / residue_freq_in_seq) * 100

summary_table.to_csv('cysteine_summary.csv')