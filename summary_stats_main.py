import pandas as pd
from Bio import SeqIO
import statistics
import re

# enter the names of all of the input files as lists
#unfiltered_rpg_files = ["trypsin_gluc_rpg.csv"]
#filtered_rpg_files = ["trypsin_gluc_rpg_len_7_35.csv"]
#fasta_files = ['at1g6600_at1g6610.fasta']


# creating the summary table
def create_summary_table(unfiltered_rpg_files, filtered_rpg_files, fasta_files, output_name):
    summary_table = pd.DataFrame()

    # TOTAL PEPTIDES GENERATED
    testname = []
    testenzyme = []
    testtotal = []
    xlist = []
    all_peptide_lengths = []
    for name in unfiltered_rpg_files:
        inputfile = pd.read_csv(name)
        ylist = inputfile['enzyme'].tolist()
        xlist.extend(ylist)
        peptide_lengths = inputfile['peptide_size'].tolist()
        all_peptide_lengths.extend(peptide_lengths)
        for enzyme in sorted(set(xlist)):
            testname.append(name)
            testenzyme.append(enzyme)
            testtotal.append(xlist.count(enzyme))

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
            if testtotal[index] != testtotal[index - 1]:
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

    # MEDIAN PEPTIDE LENGTHS
    length_dist = pd.DataFrame()
    length_dist['enzyme'] = sorted(xlist)
    length_dist['peptide_size'] = all_peptide_lengths
    medians = []
    for enzyme in sorted(set(xlist)):
        median_length = statistics.median(length_dist.loc[(length_dist.enzyme == enzyme)].peptide_size)
        medians.append(median_length)
    summary_table['median length'] = medians

    #  TOTAL FILTERED PEPTIDES GENERATED
    testname = []
    testenzyme = []
    testtotal = []
    xlist = []
    all_peptide_lengths = []
    for name in filtered_rpg_files:
        inputfile = pd.read_csv(name)
        ylist = inputfile['enzyme'].tolist()
        xlist.extend(ylist)
        peptide_lengths = inputfile['peptide_size'].tolist()
        all_peptide_lengths.extend(peptide_lengths)
        for enzyme in sorted(set(xlist)):
            testname.append(name)
            testenzyme.append(enzyme)
            testtotal.append(xlist.count(enzyme))

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
            if testtotal[index] != testtotal[index - 1]:
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
        # print(name, sum(sizes))
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

    # changing the column orders so that mean length is after total no. of unfiltered pepts
    summary_table = summary_table[['enzyme', 'total peptides', 'mean length', 'median length', 'filtered peptides', 'coverage']]

    summary_table.to_csv(output_name, index=False)

#create_summary_table(unfiltered_rpg_files, filtered_rpg_files, fasta_files, 'method_summary.csv')


##### WITH RESIDUE FILTERING - doesn't work yet for parallel combinations
# FILTERED RPG FILES ABOVE NEED TO HAVE BEEN FILTERED FOR RESIDUE FOR THIS TO WORK
# COUNTING FREQ OF RESIDUES IN THE SEQUENCES

def add_residue_coverage_column(residue, fasta_files,summary_table, output_table_name, filtered_rpg_files):
    summary_table = pd.read_csv(summary_table)
    summary_table = summary_table.sort_values('enzyme')
    residue_freq_in_seq = 0
    for file in fasta_files:
        for rec in SeqIO.parse(file, "fasta"):
            freq_in_fasta = rec.seq.count(residue)
            residue_freq_in_seq += freq_in_fasta

    # COUNTING FREQ OF RESIDUES IN FILTERED DATA
    enzymelist = []
    freqlist = []
    freq_table = pd.DataFrame()
    for name in filtered_rpg_files:
        rpg_file = pd.read_csv(name)
        xlist = rpg_file['enzyme'].tolist()
        for enzyme in sorted(set(xlist)):
            temp_table = pd.DataFrame()
            positionslist = []
            parentisoforms = []
            rpg_by_enzyme = rpg_file[rpg_file.enzyme == enzyme]
            frequency = 0
            for parent, sequence, start in zip(rpg_by_enzyme['parent'], rpg_by_enzyme['sequence'], rpg_by_enzyme['peptide_start']):
                for resstartpos in re.finditer(residue, sequence):
                    respos = start + resstartpos.start()
                    positionslist.append(respos)
                    parentisoforms.append(parent)
                    # print(enzyme, resstartpos.start(), start, respos)
            temp_table['isoform'] = parentisoforms
            temp_table['residue_position'] = positionslist
            temp_table = temp_table.drop_duplicates()
            frequency = len(temp_table)
            freqlist.append(frequency)
            enzymelist.append(enzyme)
    freq_table['residue_freq'] = freqlist
    freq_table['enzyme'] = enzymelist
    freq_table = freq_table.sort_values(['enzyme']).reset_index(drop=True)
    freq_table['cumulative_freq'] = freq_table.groupby(['enzyme'])['residue_freq'].cumsum(axis=0)
    enzyme_table = freq_table.groupby(['enzyme']).tail(1)

    # ADDING COVERAGE COLUMN TO EXISTING TABLE
    enzyme_table['residue_coverage_%'] = (enzyme_table['cumulative_freq'] / residue_freq_in_seq) * 100
    summary_table[residue+' coverage %'] = enzyme_table['residue_coverage_%'].tolist()
    summary_table.to_csv(output_table_name, index=False)
