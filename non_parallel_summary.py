import pandas as pd
from Bio import SeqIO

unfiltered_rpg_files = ["trypsin_fwd_rpg_unfiltered.csv", "trypsin_rev_rpg_unfiltered.csv", "chym_fwd_rpg_unfiltered.csv", "chym_rev_rpg_unfiltered.csv", "lysc_fwd_rpg_unfiltered.csv", "lysc_rev_rpg_unfiltered.csv", "aspn_fwd_rpg_unfiltered.csv", "aspn_rev_rpg_unfiltered.csv"]
filtered_rpg_files = ["trypsin_fwd_rpg_filtered.csv", "trypsin_rev_rpg_filtered.csv", "chym_fwd_rpg_filtered.csv", "chym_rev_rpg_filtered.csv", "lysc_fwd_rpg_filtered.csv", "lysc_rev_rpg_filtered.csv", "aspn_fwd_rpg_filtered.csv", "aspn_rev_rpg_filtered.csv"]
fasta_files = ['araport11_proteome_fwd.fasta', 'araport11_proteome_rev.fasta']

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
    coverage_table = pd.DataFrame(columns=['input file', 'enzyme', 'peptide length', 'protein length', 'coverage'])
    temp_table = pd.DataFrame(columns=['input file', 'enzyme', 'peptide length', 'protein length'])
    protein_len_sum = 0
    for name in fasta_files:
        sizes = [len(rec) for rec in SeqIO.parse(name, "fasta")]
        print(name, sum(sizes))
        protein_len_sum += sum(sizes)
    for name in filtered_rpg_files:
        rpg_file = pd.read_csv(name)
        a = pd.DataFrame(rpg_file.groupby(['enzyme'], as_index=False)['peptide_size'].sum())
        pep_len = a['peptide_size'].tolist()
        enzymes = a['enzyme'].tolist()
        temp_table['peptide length'] = pep_len
        temp_table['enzyme'] = enzymes
        temp_table['protein length'] = protein_len_sum
        temp_table['input file'] = name
        coverage_table = coverage_table.append(temp_table, sort=True)
    coverage_table['coverage'] = (coverage_table['peptide length'] / coverage_table['protein length']) * 100
    coverage_summary_table = pd.DataFrame(coverage_table.groupby(['enzyme'], as_index=False)['peptide length'].sum())
    total_protein_length = sum(set(coverage_table['protein length'].tolist()))
    coverage_summary_table['protein length'] = total_protein_length
    coverage_summary_table['coverage'] = (coverage_summary_table['peptide length'] / coverage_summary_table[
        'protein length']) * 100
    # write.csv(coverage_summary_table, 'coverage_summary_table.csv')
    summary_table['coverage'] = coverage_summary_table['coverage']
    # MEAN PEPTIDE LENGTHS
    summary_table['mean length'] = coverage_summary_table['protein length'] / summary_table['total peptides']
    # summary_table
    # changing the column orders so that mena length is after total no. of unfiltered pepts
    summary_table = summary_table[['enzyme', 'total peptides', 'mean length', 'filtered peptides', 'coverage']]

    summary_table.to_csv(output_name)


#create_summary_table(unfiltered_rpg_files, filtered_rpg_files, fasta_files, 'single_summary.csv')