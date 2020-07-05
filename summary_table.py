import pandas as pd
from Bio import SeqIO

# enter the names of all of the input files as a list
rpg_files = ['chr1_fwd', 'chr1_rev', 'chr2_fwd', 'chr2_rev', 'chr3_fwd', 'chr3_rev', 'chr4_fwd', 'chr4_rev',
               'chr5_fwd', 'chr5_rev']
#  '_rpg_filtered.csv' or '_rpg_unfiltered.csv' can then be added

#  TOTAL PEPTIDES GENERATED - make into total_peptides_generated() function
testname = []
testenzyme = []
testtotal = []
for name in rpg_files:
    input = pd.read_csv(name + "_rpg_unfiltered.csv")
    x = input['enzyme'].tolist()
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
temp_table['next_enzyme'] = temp_table['enzyme'].shift(periods = -1)
for index, (a, b, c) in enumerate(zip(temp_table.enzyme, temp_table.next_enzyme, temp_table.cumulative_totals)):
    if a != b:
        enzymelist.append(a)
        totallist.append(c)
    else:
        continue

# creating the summary table
summary_table = pd.DataFrame()
summary_table['enzyme'] = enzymelist
summary_table['total peptides'] = totallist
summary_table

#  TOTAL FILTERED PEPTIDES GENERATED
testname = []
testenzyme = []
testtotal = []
for name in rpg_files:
    input = pd.read_csv(name + "_rpg_filtered.csv")
    x = input['enzyme'].tolist()
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
temp_table['next_enzyme'] = temp_table['enzyme'].shift(periods = -1)
for index, (a, b, c) in enumerate(zip(temp_table.enzyme, temp_table.next_enzyme, temp_table.cumulative_totals)):
    if a != b:
        totallist.append(c)
    else:
        continue
summary_table['filtered peptides'] = totallist

coverage_table = pd.DataFrame(columns = ['input file', 'enzyme', 'peptide length', 'protein length', 'coverage'])
temp_table = pd.DataFrame(columns = ['input file', 'enzyme', 'peptide length', 'protein length'])

for name in rpg_files:
    sizes = [len(rec) for rec in SeqIO.parse(name+".fasta", "fasta")]
    rpg_file = pd.read_csv(name+'_rpg_filtered.csv')
    input_name = name+'_rpg_filtered'
    a = pd.DataFrame(rpg_file.groupby(['enzyme'], as_index=False)['peptide_size'].sum())
    pep_len = a['peptide_size'].tolist()
    enzymes = a['enzyme'].tolist()
    temp_table['peptide length'] = pep_len
    temp_table['enzyme'] = enzymes
    temp_table['protein length'] = sum(sizes)
    temp_table['input file'] = input_name
    coverage_table = coverage_table.append(temp_table, sort=True)
coverage_table['coverage'] = (coverage_table['peptide length']/ coverage_table['protein length']) *100
# write.csv(coverage_table, 'coverage_table.csv')

#coverage_table

coverage_summary_table = pd.DataFrame(coverage_table.groupby(['enzyme'], as_index=False) ['peptide length'].sum())
total_protein_length = sum(set(coverage_table['protein length'].tolist())) # why use set????
coverage_summary_table['protein length'] = total_protein_length
coverage_summary_table['coverage'] = (coverage_summary_table['peptide length']/ coverage_summary_table['protein length']) *100
#write.csv(coverage_summary_table, 'coverage_summary_table.csv')

summary_table['coverage'] = coverage_summary_table['coverage']
#summary_table

#MEAN PEPTIDE LENGTHS
summary_table['mean length'] = coverage_summary_table['protein length']/summary_table['total peptides']
#summary_table

#changing the column orders so that mena length is after total no. of unfiltered pepts
summary_table = summary_table[['enzyme', 'total peptides', 'mean length', 'filtered peptides', 'coverage']]

summary_table.to_csv('summary_table.csv')