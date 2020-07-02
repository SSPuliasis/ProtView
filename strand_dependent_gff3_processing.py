import gffpandas.gffpandas as gffpd
import numpy as np


# CDS EXTRACTION FROM GFF FILE

#change this each time
input_gff_name = "chr1"

testprotein = gffpd.read_gff3(input_gff_name+".gff3")  # Read in the GFF3 file

cds_only = testprotein.filter_feature_of_type(['CDS']).to_gff3(
    input_gff_name+'_cds.gff3')  # extract only the CDS & write out to file

cds_only = gffpd.read_gff3(
    input_gff_name+"_cds.gff3")  # read in the original gff file and create a dataframe named cdsdf

cdsdf = cds_only.attributes_to_columns()
cdsdf = cdsdf.drop(['attributes'], axis=1)
cdsdf = cdsdf.drop('score', axis=1)
cdsdf = cdsdf.drop('phase', axis=1)

# unique ID for each CDS
# ASSIGNING ID TO EACH CDS
cdsdf['cds_id'] = cdsdf['seq_id']+'_'+cdsdf['start'].map(str)+'_'+cdsdf['end'].map(str)+'_'+cdsdf['strand']

##split off into separated dfs by strand here
posstrand_cdsdf = cdsdf[cdsdf['strand']=='+']
negstrand_cdsdf = cdsdf[cdsdf['strand']=='-']

#ordering by parent and start site, depending on strand
posstrand_cdsdf = posstrand_cdsdf.sort_values(by=['Parent', 'start'])  # ordering by parent and start site
negstrand_cdsdf = negstrand_cdsdf.sort_values(by = ['Parent', 'end'], ascending = False)

posstrand_cdsdf['prev_parent'] = posstrand_cdsdf['Parent'].shift()
posstrand_cdsdf['my_count'] = np.where(posstrand_cdsdf['prev_parent'] != posstrand_cdsdf['Parent'], 0, 5)  # 5 was chosen randomly, doesn't affect script
negstrand_cdsdf['prev_parent'] = negstrand_cdsdf['Parent'].shift()
negstrand_cdsdf['my_count'] = np.where(negstrand_cdsdf['prev_parent'] != negstrand_cdsdf['Parent'], 0, 5)

# ADDING PROTEIN SEQUENCE POSITIONS TO THE GFF TABLE
new_count = posstrand_cdsdf['my_count'].tolist()  # set 'no' as list
for index, og_count in enumerate(new_count):
    if og_count == 0:
        new_count[index] = 0
    elif og_count != 0:
        new_count[index] = new_count[index - 1] + 1
    # print(index,og_count,new_count)
posstrand_cdsdf['my_count'] = new_count

new_count = negstrand_cdsdf['my_count'].tolist()  # set 'no' as list
for index, og_count in enumerate(new_count):
    if og_count == 0:
        new_count[index] = 0
    elif og_count != 0:
        new_count[index] = new_count[index - 1] + 1
    # print(index,og_count,new_count)
negstrand_cdsdf['my_count'] = new_count

# end of prev cds = zero whenever no is zero
posstrand_cdsdf['end_of_previous'] = posstrand_cdsdf['end'].shift(periods=1)
posstrand_cdsdf = posstrand_cdsdf.fillna(0)
prev_end = posstrand_cdsdf['end_of_previous'].tolist()
for index, (count_value, end_value) in enumerate(zip(posstrand_cdsdf.my_count, prev_end)):
    if count_value == 0:
        prev_end[index] = 0
    # print(index, count_value ,prev_end)
posstrand_cdsdf['end_of_previous'] = list(map(int, prev_end))

negstrand_cdsdf['start_of_previous'] = negstrand_cdsdf['start'].shift(periods=1)
negstrand_cdsdf = negstrand_cdsdf.fillna(0)
prev_start = negstrand_cdsdf['start_of_previous'].tolist()
for index, (count_value, start_value) in enumerate(zip(negstrand_cdsdf.my_count, prev_start)):
    if count_value == 0:
        prev_start[index] = 0
    # print(index, count_value ,prev_start)
negstrand_cdsdf['start_of_previous'] = list(map(int, prev_start))


# Adding Intron ID column
posstrand_cdsdf["intron_id"] = posstrand_cdsdf["seq_id"] + "_" + (posstrand_cdsdf["end_of_previous"]+1).map(str) + "_" + (posstrand_cdsdf["start"]-1).map(str) + "_" + \
                     posstrand_cdsdf["strand"]
int_id = posstrand_cdsdf['intron_id'].tolist()
# make a new list to add values to
new_intron_id = []
for index, (count_value, id_value) in enumerate(zip(posstrand_cdsdf.my_count, int_id)):
    if count_value == 0:
        new_intron_id.append('na')
    else:
        new_intron_id.append(id_value)
posstrand_cdsdf['intron_id'] = new_intron_id

negstrand_cdsdf["intron_id"] = negstrand_cdsdf["seq_id"] + "_" + (negstrand_cdsdf["end"]+1).map(str) + "_" + (negstrand_cdsdf["start_of_previous"]-1).map(str) + "_" + \
                     negstrand_cdsdf["strand"]
int_id = negstrand_cdsdf['intron_id'].tolist()
# using new_count form before
# make a new list to add values to
new_intron_id = []
for index, (count_value, id_value) in enumerate(zip(negstrand_cdsdf.my_count, int_id)):
    if count_value == 0:
        new_intron_id.append('na')
    else:
        new_intron_id.append(id_value)
negstrand_cdsdf['intron_id'] = new_intron_id

# intron length & cumulative intron length
posstrand_cdsdf['intron_length'] = posstrand_cdsdf['start'] - posstrand_cdsdf[
    'end_of_previous']-1  # gives intron length as start for first cds of each isoform
intron_length = posstrand_cdsdf['intron_length'].tolist()
for index, (end, length) in enumerate(zip(posstrand_cdsdf.end_of_previous, intron_length)):
    if end == 0:
        intron_length[index] = 0
    # print(end,length,prev_end, intron_length)
posstrand_cdsdf['intron_length'] = intron_length
# can now drop end of previous column and prev parent
posstrand_cdsdf = posstrand_cdsdf.drop('end_of_previous', axis=1)
posstrand_cdsdf = posstrand_cdsdf.drop('prev_parent', axis=1)
# where value += previous value, unless value = 0
for index, length in enumerate(intron_length):
    if length == 0:
        intron_length[index] = 0
    else:
        intron_length[index] += intron_length[index - 1]
    # print(index,length,intron_length)
posstrand_cdsdf['cum_intron'] = intron_length

negstrand_cdsdf['intron_length'] = negstrand_cdsdf['start_of_previous'] - negstrand_cdsdf['end'] -1
# gives intron length = the end pos  for first cds of each isoform
intron_length = negstrand_cdsdf['intron_length'].tolist()
for index, (start, length) in enumerate(zip(negstrand_cdsdf.start_of_previous, intron_length)):
    if start == 0:
        intron_length[index] = 0
    # print(end,length,prev_end, intron_length)
negstrand_cdsdf['intron_length'] = intron_length
# can now drop end of previous column and prev parent
negstrand_cdsdf = negstrand_cdsdf.drop('start_of_previous', axis=1)
negstrand_cdsdf = negstrand_cdsdf.drop('prev_parent', axis=1)
# where value += previous value, unless value = 0
for index, length in enumerate(intron_length):
    if length == 0:
        intron_length[index] = 0
    else:
        intron_length[index] += intron_length[index - 1]
    # print(index,length,intron_length)
negstrand_cdsdf['cum_intron'] = intron_length

#column for start position of first cds for each isoform
first_start_list = []
for index, (count, start) in enumerate(zip(posstrand_cdsdf.my_count, posstrand_cdsdf.start)):
    if count == 0:
        first_start = start
        first_start_list.append(first_start)
    else:
        first_start_list.append(first_start)
posstrand_cdsdf['cds_first_start'] = first_start_list
posstrand_cdsdf['new_start'] = posstrand_cdsdf['start'] - posstrand_cdsdf['cds_first_start'] - posstrand_cdsdf['cum_intron']
posstrand_cdsdf['new_end'] = posstrand_cdsdf['end'] - posstrand_cdsdf['cds_first_start'] - posstrand_cdsdf['cum_intron']+1

first_end_list = []
for index, (count, end) in enumerate(zip(negstrand_cdsdf.my_count, negstrand_cdsdf.end)):
    if count == 0:
        first_end = end
        first_end_list.append(first_end)
    else:
        first_end_list.append(first_end)
negstrand_cdsdf['cds_first_end'] = first_end_list
negstrand_cdsdf['new_start'] = abs(negstrand_cdsdf['start'] - negstrand_cdsdf['cds_first_end'] + negstrand_cdsdf['cum_intron'])+1
negstrand_cdsdf['new_end'] = abs(negstrand_cdsdf['end'] - negstrand_cdsdf['cds_first_end'] + negstrand_cdsdf['cum_intron'])

# POSITIONS ON THE AMINO ACID SEQUENCE
# for the start sites
starts = posstrand_cdsdf['new_start']
startslist_dna = starts.tolist()
startslist_protein = []
position = 0
for position in startslist_dna:
    while position % 3 != 0:
        position = position + 1
    else:
        position = position / 3
        # print(position)
        startslist_protein.append(position)
posstrand_cdsdf['protein_start'] = startslist_protein
posstrand_cdsdf['protein_start'] = posstrand_cdsdf['protein_start'].astype(int)


starts = negstrand_cdsdf['new_start']
startslist_dna = starts.tolist()
startslist_protein = []
position = 0
for position in startslist_dna:
    while position % 3 != 0:
        position = position + 1
    else:
        position = position / 3
        # print(position)
        startslist_protein.append(position)
negstrand_cdsdf['protein_start'] = startslist_protein
negstrand_cdsdf['protein_start'] = negstrand_cdsdf['protein_start'].astype(int)

# for the end sites
ends = posstrand_cdsdf['new_end']
endslist_dna = ends.tolist()
endslist_protein = []
position = 0
for position in endslist_dna:
    while position % 3 != 0:
        position += 1
    else:
        position /= 3
        # print(position)
        endslist_protein.append(position)
posstrand_cdsdf['protein_end'] = endslist_protein
posstrand_cdsdf['protein_end'] = posstrand_cdsdf['protein_end'].astype(int)
posstrand_cdsdf = posstrand_cdsdf.drop('new_start', axis=1)  # remove columns that aren't necessary anymore
posstrand_cdsdf = posstrand_cdsdf.drop('new_end', axis=1)

ends = negstrand_cdsdf['new_end']
endslist_dna = ends.tolist()
endslist_protein = []
position = 0
for position in endslist_dna:
    while position % 3 != 0:
        position += 1
    else:
        position /= 3
        # print(position)
        endslist_protein.append(position)
negstrand_cdsdf['protein_end'] = endslist_protein
negstrand_cdsdf['protein_end'] = negstrand_cdsdf['protein_end'].astype(int)
negstrand_cdsdf = negstrand_cdsdf.drop('new_start', axis=1)  # remove columns that aren't necessary anymore
negstrand_cdsdf = negstrand_cdsdf.drop('new_end', axis=1)

# adding +1 to protein start positions to avoid overlap
posstrand_cdsdf['prev_prot_end'] = posstrand_cdsdf['protein_end'].shift()
protein_starts = posstrand_cdsdf['protein_start'].tolist()
for index, (start, prev_end) in enumerate(zip(protein_starts, posstrand_cdsdf.prev_prot_end)):
    if start == prev_end:
        protein_starts[index] +=1
    elif start ==0:
        protein_starts[index] +=1
    else:
        continue
posstrand_cdsdf['protein_start'] = protein_starts
posstrand_cdsdf = posstrand_cdsdf.drop('prev_prot_end', axis=1)

negstrand_cdsdf['prev_prot_start'] = negstrand_cdsdf['protein_start'].shift()
protein_ends = negstrand_cdsdf['protein_end'].tolist()
for index, (end, prev_start) in enumerate(zip(protein_ends, negstrand_cdsdf.prev_prot_start)):
    if end == prev_start:
        protein_ends[index] +=1
    elif end == 0:
        protein_ends[index] +=1
    else:
        continue
negstrand_cdsdf['protein_end'] = protein_ends
negstrand_cdsdf = negstrand_cdsdf.drop('prev_prot_start', axis=1)



#save the cdsdf
posstrand_cdsdf.to_csv(input_gff_name+'_+_cdsdf.csv')
negstrand_cdsdf.to_csv(input_gff_name+'_-_cdsdf.csv')