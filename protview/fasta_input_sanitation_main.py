from Bio import SeqIO
import os

def replace_underscores_in_fasta(fasta_file):
    filein = open(fasta_file, 'r')
    data = filein.read()
    filein.close
    data = data.replace('_', '-')
    filein.close()
    filein = open(fasta_file, 'wt')
    filein.write(data)
    filein.close()

def replace_arrows_in_fasta(fasta_file):
    input_name = fasta_file
    with open('corrected.fasta', 'w') as outfile:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            seq.description = seq.description.replace('>', '-')
            SeqIO.write(seq, outfile, 'fasta')
    os.remove(fasta_file)
    os.rename('corrected.fasta', input_name)