def replace_underscores_in_fasta(fasta_file):
    filein = open(fasta_file, 'r')
    data = filein.read()
    filein.close
    data = data.replace('_', '-')
    filein.close()
    filein = open(fasta_file, 'wt')
    filein.write(data)
    filein.close()