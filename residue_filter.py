# Filter peptides in the processed RPG output format (CSV) for peptides containing
# certain amino acids

#start by using cysteine as an example

#packages used
import pandas as pd


#import the rpg results/ peptide table

# read through and filter similarly to junction spanning
# write positive outcomes to a new file