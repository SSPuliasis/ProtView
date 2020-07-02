#Replacing missing "NA" values in the PeptideCutter output table


#set working directory to where the files are located, this is also where the new file will be saved

#The tidyr package is required
install.packages("tidyr")

#If tidyr is already installled, start here
library(tidyr)

#Importing and viewing the PeptideCutter results table
library(readxl)
#Change the excel file name in the script below to the file you want to import (& directory - if not in cluster folder)
PC_digest <- read_excel("AT1G67120_PC.xlsx")
View(PC_digest)
attach(PC_digest)

#Filling in the gaps
filled_in <- fill(PC_digest, `Position of cleavage site`:`Peptide mass [Da]`, .direction=c("down"))
View(filled_in)

#changing the column headers
colnames(AT1G67120_PC_table) <- c("cleavage_position", "enzyme", "sequence", "peptide_size", "mol_weight")

#Write the new table out to a CSV file
#The name of the CSV file needs to be changed every time
write.csv(filled_in, "AT1G67120_PC_table.csv")

detach(PC_digest)
