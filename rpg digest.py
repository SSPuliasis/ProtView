rpg -i AT1G67120.fasta -o AT1G67120_rpg.fasta -e 42 15 16 28 29 2 1

#To merge different results:
import shutil

newfile = 'chr5_rev_42_15_28_rpg.fasta'

firstfile = 'chr5_rev_42_rpg.fasta'
secondfile = 'chr5_rev_42_28_rpg.fasta'
thirdfile = 'chr5_rev_42_15_rpg.fasta'

with open(newfile, "wb") as wfd:
    for f in [firstfile, secondfile, thirdfile]:
        with open(f, "rb") as fd:
            shutil.copyfileobj(fd, wfd)

print("The content is merged successfully")