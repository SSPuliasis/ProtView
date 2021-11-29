import argparse
import replace_underscores_main as replace_underscores_main

parser = argparse.ArgumentParser(
    description="replaces '_' with '-' before carrying out a digest", add_help=False)

parser.add_argument('fasta_file_name',
                    help='the fasta file containing the sequences')

def main():
    args = parser.parse_args()

    if args.fasta_file_name:
        print('replacing underscores in {}'.format(args.fasta_file_name))
        replace_underscores_main.replace_underscores_in_fasta(args.fasta_file_name)
        print('underscores replaced in {}'.format(args.fasta_file_name))