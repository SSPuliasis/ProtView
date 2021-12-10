import argparse
import protview.fasta_input_sanitation as fasta_input_sanitation

parser = argparse.ArgumentParser(
    description="replaces '_' and '>' before carrying out a digest to avoid downstream issues", add_help=False)

parser.add_argument('fasta_file_name',
                    help='the fasta file containing the sequences')

def main():
    args = parser.parse_args()

    if args.fasta_file_name:
        print('replacing underscores in {}'.format(args.fasta_file_name))
        fasta_input_sanitation.replace_underscores_in_fasta(args.fasta_file_name)
        print('underscores replaced in {}'.format(args.fasta_file_name))
        print('replacing > in {}'.format(args.fasta_file_name))
        fasta_input_sanitation.replace_arrows_in_fasta(args.fasta_file_name)
        print('> replaced in {}'.format(args.fasta_file_name))