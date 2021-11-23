import argparse
import protview.unique_count_main as unique_count_main

parser = argparse.ArgumentParser(
    description='counts the number of peptides per enzyme that are only generated for one isoform per digest and '
                'adds column onto an existing statistics table', add_help=False)

parser.add_argument('table_name',
                    help='the general or junction csv summary table for the column to be appended to')
parser.add_argument('rpg_files', nargs='+',
                    help='the filtered digest results in csv format')
def main():
    args = parser.parse_args()

    if args.table_name and args.rpg_files:
        print('counting isoform-unique peptides in {}'.format(args.rpg_files))
        unique_count_main.add_unique_pept_column(args.table_name, args.rpg_files)
        print('isoform-unique peptides column has been added to {}'.format(args.table_name))


