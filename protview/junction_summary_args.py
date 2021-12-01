import argparse
import protview.junction_summary_main as junction_summary_main

parser = argparse.ArgumentParser(
    description="summary statistics for peptides covering exon-exon junctions", add_help=False)

parser.add_argument('-pept', '--junction_spanning_peptides', metavar='', nargs='+',
                    help='csv files containing the peptides to be filtered')
parser.add_argument('-cds', '--cds_files', metavar='', nargs='+',
                    help='both csv files containing extracted coding sequences for each DNA strand, ending in +/-_cdsdf.csv')
parser.add_argument('-out', '--output_name', default='junction_summary.csv', metavar='',
                    help='desired name of output csv summary file, default: junction_summary.csv')

def main():
    args = parser.parse_args()

    if args.junction_spanning_peptides and args.cds_files and args.output_name:
        print('calculating junction summary statistics for {}'.format(args.junction_spanning_peptides))
        junction_summary_main.junction_summary_stats(args.junction_spanning_peptides, args.cds_files, args.output_name)
        print('output: {}'.format(args.output_name))

