import argparse
import junction_summary_stats

if __name__ == '__main__':
    #print("junction_summary_stats has been imported")

    parser = argparse.ArgumentParser(
        description="summary statistics for peptides covering exon-exon junctions")

    parser.add_argument('-pept', '--junction_spanning_peptides', metavar='', nargs='+',
                        help='csv files containing the peptides to be filtered')
    parser.add_argument('-cds', '--cds_files', metavar='', nargs='+',
                        help='both csv files containing extracted coding sequences for each DNA strand, ending in +/-_cdsdf.csv')
    parser.add_argument('-out', '--output_name', default='junction_summary.csv', metavar='',
                        help='desired name of output csv summary file, default: junction_summary.csv')

    args = parser.parse_args()

    if args.junction_spanning_peptides and args.cds_files and args.output_name:
        print('calculating junction summary statistics for {}'.format(args.junction_spanning_peptides))
        junction_summary_stats.junction_summary_stats(args.junction_spanning_peptides, args.cds_files, args.output_name)
        print('output: {}'.format(args.output_name))

