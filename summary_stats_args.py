import argparse
import summary_table

if __name__ == '__main__':
    #print("summary_table funtions have been imported")

    parser = argparse.ArgumentParser(
        description="generates a table of summary statistics")

    parser.add_argument('-fasta', '--fasta_files', nargs='+',
                        help='fasta sequence files that generated the peptides', metavar='')
    parser.add_argument('-u', '--unfiltered_rpg_files', nargs='+',
                        help='unfiltered rpg peptides in csv format', metavar='')
    parser.add_argument('-f', '--filtered_rpg_files', nargs='+',
                        help='filtered rpg peptides in csv format', metavar='')
    parser.add_argument('-o', '--output_csv_name',
                        help='name of the output file, ending in .csv. Default = summary_table.csv',
                        metavar='', default='summary_table.csv')
    parser.add_argument('-r', '--residues',
                        help='residues to calculate the coverage of, eg C',
                        metavar='', nargs='+')

    args = parser.parse_args()

    if args.fasta_files and args.unfiltered_rpg_files and args.filtered_rpg_files:
        print('generating summary statistics, input fasta: {}, unfiltered peptides: {}, filtered peptides: {}'.format(
            args.fasta_files, args.unfiltered_rpg_files, args.filtered_rpg_files))
        summary_table.create_summary_table(args.unfiltered_rpg_files, args.filtered_rpg_files, args.fasta_files,
                                           args.output_csv_name)
        if args.residues:
            for res in args.residues:
                summary_table.add_residue_coverage_column(res, args.fasta_files, args.output_csv_name,
                                                          args.output_csv_name, args.filtered_rpg_files)
        print('output: {}'.format(args.output_csv_name))

