import argparse
import junction_spanning_script

if __name__ == '__main__':
    print("junction_spanning has been imported")

    parser = argparse.ArgumentParser(
        description="filter for peptides that cover exon-exon junctions")

    parser.add_argument('rpg_file',
                        help='csv file containing the peptides to be filtered')
    parser.add_argument('cds_files', nargs =2,
                        help='both csv files containing extracted coding sequences for each DNA strand, ending in +/-_cdsdf.csv')
    parser.add_argument('output_name',
                        help='desired name of output csv file, ProtView will automatically differentiate between DNA strands')

    args = parser.parse_args()

    if args.rpg_file and args.cds_files and args.output_name:
        print('filtering {} for junction covering peptides'.format(args.rpg_file))
        junction_spanning_script.junction_spanning(args.cds_files, args.rpg_file, args.output_name)
        print('output: +_{}, -_{}'.format(args.output_name, args.output_name))


