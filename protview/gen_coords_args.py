import argparse
import modules.gen_coords_main as gen_coords_main

parser = argparse.ArgumentParser(
    description="filter for peptides that cover exon-exon junctions", add_help=False)

parser.add_argument('rpg_file',
                    help='csv file containing the peptides and their proteomic coorindates')
parser.add_argument('cds_files', nargs=2,
                    help='both csv files containing extracted coding sequences for each DNA strand, ending in +/-_cdsdf.csv')
def main():
    args = parser.parse_args()

    if args.rpg_file and args.cds_files:
        print('calculating genomic coordinates for peptides in {}'.format(args.rpg_file))
        gen_coords_main.calculate_gen_coords(args.rpg_file, args.cds_files)
        file_name = args.rpg_file.replace('.csv', '')
        print('output: {}_+_positions_df.csv, {}_-_positions_df.csv'.format(file_name, file_name))