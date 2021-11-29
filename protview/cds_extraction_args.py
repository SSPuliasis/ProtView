import argparse
import cds_extraction_main as cds_extraction_main

parser = argparse.ArgumentParser(
    description="extract CDS information from a gff3 file", add_help=False)

parser.add_argument('input_file',
                    help='input file name of gff3 file')

def main():
    args = parser.parse_args()

    if args.input_file:
        input_name = args.input_file.replace('.gff3', '')
        print('extracting CDS information, input: {}.gff3'.format(input_name))
        cds_extraction_main.process_gff3(input_name)
        print('output: {}_+_cdsdf.csv, {}_-_cdsdf.csv'.format(input_name, input_name))

#print(args.input_file)