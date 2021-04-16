import argparse
import cds_extraction

if __name__ == '__main__':
    print("strand_dependent_gff3_processing has been imported")

    parser = argparse.ArgumentParser(
        description="extract CDS information from a gff3 file")

    parser.add_argument('input_file',
                        help='input file name of gff3 file')

    args = parser.parse_args()

    if args.input_file:
        input_name = args.input_file.replace('.gff3', '')
        print('extracting CDS information, input: {}.gff3'.format(input_name))
        cds_extraction.process_gff3(input_name)
        print('output: {}_+_cdsdf.csv, {}_-_cdsdf.csv'.format(input_name, input_name))
