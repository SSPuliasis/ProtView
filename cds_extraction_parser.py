import argparse
import cds_extraction

if __name__ == '__main__':
    print("strand_dependent_gff3_processing has been imported")

    parser = argparse.ArgumentParser(
        description="extract CDS information from a gff3 file")

    parser.add_argument('input_file',
                        help='input file name of gff3 file excluding the ".gff3" extension')

    args = parser.parse_args()

    if args.input_file:
        print('extracting CDS information, input: {}.gff3'.format(args.input_file))
        cds_extraction.process_gff3(args.input_file)
        print('output: {}_+_cdsdf.csv, {}_-_cdsdf.csv'.format(args.input_file, args.input_file))
