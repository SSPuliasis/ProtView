# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 12:53:02 2021

@author: SSPuliasis
"""
import argparse
import modules.rpg_output_processing_main as rpg_output_processing_main

parser = argparse.ArgumentParser(
        description="process peptides generated in silico by RPG", add_help =False)

parser.add_argument('input_file',
                    help='input file name of fasta file')
parser.add_argument('-mc', '--miscleavage', help='miscleavage value, default =0',
                    metavar='', type=int, default=0)
parser.add_argument('-e', '--enzymes', nargs='+',
                    help='enzymes to create parallel digest with',
                    metavar='')
parser.add_argument('-r', '--residue', help='residue to be filtered for',
                    metavar='')
parser.add_argument('-min', '--minlen',
                    help='minimum peptide length to filter for, default = 7',
                    type=int, metavar='', default=7)
parser.add_argument('-max', '--maxlen',
                    help='maximum peptide length to filter for, default = 35',
                    type=int, metavar='', default=35)

def main():
    args = parser.parse_args()

    print("processing peptides")
    #process_rpg_output
    input_name = args.input_file.replace('.fasta', '')

    output_file_name_csv = str("{}.csv".format(args.input_file))
    print("process_rpg_output - input:{}, output:{}".format(args.input_file,
          output_file_name_csv))
    recent_file_name = str(input_name)
    recent_file_name_csv = str("{}.csv".format(input_name))
    rpg_output_processing_main.process_rpg_output(recent_file_name)

    #only create MC is a value is provided
    if args.miscleavage > 0:
        output_file_name = str("{}_mc{}".format(input_name, args.miscleavage))
        output_file_name_csv = str("{}.csv".format(output_file_name))
        #rpg_output_processing_main.create_miscleavage(recent_file_name_csv, args.miscleavage, output_name_csv)
        print("creating miscleavage <= {} - input:{}, output:{}".format(
                args.miscleavage, recent_file_name_csv, output_file_name_csv))
        rpg_output_processing_main.create_miscleavage(recent_file_name_csv,
                                                 args.miscleavage,
                                                 output_file_name_csv)
        recent_file_name = output_file_name
        recent_file_name_csv = output_file_name_csv
    else:
        print("no-miscleavage")

    #create parallel
    if args.enzymes:
        output_file_name = str("{}_parallel".format(recent_file_name))#correct
        output_file_name_csv = str("{}.csv".format(output_file_name))#correct
        print("creating a parallel digest with {}, input:{}, output:{}.".format(
                args.enzymes, recent_file_name_csv, output_file_name_csv))
        # for enz in args.enzymes:
            # print(enz)
        # print(type(args.enzymes))
        rpg_output_processing_main.create_parallel_digest(recent_file_name_csv,
                                                     output_file_name_csv,
                                                     args.enzymes)
        recent_file_name = output_file_name
        recent_file_name_csv = output_file_name_csv

    #filter by residue by min & max lens (defaults if necessary)
    if args.minlen and args.maxlen:
        output_file_name = str("{}_len_{}_{}".format(recent_file_name, args.minlen,
                               args.maxlen))
        output_file_name_csv = str("{}.csv".format(output_file_name))
        print("filtering for peptide length of {}-{} amino acids, input:{}, output:{}".format(
                args.minlen, args.maxlen, recent_file_name_csv, output_file_name_csv))
        rpg_output_processing_main.filter_by_length(args.minlen, args.maxlen,
                                               recent_file_name_csv,
                                               output_file_name_csv)
        recent_file_name = output_file_name
        recent_file_name_csv = output_file_name_csv

    #filter for residue - this is an extra option
    if args.residue:
        output_file_name = str("{}_{}".format(recent_file_name, args.residue))
        output_file_name_csv = str("{}.csv".format(output_file_name))
        print("filtering for peptides containing {}, input:{}, output:{}".format(
                args.residue, recent_file_name_csv, output_file_name_csv))
        rpg_output_processing_main.filter_for_residue(args.residue, recent_file_name_csv,
                                                 output_file_name_csv)

    #merge files - should probably be moved


