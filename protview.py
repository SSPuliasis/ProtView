import argparse
import sys

import modules.cds_extraction_args as cds_extraction
import modules.gen_coords_args as gen_coords
import modules.junction_peptides_args as junction_peptides
import modules.junction_summary_args as junction_summary
import modules.rpg_output_processing_args as rpg_output_processing
import modules.summary_stats_args as summary_stats
import modules.unique_count_args as unique_count

# make a parser in this with help options etc
parser = argparse.ArgumentParser(prog='ProtView',
                                 description="ProtView is a a tool for comparing proteases and their combinations \n"
                                             "in silico, in a proteomic and transcriptomic context")

subparsers = parser.add_subparsers()

# cds_extraction parser
cds_extraction_subparser = subparsers.add_parser("cds_extraction",
                                                 parents=[cds_extraction.parser])
cds_extraction_subparser.set_defaults(func="cds_extraction")

# gen coordinate parser
gen_coords_subparser = subparsers.add_parser("gen_coords",
                                             parents=[gen_coords.parser])
gen_coords_subparser.set_defaults(func="gen_coords")

# junction peptides parser
junction_peptides_subparser = subparsers.add_parser("junction_peptides",
                                                    parents=[junction_peptides.parser])
junction_peptides_subparser.set_defaults(func="junction_peptides")

# junction summary parser
junction_summary_subparser = subparsers.add_parser("junction_summary",
                                                   parents=[junction_summary.parser])
junction_summary_subparser.set_defaults(func="junction_summary")

# rpg_output_processing parser
rpg_output_processing_subparser = subparsers.add_parser("rpg_output_processing",
                                                        parents=[rpg_output_processing.parser])
rpg_output_processing_subparser.set_defaults(func="rpg_output_processing")

# summary stats parser
summary_stats_subparser = subparsers.add_parser("summary_stats",
                                                parents=[summary_stats.parser])
summary_stats_subparser.set_defaults(func="summary_stats")

# unique count parser
unique_count_subparser = subparsers.add_parser("unique_count",
                                               parents=[unique_count.parser])
unique_count_subparser.set_defaults(func="unique_count")


def main():
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if args.func == "cds_extraction":
        cds_extraction.parser = parser
        cds_extraction.main()
    elif args.func == "gen_coords":
        gen_coords.parser = parser
        gen_coords.main()
    elif args.func == "junction_peptides":
        junction_peptides.parser = parser
        junction_peptides.main()
    elif args.func == "junction_summary":
        junction_summary.parser = parser
        junction_summary.main()
    elif args.func == "rpg_output_processing":
        rpg_output_processing.parser = parser
        rpg_output_processing.main()
    elif args.func == "summary_stats":
        summary_stats.parser = parser
        summary_stats.main()
    elif args.func == "unique_count":
        unique_count.parser = parser
        unique_count.main()


if __name__ == '__main__':
    main()
