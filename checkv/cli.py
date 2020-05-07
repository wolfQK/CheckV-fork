import argparse
import sys

import checkv


def cli():
    parser = argparse.ArgumentParser(
        description="Assess the quality of viral genomes from metagenomes & viromes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"{parser.prog} v{checkv.__version__}"
    )
    subparsers = parser.add_subparsers()

    contamination_parser = subparsers.add_parser(
        "contamination",
        help="estimate host contamination for integrated proviruses",
        description="Estimate host contamination for integrated proviruses",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.contamination.fetch_arguments(contamination_parser)

    completeness_parser = subparsers.add_parser(
        "completeness",
        help="estimate completeness for genome fragments",
        description="Estimate completeness for genome fragments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.completeness.fetch_arguments(completeness_parser)

    repeats_parser = subparsers.add_parser(
        "repeats",
        help="identify terminal repeats & estimate genome copy #",
        description="Identify terminal repeats & estimate genome copy #",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.repeats.fetch_arguments(repeats_parser)

    summary_parser = subparsers.add_parser(
        "quality_summary",
        help="summarize results across modules",
        description="Summarize results across modules",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.quality_summary.fetch_arguments(summary_parser)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "repeats":
            repeats_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "completeness":
            completeness_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "contamination":
            contamination_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "quality_summary":
            summary_parser.print_help()
            sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)
