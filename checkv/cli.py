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

    circularity_parser = subparsers.add_parser(
        "circularity",
        help="identify circular genomes.",
        description="Identify circular genomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.dtr.fetch_arguments(circularity_parser)


    completeness_parser = subparsers.add_parser(
        "completeness",
        help="estimate completeness for genome fragments.",
        description="Estimate completeness for genome fragments.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.completeness.fetch_arguments(completeness_parser)
    contamination_parser = subparsers.add_parser(
        "contamination",
        help="estimate contamination for integrated prophages.",
        description="Estimate contamination for integrated prophages.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    checkv.contamination.fetch_arguments(contamination_parser)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "circularity":
            circularity_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "completeness":
            completeness_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "contamination":
            contamination_parser.print_help()
            sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)
