import argparse
import sys
import checkv


def cli():
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""CheckV v{checkv.__version__}: assessing the quality of metagenome-assembled viral genomes
https://bitbucket.org/berkeleylab/checkv

usage: checkv <program> [options]

programs:
    end_to_end          run full pipeline to estimate completeness, contamination, and identify closed genomes
    contamination       identify and remove host contamination on integrated proviruses
    completeness        estimate completeness for genome fragments
    complete_genomes    identify complete genomes based on terminal repeats and flanking host regions
    quality_summary     summarize results across modules
    download_database   download the latest version of CheckV's database
    update_database     update CheckV's database with your own complete genomes""",
    )

    subparsers = parser.add_subparsers(help=argparse.SUPPRESS)

    end_to_end_parser = subparsers.add_parser(
        "end_to_end",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Run full pipeline to estimate completeness, contamination, and identify closed genomes
\nusage: checkv end_to_end <input> <output> [options]""",
    )
    checkv.end_to_end.fetch_arguments(end_to_end_parser)

    download_database_parser = subparsers.add_parser(
        "download_database",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Download the latest version of CheckV's database
\nusage: checkv download_database <destination>""",
    )
    checkv.download_database.fetch_arguments(download_database_parser)

    update_database_parser = subparsers.add_parser(
        "update_database",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Update CheckV's database with your own complete genomes
\nusage: checkv update_database <source_db> <dest_db> <genomes> [options]""",
    )
    checkv.update_database.fetch_arguments(update_database_parser)

    contamination_parser = subparsers.add_parser(
        "contamination",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Estimate host contamination for integrated proviruses
\nusage: checkv contamination <input> <output> [options]""",
    )
    checkv.contamination.fetch_arguments(contamination_parser)

    completeness_parser = subparsers.add_parser(
        "completeness",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Estimate completeness for genome fragments
\nusage: checkv completeness <input> <output> [options]""",
    )
    checkv.completeness.fetch_arguments(completeness_parser)

    complete_genomes_parser = subparsers.add_parser(
        "complete_genomes",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Identify complete genomes based on terminal repeats and flanking host regions
\nusage: checkv complete_genomes <input> <output> [options]""",
    )
    checkv.complete_genomes.fetch_arguments(complete_genomes_parser)

    summary_parser = subparsers.add_parser(
        "quality_summary",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Summarize results across modules
\nusage: checkv quality_summary <input> <output> [options]""",
    )
    checkv.quality_summary.fetch_arguments(summary_parser)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "download_database":
            download_database_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "update_database":
            update_database_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "complete_genomes":
            complete_genomes_parser.print_help()
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
        elif sys.argv[1] == "end_to_end":
            end_to_end_parser.print_help()
            sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)
