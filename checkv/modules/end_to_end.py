import os
import shutil
import checkv

def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="end_to_end")
    parser.add_argument(
        "input", type=str, help="Input nucleotide sequences in FASTA format"
    )
    parser.add_argument("output", type=str, help="Output directory")
    parser.add_argument(
        "-d",
        dest="db",
        type=str,
        required=False,
        metavar="PATH",
        help="Reference database path. By default the CHECKVDB environment variable is used",
    )
    parser.add_argument(
        "--remove_tmp",
        action="store_true",
        default=False,
        help="Delete intermediate files from the output directory",
    )
    parser.add_argument(
        "-t",
        dest="threads",
        type=int,
        default=1,
        metavar="INT",
        help="Number of threads to use for Prodigal and DIAMOND",
    )
    parser.add_argument(
        "--restart",
        action="store_true",
        default=False,
        help="Overwrite existing intermediate files. By default CheckV continues where program left off",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def main(args):

    # if restart flag, delete tmp only 1x
    args["tmp"] = os.path.join(args["output"], "tmp")
    if args["restart"] and os.path.exists(args["tmp"]):
        shutil.rmtree(args["tmp"])
    args["restart"] = False

    checkv.contamination.main(args)
    checkv.completeness.main(args)
    checkv.complete_genomes.main(args)
    checkv.quality_summary.main(args)

    if args["remove_tmp"]:
        shutil.rmtree(args["tmp"])
