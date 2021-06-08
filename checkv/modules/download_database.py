import os
import shutil
import time
import urllib.request
import checkv
from checkv import utility


class DatabaseDownloader:
    def __init__(self, destination):
        self.url = "https://portal.nersc.gov/CheckV/"
        self.destination = destination
        self.version = (
            urllib.request.urlopen(self.url + "CURRENT_RELEASE.txt")
            .read()
            .decode("utf-8")
            .strip()
        )
        self.filename = self.version + ".tar.gz"
        self.output_file = os.path.join(self.destination, self.filename)

    def download(self):
        database_url = self.url + self.filename
        with urllib.request.urlopen(database_url) as response:
            with open(self.output_file, "wb") as fout:
                shutil.copyfileobj(response, fout)

    def extract(self):
        shutil.unpack_archive(self.output_file, self.destination, "gztar")
        os.remove(self.output_file)


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="download_database")
    parser.add_argument(
        "destination",
        type=str,
        help="Directory where the database will be downloaded to.",
    )
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Suppress logging messages",
    )


def main(args):
    program_start = time.time()
    logger = utility.get_logger(args["quiet"])
    if not os.path.exists(args["destination"]):
        os.makedirs(args["destination"])

    logger.info(f"\nCheckV v{checkv.__version__}: download_database")

    logger.info("[1/3] Checking latest version of CheckV's database...")
    db = DatabaseDownloader(args["destination"])

    logger.info(f"[2/3] Downloading '{db.version}'...")
    db.download()

    logger.info(f"[3/3] Extracting '{db.version}'...")
    db.extract()

    logger.info("Run time: %s seconds" % round(time.time() - program_start, 2))
    logger.info("Peak mem: %s GB" % round(utility.max_mem_usage(), 2))
