import argparse
from .tester import miRpyTest
from .downloader import download_mirbase_gff
from .gfftools import subset_gff_by_criteria

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "test":
        return miRpyTest(args.bams, n=args.num, region=args.region)

    if args.cmd == "download":
        url = args.url or None
        try:
            out_path = download_mirbase_gff(args.dest, url=url or None)
            print(f"Downloaded miRBase GFF to: {out_path}")
            return 0
        except FileExistsError as e:
            print(f"[ERROR] {e}")
            return 2
        except Exception as e:
            print(f"[ERROR] {e}")
            return 1

    parser.error("Unknown command")
    return 2


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mirpy",
        description="Count mature miRNA with bamnostic."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # Sanity checking of BAM files
    t = sub.add_parser(
        "test",
        help="Print first N reads from each BAM file (bamnostic-based)."
    )
    t.add_argument("bams", nargs="+", help="One or more BAM files to test.")
    t.add_argument("-n", "--num", type=int, default=10, help="Number of reads per BAM.")
    t.add_argument("--region", help="Optional region like 'chr1:1000-2000'.")

    # Support for downloading the latest miRBase GFF3
    d = sub.add_parser("download", help="Download miRBase GFF3 to a path.")
    d.add_argument(
        "dest",
        help="Destination file path, e.g. data/hsa.gff3.",
    )
    d.add_argument(
        "--url",
        default=None,
        help="Optional alternate URL (defaults to miRBase CURRENT hsa.gff3).",
    )

    # Subsetting and processing GFF3 (tab-separated miRNA-annotations)
    g = sub.add_parser(
        "gff-subset",
        help="Subset a GFF3 file based on a column-matching criterion (default: keep only miRNA features).",
    )
    g.add_argument(
        "--in",
        dest="in_path",
        required=True,
        help="Input GFF3 (.gff3 or .gff3.gz).",
    )
    g.add_argument(
        "--out",
        dest="out_path",
        required=True,
        help="Output GFF3 (.gff3 or .gff3.gz).",
    )
    g.add_argument(
        "--col",
        dest="col_index",
        type=int,
        default=2,
        help="Zero-based column index to test (default=2 for the 3rd column).",
    )
    g.add_argument(
        "--criteria",
        nargs="+",
        default=["miRNA"],
        help="Value(s) to keep in the selected column. Default: 'miRNA'.",
    )
    g.add_argument(
        "--split-char",
        dest="split_char",
        default="\t",
        help="Character to split columns on (default: tab).",
    )
    g.add_argument(
        "--case-insensitive",
        dest="case_insensitive",
        action="store_true",
        help="Ignore case when matching criteria.",
    )


    return p

if __name__ == "__main__":
    raise SystemExit(main())
