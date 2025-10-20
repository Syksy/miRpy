import argparse
from .tester import miRpyTest

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "test":
        return miRpyTest(args.bams, n=args.num, region=args.region)

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

    # Support for downloading the latest miRBase GFF3 where user
    d = sub.add_parser("download", help="Download miRBase hairpin GFF3 to a path.")
    d.add_argument(
        "dest",
        help="Destination file path, e.g. data/hairpin.gff3.gz (or .gff3 to auto-decompress).",
    )
    d.add_argument(
        "--url",
        default=None,
        help="Optional alternate URL (defaults to miRBase CURRENT hairpin.gff3.gz).",
    )

    return p

if __name__ == "__main__":
    raise SystemExit(main())
