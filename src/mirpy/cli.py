import argparse
from .tester import miRpyTest
from .downloader import download_mirbase_gff

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

    # Support for downloading the latest miRBase GFF3 where user
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

    return p

if __name__ == "__main__":
    raise SystemExit(main())
