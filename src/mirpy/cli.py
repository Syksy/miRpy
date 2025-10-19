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

    t = sub.add_parser(
        "test",
        help="Print first N reads from each BAM file (bamnostic-based)."
    )
    t.add_argument("bams", nargs="+", help="One or more BAM files to test.")
    t.add_argument("-n", "--num", type=int, default=10, help="Number of reads per BAM.")
    t.add_argument("--region", help="Optional region like 'chr1:1000-2000'.")

    return p

if __name__ == "__main__":
    raise SystemExit(main())
