import argparse
import importlib.resources as res
import shutil

from .viewer import view_bam_head
from .downloader import download_mirbase_gff
from .gfftools import subset_gff_by_criteria
from .count import count_matrix

from pathlib import Path

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Test by taking the top-view (head) of BAM files
    if args.cmd in ["view", "head"]:
        return view_bam_head(args.bams, n=args.num, region=args.region)

    # Download miRBase GFF3 annotations
    elif args.cmd == "download":
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

    # Handling GFF3 subsetting
    elif args.cmd in ["subset", "gff-subset"]:
        try:
            n = subset_gff_by_criteria(
                in_path=args.in_path,
                out_path=args.out_path,
                col_index=args.col_index,
                split_char=args.split_char,
                criteria=args.criteria,
                case_insensitive=args.case_insensitive,
            )
            print(
                f"Wrote {n} rows matching {args.criteria} "
                f"(column {args.col_index}) to {args.out_path}"
            )
            return 0
        except Exception as e:
            print(f"[ERROR] {e}")
            return 1

    # Example *.bam / *.bai (1k reads from 3 samples)
    elif args.cmd == "examples":
        try:
            data_root = res.files("mirpy.data")
            examples = [f for f in data_root.iterdir() if f.name.endswith((".bam", ".bai"))]

            if not examples:
                print("No example BAM files found in mirpy.data.")
                return 0

            if args.out_dir:
                out_path = Path(args.out_dir)
                out_path.mkdir(parents=True, exist_ok=True)

                for f in examples:
                    dest = out_path / f.name
                    with f.open("rb") as src, open(dest, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                print(f"Copied {len(examples)} example files to: {out_path.resolve()}")
            else:
                print("Available example files packaged with mirpy:")
                for f in examples:
                    print(f"  {f.name}")

            return 0
        except Exception as e:
            print(f"[ERROR] {e}")
            return 1

    # Count matched hits
    elif args.cmd == "count":
        return count_matrix(
            bam_paths=args.bams,
            gff_path=args.gff,
            out_path=args.out,
            level=args.level,
            metric=args.metric,
            shift=args.shift,
            max_nh=args.max_nh,
            mode=args.mode,
            log_level=args.log_level,
            log_reads=args.log_reads,
        )
    else:
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
        "view",
        help="Print first N reads from each BAM file (bamnostic-based)."
    )
    t.add_argument(
        "bams",
        nargs="+",
        help="One or more BAM files to test."
    )
    t.add_argument(
        "-n", "--num",
        type=int,
        default=10,
        help="Number of reads per BAM."
    )
    t.add_argument(
        "-r", "--region",
        help="Optional region like 'chr1:1000-2000'."
    )

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

    # Example *.bam and *.bai files
    e = sub.add_parser(
        "examples",
        help="List or export example BAM/BAM index files packaged with miRpy."
    )
    e.add_argument(
        "--out",
        dest="out_dir",
        help="Optional directory to copy example BAM files into. "
             "If omitted, only lists available examples."
    )

    # count (matrix over multiple BAMs)
    c = sub.add_parser(
        "count",
        help="Count mature/precursor miRNA across one or more name-sorted BAMs; produces an output matrix."
    )
    c.add_argument(
        "bams",
        nargs="+",
        help="One or more BAM files or glob patterns (e.g., sample*.namesorted.bam). Must be QNAME-sorted."
    )
    c.add_argument(
        "--gff",
        required=True,
        help="miRBase-like GFF3 (with miRNA + primary_transcript). See also command 'download'."
    )
    c.add_argument(
        "--out",
        required=True,
        help="Output TSV matrix path."
    )
    c.add_argument(
        "--level",
        choices=["mature", "precursor"],
        default="mature",
        help="Row level: 'mature' (default) or 'precursor'."
    )
    c.add_argument(
        "--metric",
        choices=["exact", "approx", "nonspecific", "exact_cpm", "approx_cpm", "nonspecific_cpm"],
        default="exact",
        help="Which metric to report in the matrix (counts or CPM)."
    )
    c.add_argument(
        "--shift",
        type=int,
        default=4,
        help="Max |shift_start|+|shift_end| for approximate matches (default 4)."
    )
    c.add_argument(
        "--max-nh",
        type=int,
        default=50,
        help="Skip reads with NH greater than this (default 50)."
    )
    c.add_argument(
        "--mode",
        choices=["auto", "qname", "nh-bucket"],
        default="auto",
        help="How to read BAMs: 'qname' expects name-sorted; 'nh-bucket' streams unsorted using NH to know when a read is complete; "
             "'auto' tries qname else nh-bucket."
    )
    c.add_argument(
        "--multi",
        choices=["unique", "fractional"], # TODO: mature-fractional
        default="unique",
        help="Strategy on how to handle multiple mapped reads; 'unique' conserves only uniquely mapped reads, 'fractional' "
             "assigns proportional abundance to each mapped location as per 1 / mapped_locations."
    )
    # Debugging assistance
    c.add_argument(
        "--log-level",
        default="INFO",
        choices=["ERROR", "WARNING", "INFO", "DEBUG"],
        help="Logging verbosity (default: INFO)."
    )
    c.add_argument(
        "--log-reads",
        type=int,
        default=0,
        help="When DEBUG, log details for the first N reads per BAM (default: 0)."
    )
    return p

if __name__ == "__main__":
    raise SystemExit(main())
