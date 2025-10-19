import argparse


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mirpy",
        description="Count mature miRNA hits from SAM/BAM against miRBase-like GFF annotations.",
    )
    # Required arguments
    p.add_argument("--gff", required=True, help="Pre-miRNA coordinates (miRBase GFF3).")
    p.add_argument("--out", required=True, help="Output TSV path.")
    # Optional arguments; either --sam or --bam is required
    p.add_argument("--sam", required=True, help="Input samples in SAM format.")
    p.add_argument("--bam", required=True, help="Input samples in BAM format (binary version of SAM).")

    args = p.parse_args(argv)

    return 0
