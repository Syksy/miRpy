import pytest
from mirpy import cli

def test_cli_argument_parsing(monkeypatch):
    # Example fake command line input (as if typed into terminal)
    argv = [
        "--gff", "data/mirbase.gff3",
        "--out", "results/output.tsv",
        "--sam", "alignments/sample.sam",
        "--bam", "alignments/sample.bam",
    ]

    # Call main() directly, so that the argparse will parse this list
    result = cli.main(argv)

    # The function should return 0 (success)
    assert result == 0
    