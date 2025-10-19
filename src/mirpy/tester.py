from __future__ import annotations
import os
import bamnostic as bn


def _get_read_name(aln) -> str:
    # Different libs/files expose different attributes
    for attr in ("query_name", "qname", "read_name"):
        v = getattr(aln, attr, None)
        if v:
            return v
    return "NA"


def miRpyTest(bams: list[str], n: int = 10, region: str | None = None) -> int:
    """
    Print the first N mapped reads from each BAM file using bamnostic.

    If --region is given, a BAM index (.bai) must be present; otherwise we stream
    through the file sequentially (no index required).
    Output is TSV: read_name, locus, MAPQ[, NH if present]
    """
    for bam in bams:
        try:
            bf = bn.AlignmentFile(bam, "rb")
        except Exception as e:
            print(f"[ERROR] Could not open {bam}: {e}")
            return 1

        print(f"== {bam} ==")

        try:
            if region:
                bai_candidates = [bam + ".bai", os.path.splitext(bam)[0] + ".bai"]
                has_index = any(os.path.exists(p) for p in bai_candidates)
                if not has_index:
                    print(f"[ERROR] Region queries require an index (.bai). Not found next to {bam}.")
                    return 2
                try:
                    it = bf.fetch(region=region)
                except KeyError as ke:
                    print(f"[ERROR] Could not resolve region '{region}' in header ({ke}). "
                          f"Check contig names via bf.references.")
                    return 2
                except Exception as e:
                    print(f"[ERROR] Region fetch failed for '{region}': {e}")
                    return 2
            else:
                it = iter(bf)

            printed = 0
            for aln in it:
                if getattr(aln, "is_unmapped", False):
                    continue

                name = _get_read_name(aln)
                rname = getattr(aln, "reference_name", "")
                start_1b = (getattr(aln, "pos", 0) or 0) + 1  # bamnostic uses 0-based pos
                end_1b = getattr(aln, "reference_end", start_1b)
                strand = "-" if getattr(aln, "is_reverse", False) else "+"
                mapq = getattr(aln, "mapq", "")

                # NH tag (optional)
                nh_s = ""
                try:
                    nh = aln.opt("NH")
                    nh_s = f"\tNH={nh}"
                except Exception:
                    pass

                # TSV line (read_name  locus  MAPQ  [NH])
                print(f"{name}\t{rname}:{start_1b}-{end_1b}({strand})\tMAPQ={mapq}{nh_s}")

                printed += 1
                if printed >= n:
                    break

            if printed == 0:
                if region:
                    print("[info] No mapped reads found in region (or region outside data).")
                else:
                    print("[info] No mapped reads found.")
        finally:
            bf.close()

    return 0
