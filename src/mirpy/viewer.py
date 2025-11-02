from __future__ import annotations
import os
import bamnostic as bn
from pathlib import Path
import glob


def _get_read_name(aln) -> str:
    # Different libs/files expose different attributes
    for attr in ("query_name", "qname", "read_name"):
        v = getattr(aln, attr, None)
        if v:
            return v
    return "NA"

def _get_read_full(aln) -> str:
    # Iterate across all available fields
    out = {}
    for attr in dir(aln):
        if attr.startswith("_"):
            continue  # skip internals
        try:
            val = getattr(aln, attr)
        except Exception:
            continue
        # Keep only primitive-like values (avoid methods)
        if callable(val):
            continue
        out[attr] = val
    # bamnostic stores optional tags like NM:i, AS:i, etc in aln.tags (a dict or list of tuples)
    tags = {}
    try:
        # Newer bamnostic exposes aln.tags as dict-like; older as list of tuples
        t = getattr(aln, "tags", None)
        if isinstance(t, dict):
            tags = t
        elif isinstance(t, list):
            tags = dict(t)
        else:
            # try the opt() interface
            if hasattr(aln, "opt"):
                for tag in ("AS", "NM", "NH", "MD", "XS", "YT", "XN", "XO", "XG"):
                    try:
                        tags[tag] = aln.opt(tag)
                    except Exception:
                        pass
    except Exception:
        pass

    if tags:
        out["tags"] = tags

    return out

def _print_read_full(aln, delimiter: str = " | ") -> str:
    outs = _get_read_full(aln)
    parts = []
    for k, v in sorted(outs.items()):
        if k == "tags" and isinstance(v, dict):
            tags_str = ", ".join(f"{tk}:{tv}" for tk, tv in sorted(v.items()))
            parts.append(f"tags=[{tags_str}]")
        else:
            parts.append(f"{k}={v}")
    return  delimiter.join(parts)

def _expand_bam_patterns(bams: list[str]) -> list[str]:
    """Expand glob patterns (sample*.bam) cross-platform; keep order; de-dupe."""
    seen = set()
    out: list[str] = []
    for pat in bams:
        # If pattern contains wildcards, expand; otherwise treat as literal path
        matches = glob.glob(pat) if any(ch in pat for ch in "*?[]") else ([pat] if os.path.exists(pat) else [])
        if not matches:
            print(f"[WARNING] No BAMs matched: {pat}")
        for m in sorted(matches):
            if m not in seen:
                seen.add(m)
                out.append(m)
    return out


def view_bam_head(bams: list[str], n: int = 10, region: str | None = None) -> int:
    """
    Print the first N mapped reads from each BAM file using bamnostic.
    Supports wildcards like *.bam on Windows.

    If --region is given, a BAM index (.bai) must be present; otherwise we stream
    through the file sequentially (no index required).
    Output is TSV: read_name, locus, MAPQ[, NH if present]
    """
    bam_paths = _expand_bam_patterns(bams)
    if not bam_paths:
        print("[ERROR] No BAM files found.")
        return 1

    for bam in bam_paths:
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
