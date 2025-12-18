from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import gzip
import os
import glob
import logging
import bamnostic as bn
import traceback
import psutil
import sys
import gc
from .miRpyClasses import AlignmentData, Mature


def _extract_alignment_data(aln) -> Optional[AlignmentData]:
    """Extract only necessary data from bamnostic alignment to reduce memory."""
    if getattr(aln, "is_unmapped", False):
        return None

    chr_ = getattr(aln, "reference_name", None)
    if chr_ is None:
        return None

    # bamnostic uses 'pos' (0-based) instead of 'reference_start'
    start_1b = (getattr(aln, "pos", 0) or 0) + 1
    end_1b = getattr(aln, "reference_end", start_1b)
    is_reverse = getattr(aln, "is_reverse", False)

    try:
        nm = aln.opt("NM")
    except (KeyError, AttributeError):
        nm = 1  # Assume non-perfect match if NM tag missing

    return AlignmentData(chr_, start_1b, end_1b, is_reverse, nm)


def _get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024 # Current memory usage in MB


def _make_logger(level: str) -> logging.Logger:
    lvl = getattr(logging, level.upper(), logging.INFO)
    logger = logging.getLogger("mirpy.map")
    if not logger.handlers:
        handler = logging.StreamHandler()
        fmt = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    logger.setLevel(lvl)
    return logger


def _open_text_auto(path: str | Path):
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "rt", encoding="utf-8", errors="replace")


def _parse_attrs(attr_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for kv in attr_field.strip().split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def _load_mature_only_gff(
    gff_path: str | Path,
    logger: logging.Logger | None = None,
) -> Dict[str, Dict[str, List[Mature]]]:
    """Load mature miRNA tab separated TSV database from e.g. miRBase"""
    # Index of all mature miRNAs as per Mature-class
    mature_index: Dict[str, Dict[str, List[Mature]]] = {}

    with _open_text_auto(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start_s, end_s, _score, strand, _phase, attrs = cols
            #if feature != "miRNA":
            #    # assume prefiltered but skip any stray non-mature features
            #    continue

            A = _parse_attrs(attrs)
            name = A.get("Name")
            if not name:
                continue
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue

            m = Mature(
                name=name,
                chr=chrom,
                strand=strand,
                start=start,
                end=end,
            )
            mature_index.setdefault(chrom, {}).setdefault(strand, []).append(m)

    # Make sure the mature miRNA GFFs are sorted for coordinates,
    # as linear search relies on them being in ascending order
    for chr_, strands in mature_index.items():
        for strand, matures in strands.items():
            matures.sort(key=lambda x: x.start)

    # Info/debugging of the loaded miRNA GFF
    if logger:
        n_matures = sum(len(v2) for v in mature_index.values() for v2 in v.values())
        logger.info(f"GFF (mature-only) loaded: {n_matures} mature miRNAs")
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("GFF overview (mature-only):")
            for chr_ in sorted(mature_index.keys()):
                for st, matures in mature_index[chr_].items():
                    logger.debug(f"  chr={chr_!r} strand={st}: {len(matures)} matures")
            # List all candidate mature miRNA in debug
            all_matures = []
            for chr_ in sorted(mature_index.keys()):
                for st, matures in mature_index[chr_].items():
                    all_matures.extend(matures)
            logger.debug("Head and tail of candidate mature miRNAs from GFF:")
            for m in all_matures[:5] + all_matures[-5:]:
                logger.debug(
                    f"  {m.name}\t{m.chr}:{m.start}-{m.end}({m.strand})"
                )

    return mature_index


def _get_read_name(aln) -> str:
    for attr in ("query_name", "qname", "read_name"):
        v = getattr(aln, attr, None)
        if v:
            return v
    return ""


def _expand_bam_patterns(bams: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for pat in bams:
        matches = glob.glob(pat) if any(ch in pat for ch in "*?[]") else (
            [pat] if os.path.exists(pat) else []
        )
        for m in sorted(matches):
            if m not in seen:
                seen.add(m)
                out.append(m)
        if not matches:
            print(f"[WARNING] No BAMs matched: {pat}")
    return out


def _map_one_bam(
    bam_path: str | Path,
    mature_index: Dict[str, Dict[str, List[Mature]]],
    *,
    max_shift: int,
    max_nh: int,
    multi: str = "unique",  # "unique" or "fractional"
    dist: str = "strict",   # "strict" or "containment"
    logger: logging.Logger | None = None,
    log_reads: int = 1,
) -> Tuple[Dict[str, List[float]], int]:
    """
    Stream an unsorted BAM, buffering per read until we've seen NH alignments.
    Same mapping logic as _map_qname_sorted, but grouping reads by QNAME without
    requiring name-sorted input.
    """
    if logger:
        logger.info(f"Mapping (NH-bucket) in {bam_path}")
        if logger.isEnabledFor(logging.DEBUG):
            try:
                with bn.AlignmentFile(str(bam_path), "rb") as btmp:
                    refs = list(getattr(btmp, "references", []))
                    logger.debug(f"BAM references: {refs}")
                    gff_contigs = set(mature_index.keys())
                    bam_contigs = set(refs)
                    only_gff = sorted(gff_contigs - bam_contigs)
                    only_bam = sorted(bam_contigs - gff_contigs)
                    logger.debug(f"Contigs in GFF not in BAM (first 20): {only_gff[:20]}")
                    logger.debug(f"Contigs in BAM not in GFF (first 20): {only_bam[:20]}")
            except Exception:
                pass

    mature_counts: Dict[str, List[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    aligned_reads = 0
    matched_reads = 0

    # debugging counters
    reason_nh_too_big = 0
    reason_no_chr = 0
    reason_no_strand = 0
    reason_no_mature_overlap = 0
    reason_distance_too_big = 0
    reason_multi_conflict = 0

    # Memory efficiency update using AlignmentData-class
    buckets: Dict[str, List[AlignmentData]] = {}
    nh_by_read: Dict[str, Optional[int]] = {}
    reads_logged = 0
    progress_every = 1000000

    # Define internal function for processing buckets
    def process_bucket(qname: str):
        nonlocal reads_logged, matched_reads
        nonlocal reason_nh_too_big, reason_no_chr, reason_no_strand, reason_no_mature_overlap, \
            reason_distance_too_big, reason_multi_conflict

        bucket = buckets.pop(qname, [])
        nhpop = nh_by_read.pop(qname, None)
        genomic_matches = nhpop if isinstance(nhpop, int) else len(bucket)
        if genomic_matches > max_nh:
            reason_nh_too_big += 1
            return

        # Sets are more memory efficient with minimum alignment info used
        exact_set = set()
        approx_set = set()
        hit_any_matures = False

        for aln_data in bucket:
            if aln_data.chr not in mature_index:
                reason_no_chr += 1
                continue

            strand = "-" if aln_data.is_reverse else "+"
            if strand not in mature_index[aln_data.chr]:
                reason_no_strand += 1
                continue

            non_pm = 0 if (aln_data.nm == 0) else 1

            # Mature-objects based on chr and strand ('+' or '-')
            mature_list = mature_index[aln_data.chr][strand]

            # Loop over annotated mature miRNAs
            for m in mature_list:
                # Linear search optimization;
                # if mature.start > aln_data.end,
                # all subsequent matures will also start after our read ends so we can quit
                if m.start > aln_data.end_1b:
                    break

                # Linear search optimization;
                # If this specific miRNA ends before our read starts it can be skipped
                # (might be not very efficient if-checking, to be tested)
                if m.end < aln_data.start_1b:
                    continue

                # Distance calculation coded here to avoid function call overheads
                if dist == "strict":
                    d = abs(m.start - aln_data.start_1b) + abs(m.end - aln_data.end_1b)
                else: # Should be 'containment' due to sanity checks
                    d = max(0, m.start - aln_data.start_1b) + max(0, aln_data.end_1b - m.end)

                if d <= max_shift:
                    if non_pm == 0 and d == 0:
                        exact_set.add(m.name)
                    approx_set.add(m.name)
                    hit_any_matures = True
                else:
                    reason_distance_too_big += 1

        if not hit_any_matures:
            reason_no_mature_overlap += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                a0 = bucket[0] if bucket else None
                if a0 is not None:
                    st = "-" if a0.is_reverse else "+"
                    logger.debug(
                        f"{qname}: no mature overlap. "
                        f"chr={a0.chr!r}, strand={st}, start/end=({a0.start_1b}, {a0.end_1b})"
                    )
                reads_logged += 1
            return

        matched = list(approx_set) # All matches within shift tolerance
        K = len(matched)
        if K == 0:
            reason_no_mature_overlap += 1
            return

        if multi == "unique" and K > 1:
            reason_multi_conflict += 1
            return

        if multi == "fractional" and K > 1:
            w = 1.0 / float(K)
        else:
            w = 1.0

        for name in matched:
            counts = mature_counts[name]
            if name in exact_set:
                counts[0] += w  # exact
                counts[1] += w  # approx (exact is a subset of approx)
            else:
                counts[1] += w  # approx only
            counts[2] += w      # nonspecific always gets weight

        # Count the read as being matched as it passed through all criteria and wasn't returned
        matched_reads += 1


    # Main loop reading the aligned reads in the BAM file
    try:
        with bn.AlignmentFile(str(bam_path), "rb") as bam:
            for aln in bam:
                if getattr(aln, "is_unmapped", False):
                    continue
                aligned_reads += 1

                if logger and aligned_reads % progress_every == 0:
                    current_mem = _get_memory_usage()
                    logger.info(f"Processed {aligned_reads:,} alignments... (Memory: {current_mem:.1f} MB, buffer: {len(buckets)} reads)")

                qn = _get_read_name(aln)

                # Extract minimal data (more memory efficient that prior wasteful approach)
                aln_data = _extract_alignment_data(aln)
                if aln_data is None:  # Skip if extraction failed
                    continue

                if qn not in buckets:
                    buckets[qn] = []
                    try:
                        nh_by_read[qn] = aln.opt("NH")
                    except (KeyError, AttributeError):
                        nh_by_read[qn] = None

                buckets[qn].append(aln_data)

                nh = nh_by_read.get(qn)
                if isinstance(nh, int) and len(buckets[qn]) >= nh:
                    process_bucket(qn)
        for qn in list(buckets.keys()):
            process_bucket(qn)

    except Exception as e:
        if logger:
            logger.error(f"Error in NH-bucket processing loop at alignment {aligned_reads:,}")
            logger.error(f"Buffer size at error: {len(buckets):,} reads")
            logger.error(f"Memory at error: {_get_memory_usage():.1f} MB")
            logger.error(f"Error: {e}")
            logger.error(traceback.format_exc())
        raise

    if logger:
        logger.info(
            f"Done {bam_path}: aligned={aligned_reads}, features_matched={len(mature_counts)}, "
            f"total_matched_reads={matched_reads}"
        )
        logger.debug(
            f"Reasons for non-matching / dropped reads in {bam_path}: "
            f"nh_too_big={reason_nh_too_big}, "
            f"no_chr={reason_no_chr}, no_strand={reason_no_strand}, "
            f"no_mature_overlap={reason_no_mature_overlap}, "
            f"distance_too_big={reason_distance_too_big}, "
            f"multi_mature_conflict={reason_multi_conflict}"
        )

    return dict(mature_counts), aligned_reads


def map_matrix(
    bam_paths: List[str],
    gff_path: str | Path,
    out_path: str | Path,
    *,
    metric: str = "exact",       # exact | approx | nonspecific | exact_cpm | approx_cpm | nonspecific_cpm
    shift: int = 4,
    max_nh: int = 50,
    multi: str = "fractional",   # "unique" or "fractional"
    dist: str = "strict",        # "strict" or "containment"
    log_level: str = "INFO",
    log_reads: int = 0,
) -> int:
    """
    Main function for building a count-like matrix with rows = mature miRNAs, columns = samples.
    """
    logger = _make_logger(log_level)

    # System info for debugging
    logger.debug(f"Python version: {sys.version}")
    logger.debug(f"Starting memory: {_get_memory_usage():.1f} MB")
    logger.debug(f"Available memory: {psutil.virtual_memory().available / 1024 / 1024:.1f} MB")
    logger.debug(f"Total memory: {psutil.virtual_memory().total / 1024 / 1024:.1f} MB")

    # Load mature annotations and the optimization read locations for faster processing
    mature_index = _load_mature_only_gff(gff_path, logger=logger)

    # Expand BAM patterns
    bam_list = _expand_bam_patterns(bam_paths)
    if not bam_list:
        logger.error("No BAMs found.")
        return 1

    logger.info(
        f"{len(bam_list)} BAM(s) to process; shift={shift}, "
        f"max_nh={max_nh}, multi={multi}"
    )
    sample_names = [Path(b).stem for b in bam_list]

    per_sample_counts: List[Dict[str, List[float]]] = []
    per_sample_aligned: List[int] = []

    # Prepare output paths
    outp = Path(out_path)
    outp.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = outp.with_suffix('.tmp.tsv')

    # Metric selection
    metric = metric.lower()
    cpm = metric.endswith("_cpm")
    base = metric.replace("_cpm", "")
    metricid = {"exact": 0, "approx": 1, "nonspecific": 2}.get(base)
    if metricid is None:
        logger.error(
            "--metric must be one of: exact, approx, nonspecific, "
            "exact_cpm, approx_cpm, nonspecific_cpm"
        )
        return 2
    # Shift distance selection
    if dist not in ["strict", "containment"]:
        logger.error(
            "--dist must be one of: strict, containment"
        )
        return 3

    for bamf, b in enumerate(bam_list, 1):
        logger.info(f"Processing BAM {bamf}/{len(bam_list)}: {b}")
        try:
            mc, aligned = _map_one_bam(
                b,
                mature_index,
                max_shift=shift,
                max_nh=max_nh,
                multi=multi,
                dist=dist,
                logger=logger,
                log_reads=log_reads,
            )
            per_sample_counts.append(mc)
            per_sample_aligned.append(aligned)

            # Clean up with garbage collection after processing each BAM
            gc.collect()

            # Write temporary matrix after each BAM
            logger.info(f"Writing temporary matrix to {tmp_path}...")
            features_so_far = sorted({m for mc in per_sample_counts for m in mc.keys()})
            samples_so_far = sample_names[:bamf]

            with open(tmp_path, "w", encoding="utf-8") as fh:
                fh.write("feature\t" + "\t".join(samples_so_far) + "\n")
                for feat in features_so_far:
                    vals: List[str] = []
                    for sampid, mc in enumerate(per_sample_counts):
                        cnt = mc.get(feat, [0.0, 0.0, 0.0])[metricid]
                        if cpm:
                            aligned_val = per_sample_aligned[sampid] or 1
                            val = 1_000_000 * cnt / aligned_val
                            vals.append(f"{val:.4f}")
                        else:
                            vals.append(f"{cnt:.4f}")
                    fh.write(feat + "\t" + "\t".join(vals) + "\n")

            logger.info(f"Temporary matrix saved ({len(features_so_far)} features with {len(samples_so_far)} samples)")
            logger.info(f"Successfully processed {bamf}/{len(bam_list)} BAMs")
            logger.info(f"Current memory: {_get_memory_usage():.1f} MB")

            if logger.isEnabledFor(logging.DEBUG):
                nonzero_exact = sum(1 for v in mc.values() if v[0] > 0)
                nonzero_approx = sum(1 for v in mc.values() if v[1] > 0)
                nonzero_ns = sum(1 for v in mc.values() if v[2] > 0)
                logger.debug(
                    f"{b}: aligned={aligned}, features={len(mc)}, "
                    f"nonzero_exact={nonzero_exact}, "
                    f"nonzero_approx={nonzero_approx}, "
                    f"nonzero_nonspecific={nonzero_ns}, "
                    f"multi={multi}"
                )
        except Exception as e:
            logger.error(f"{b}: {e}")
            logger.debug("Traceback after Exception on _map_one_bam:\n" + traceback.format_exc())
            return 1

    # Feature set is all observed mature miRNA names
    features = sorted({m for mc in per_sample_counts for m in mc.keys()})
    if not features:
        logger.warning(
            "No mature miRNAs matched any BAM. Check genome build (chr1 vs 1), "
            "strand orientation, or shift/NH/multi settings."
        )

    # Write matrix
    outp = Path(out_path)
    outp.parent.mkdir(parents=True, exist_ok=True)

    with open(outp, "w", encoding="utf-8") as fh:
        fh.write("feature\t" + "\t".join(sample_names) + "\n")
        for feat in features:
            vals: List[str] = []
            for sampid, mc in enumerate(per_sample_counts):
                cnt = mc.get(feat, [0.0, 0.0, 0.0])[metricid]
                if cpm:
                    aligned = per_sample_aligned[sampid] or 1
                    val = 1_000_000 * cnt / aligned
                    vals.append(f"{val:.4f}")
                else:
                    vals.append(f"{cnt:.4f}")
            fh.write(feat + "\t" + "\t".join(vals) + "\n")

    logger.info(f"Wrote matrix to {outp} with {len(features)} features and {len(sample_names)} samples")

    try:
        tmp_path.unlink()
        logger.info(f"Removed temporary file: {tmp_path}")
    except Exception as e:
        logger.warning(f"Could not remove temporary file {tmp_path}: {e}")

    logger.info(f"Final memory usage: {_get_memory_usage():.1f} MB")
    return 0

