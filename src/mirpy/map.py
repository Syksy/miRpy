from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import gzip
import os
import glob
import logging
import bamnostic as bn
import traceback
from miRpyClasses import AlignmentData, Mature

def _make_logger(level: str) -> logging.Logger:
    lvl = getattr(logging, level.upper(), logging.INFO)
    logger = logging.getLogger("mirpy.map")
    if not logger.handlers:
        handler = logging.StreamHandler()
        fmt = logging.Formatter("[%(levelname)s] %(message)s")
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


def load_mature_only_gff(
    gff_path: str | Path,
    logger: logging.Logger | None = None,
) -> Dict[str, Dict[str, List[Mature]]]:
    """
    Load a *mature-only* miRBase-like GFF3.

    Assumes the GFF has been prefiltered to feature == "miRNA".
    Returns an index: mature_index[chr][strand] -> List[Mature].
    """
    mature_index: Dict[str, Dict[str, List[Mature]]] = {}

    with _open_text_auto(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start_s, end_s, _score, strand, _phase, attrs = cols
            if feature != "miRNA":
                # assume prefiltered but skip any stray non-mature features
                continue

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

    # Sort by genomic start for each chr/strand
    for chr_, strands in mature_index.items():
        for st, matures in strands.items():
            matures.sort(key=lambda x: x.start)

    if logger:
        n_matures = sum(len(v2) for v in mature_index.values() for v2 in v.values())
        logger.info(f"GFF (mature-only) loaded: {n_matures} mature miRNAs")
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("GFF overview (mature-only):")
            for chr_ in sorted(mature_index.keys()):
                for st, matures in mature_index[chr_].items():
                    logger.debug(f"  chr={chr_!r} strand={st}: {len(matures)} matures")
            # List all candidate mature miRNA in debug
            logger.debug("Full list of candidate mature miRNAs from GFF:")
            for chr_ in sorted(mature_index.keys()):
                for st, matures in mature_index[chr_].items():
                    for m in matures:
                        logger.debug(
                            f"  {m.name}\t{m.chr}:{m.start}-{m.end}({m.strand})"
                        )

    return mature_index


def _is_qname_sorted(bam: bn.AlignmentFile, sample: int = 2000) -> bool:
    """Heuristic: ensure non-decreasing query_name for first N records."""
    it = iter(bam)
    last: Optional[str] = None
    checked = 0
    for aln in it:
        qn = getattr(aln, "qname", None) or getattr(aln, "query_name", None) or ""
        if last is not None and qn < last:
            return False
        last = qn
        checked += 1
        if checked >= sample:
            break
    return True


def _aln_locus_1b(aln) -> Tuple[int, int]:
    start_1b = (getattr(aln, "pos", 0) or 0) + 1
    end_1b = getattr(aln, "reference_end", start_1b)
    return start_1b, end_1b


# Debugging; this compares read against absolute known miRNA site start/end boundaries
def _distance_to_mature(start_1b: int, end_1b: int, m: Mature) -> int:
    return abs(m.start - start_1b) + abs(m.end - end_1b)


# Debugging; this compares against how much bases the aligned read crosses over the known miRNA boundaries
def _distance_out_mature(start_1b: int, end_1b: int, m: Mature) -> int:
    return max(0, m.start - start_1b) + max(0, end_1b - m.end)


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
    shift: int,
    max_nh: int,
    mode: str = "qname",    # "qname", "nh-bucket", or "auto"
    multi: str = "unique",  # "unique" or "fractional"
    logger: logging.Logger | None = None,
    log_reads: int = 0,
) -> Tuple[Dict[str, List[float]], int]:
    """
    Map per mature for a single BAM (mature-only model).

    Returns (mature_counts, aligned_reads).

    mature_counts[name] -> [exact, approx, nonspecific] (floats, to support fractional mapping)
    aligned_reads        -> number of mapped alignments observed
    """
    try:
        test = bn.AlignmentFile(str(bam_path), "rb")
    except Exception as e:
        raise RuntimeError(f"Could not open BAM: {bam_path}: {e}")

    try:
        if mode == "qname":
            if not _is_qname_sorted(test):
                raise RuntimeError(
                    "BAM is not name-sorted (QNAME). "
                    "Use --mode nh-bucket or sort with: samtools sort -n -o namesorted.bam input.bam"
                )
            return _map_qname_sorted(
                bam_path,
                mature_index,
                shift=shift,
                max_nh=max_nh,
                multi=multi,
                logger=logger,
                log_reads=log_reads,
            )
        elif mode == "nh-bucket":
            return _map_unsorted_nh_bucket(
                bam_path,
                mature_index,
                shift=shift,
                max_nh=max_nh,
                multi=multi,
                logger=logger,
                log_reads=log_reads,
            )
        elif mode == "auto":
            if _is_qname_sorted(test):
                return _map_qname_sorted(
                    bam_path,
                    mature_index,
                    shift=shift,
                    max_nh=max_nh,
                    multi=multi,
                    logger=logger,
                    log_reads=log_reads,
                )
            return _map_unsorted_nh_bucket(
                bam_path,
                mature_index,
                shift=shift,
                max_nh=max_nh,
                multi=multi,
                logger=logger,
                log_reads=log_reads,
            )
        else:
            raise RuntimeError(f"Unknown mode: {mode}")
    finally:
        try:
            test.close()
        except Exception:
            pass


def _map_qname_sorted(
    bam_path: str | Path,
    mature_index: Dict[str, Dict[str, List[Mature]]],
    *,
    shift: int,
    max_nh: int,
    multi: str = "unique",
    logger: logging.Logger | None = None,
    log_reads: int = 0,
) -> Tuple[Dict[str, List[float]], int]:
    """
    Stream a QNAME-sorted BAM and map reads at the mature level.

    multi:
      - "unique"     : keep only reads matching exactly one mature miRNA.
      - "fractional" : split read weight evenly across all matching matures.
    """
    if logger:
        logger.info(f"Mapping (QNAME) in {bam_path}")
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

    mature_counts: Dict[str, List[float]] = {}
    aligned_reads = 0

    # debugging counters
    reason_nh_too_big = 0
    reason_no_chr = 0
    reason_no_strand = 0
    reason_no_mature_overlap = 0
    reason_distance_too_big = 0
    reason_multi_conflict = 0

    def ensure(name: str):
        if name not in mature_counts:
            mature_counts[name] = [0.0, 0.0, 0.0]

    current_qname: Optional[str] = None
    bucket: List = []
    reads_logged = 0
    progress_every = 100000

    def process_bucket():
        nonlocal reads_logged
        nonlocal reason_nh_too_big, reason_no_chr, reason_no_strand
        nonlocal reason_no_mature_overlap, reason_distance_too_big, reason_multi_conflict

        if not bucket:
            return

        # Determine NH
        nh = None
        for aln in bucket:
            try:
                nh = aln.opt("NH")
                break
            except Exception:
                pass
        genomic_matches = nh if isinstance(nh, int) else len(bucket)
        if genomic_matches > max_nh:
            reason_nh_too_big += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                logger.debug(f"{_get_read_name(bucket[0])}: skip NH={genomic_matches} > max_nh={max_nh}")
                reads_logged += 1
            return

        exact_tmp: Dict[str, None] = {}
        approx_tmp: Dict[str, None] = {}
        matched_names: Dict[str, None] = {}
        sum_matches = 0
        hit_any_matures = False

        for aln in bucket:
            if getattr(aln, "is_unmapped", False):
                continue
            chr_ = getattr(aln, "reference_name", None)
            if chr_ is None:
                continue
            if chr_ not in mature_index:
                reason_no_chr += 1
                continue
            strand = "-" if getattr(aln, "is_reverse", False) else "+"
            if strand not in mature_index[chr_]:
                reason_no_strand += 1
                continue

            start_1b, end_1b = _aln_locus_1b(aln)
            nm = None
            try:
                nm = aln.opt("NM")
            except Exception:
                pass
            non_pm = 0 if (nm == 0) else 1

            cand_matures = mature_index[chr_][strand]
            local_hit = False
            for m in cand_matures:
                # overlap check
                if end_1b < m.start or start_1b > m.end:
                    continue
                #d = _distance_to_mature(start_1b, end_1b, m)
                d = _distance_out_mature(start_1b, end_1b, m)
                if d <= shift:
                    # Perfect hit
                    if non_pm == 0 and d == 0:
                        exact_tmp[m.name] = None
                    # Otherwise approximate hit, all other logic appended on top of perfect hit
                    approx_tmp[m.name] = None
                    matched_names[m.name] = None
                    sum_matches += 1
                    local_hit = True
                else:
                    reason_distance_too_big += 1
            if local_hit:
                hit_any_matures = True

        if not hit_any_matures:
            reason_no_mature_overlap += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                a0 = next((a for a in bucket if not getattr(a, "is_unmapped", False)), None)
                if a0 is not None:
                    chr_ = getattr(a0, "reference_name", None)
                    st = "-" if getattr(a0, "is_reverse", False) else "+"
                    logger.debug(
                        f"{_get_read_name(a0)}: no mature overlap at chr={chr_!r}, strand={st}, "
                        f"start/end={_aln_locus_1b(a0)}"
                    )
                reads_logged += 1
            return

        matched = list(matched_names.keys())
        K = len(matched)
        if K == 0:
            reason_no_mature_overlap += 1
            return

        if multi == "unique" and K > 1:
            reason_multi_conflict += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                logger.debug(
                    f"{_get_read_name(bucket[0])}: multi-mature conflict {matched}, "
                    f"sum_matches={sum_matches}, NH={genomic_matches}"
                )
                reads_logged += 1
            return

        # weight per mature
        if multi == "fractional" and K > 1:
            w = 1.0 / float(K)   # split across matching matures
        else:
            w = 1.0

        for name in matched:
            ensure(name)
            if name in exact_tmp:
                mature_counts[name][0] += w  # exact
                mature_counts[name][1] += w  # approx (exact is a subset of approx)
            elif name in approx_tmp:
                mature_counts[name][1] += w  # approx only
            mature_counts[name][2] += w      # nonspecific always gets weight

    with bn.AlignmentFile(str(bam_path), "rb") as bam:
        for aln in bam:
            if getattr(aln, "is_unmapped", False):
                continue
            aligned_reads += 1

            if logger and aligned_reads % progress_every == 0:
                logger.info(f"Processed {aligned_reads:,} alignments...")

            qn = _get_read_name(aln)
            if current_qname is None:
                current_qname = qn
            if qn != current_qname:
                process_bucket()
                bucket = [aln]
                current_qname = qn
            else:
                bucket.append(aln)

    process_bucket()

    if logger:
        logger.info(
            f"Done {bam_path}: aligned={aligned_reads}, features_matched={len(mature_counts)}"
        )
        logger.debug(
            f"Reasons (QNAME) for non-matching / dropped reads in {bam_path}: "
            f"nh_too_big={reason_nh_too_big}, "
            f"no_chr={reason_no_chr}, no_strand={reason_no_strand}, "
            f"no_mature_overlap={reason_no_mature_overlap}, "
            f"distance_too_big={reason_distance_too_big}, "
            f"multi_mature_conflict={reason_multi_conflict}"
        )

    return mature_counts, aligned_reads


def _map_unsorted_nh_bucket(
    bam_path: str | Path,
    mature_index: Dict[str, Dict[str, List[Mature]]],
    *,
    shift: int,
    max_nh: int,
    multi: str = "unique",
    logger: logging.Logger | None = None,
    log_reads: int = 0,
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

    mature_counts: Dict[str, List[float]] = {}
    aligned_reads = 0

    # debugging counters
    reason_nh_too_big = 0
    reason_no_chr = 0
    reason_no_strand = 0
    reason_no_mature_overlap = 0
    reason_distance_too_big = 0
    reason_multi_conflict = 0

    def ensure(name: str):
        if name not in mature_counts:
            mature_counts[name] = [0.0, 0.0, 0.0]

    buckets: Dict[str, List] = defaultdict(list)
    nh_by_read: Dict[str, Optional[int]] = {}
    reads_logged = 0
    progress_every = 100000

    def process_bucket(qn: str):
        nonlocal reads_logged
        nonlocal reason_nh_too_big, reason_no_chr, reason_no_strand, reason_no_mature_overlap, \
            reason_distance_too_big, reason_multi_conflict

        bucket = buckets.pop(qn, [])
        nh = nh_by_read.pop(qn, None)
        genomic_matches = nh if isinstance(nh, int) else len(bucket)
        if genomic_matches > max_nh:
            reason_nh_too_big += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                logger.debug(f"{qn}: skip NH={genomic_matches} > max_nh={max_nh}")
                reads_logged += 1
            return

        exact_tmp: Dict[str, None] = {}
        approx_tmp: Dict[str, None] = {}
        matched_names: Dict[str, None] = {}
        sum_matches = 0
        hit_any_matures = False

        for aln in bucket:
            if getattr(aln, "is_unmapped", False):
                continue
            chr_ = getattr(aln, "reference_name", None)
            if chr_ is None:
                continue
            if chr_ not in mature_index:
                reason_no_chr += 1
                continue
            strand = "-" if getattr(aln, "is_reverse", False) else "+"
            if strand not in mature_index[chr_]:
                reason_no_strand += 1
                continue

            start_1b, end_1b = _aln_locus_1b(aln)
            nm = None
            try:
                nm = aln.opt("NM")
            except Exception:
                pass
            non_pm = 0 if (nm == 0) else 1

            cand_matures = mature_index[chr_][strand]
            local_hit = False
            for m in cand_matures:
                if end_1b < m.start or start_1b > m.end:
                    continue
                #d = _distance_to_mature(start_1b, end_1b, m)
                d = _distance_out_mature(start_1b, end_1b, m)
                if d <= shift:
                    # Perfect hit
                    if non_pm == 0 and d == 0:
                        exact_tmp[m.name] = None
                    # Otherwise approximate hit, all other logic appended on top of perfect hit
                    approx_tmp[m.name] = None
                    matched_names[m.name] = None
                    sum_matches += 1
                    local_hit = True
                else:
                    reason_distance_too_big += 1
            if local_hit:
                hit_any_matures = True

        if not hit_any_matures:
            reason_no_mature_overlap += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                a0 = next((a for a in bucket if not getattr(a, "is_unmapped", False)), None)
                if a0 is not None:
                    chr_ = getattr(a0, "reference_name", None)
                    st = "-" if getattr(a0, "is_reverse", False) else "+"
                    logger.debug(
                        f"{_get_read_name(a0)}: no mature overlap. "
                        f"chr={chr_!r}, strand={st}, start/end={_aln_locus_1b(a0)}"
                    )
                reads_logged += 1
            return

        matched = list(matched_names.keys())
        K = len(matched)
        if K == 0:
            reason_no_mature_overlap += 1
            return

        if multi == "unique" and K > 1:
            reason_multi_conflict += 1
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                logger.debug(
                    f"{qn}: multi-mature conflict {matched}, "
                    f"sum_matches={sum_matches}, NH={genomic_matches}"
                )
                reads_logged += 1
            return

        if multi == "fractional" and K > 1:
            w = 1.0 / float(K)
        else:
            w = 1.0

        for name in matched:
            ensure(name)
            if name in exact_tmp:
                mature_counts[name][0] += w
                mature_counts[name][1] += w
            elif name in approx_tmp:
                mature_counts[name][1] += w
            mature_counts[name][2] += w

    for aln in bn.AlignmentFile(str(bam_path), "rb"):
        if getattr(aln, "is_unmapped", False):
            continue
        aligned_reads += 1
        if logger and aligned_reads % progress_every == 0:
            logger.info(f"Processed {aligned_reads:,} alignments...")

        qn = _get_read_name(aln)
        buckets[qn].append(aln)

        if qn not in nh_by_read or nh_by_read[qn] is None:
            try:
                nh_by_read[qn] = aln.opt("NH")
            except Exception:
                nh_by_read.setdefault(qn, None)

        nh = nh_by_read.get(qn)
        if isinstance(nh, int) and len(buckets[qn]) >= nh:
            process_bucket(qn)

    for qn in list(buckets.keys()):
        process_bucket(qn)

    if logger:
        logger.info(
            f"Done {bam_path}: aligned={aligned_reads}, features_matched={len(mature_counts)}"
        )
        logger.debug(
            f"Reasons (NH-bucket) for non-matching / dropped reads in {bam_path}: "
            f"nh_too_big={reason_nh_too_big}, "
            f"no_chr={reason_no_chr}, no_strand={reason_no_strand}, "
            f"no_mature_overlap={reason_no_mature_overlap}, "
            f"distance_too_big={reason_distance_too_big}, "
            f"multi_mature_conflict={reason_multi_conflict}"
        )

    return mature_counts, aligned_reads


def map_matrix(
    bam_paths: List[str],
    gff_path: str | Path,
    out_path: str | Path,
    *,
    metric: str = "exact",   # exact | approx | nonspecific | exact_cpm | approx_cpm | nonspecific_cpm
    shift: int = 4,
    max_nh: int = 50,
    mode: str = "auto",
    multi: str = "fractional",   # "unique" or "fractional"
    log_level: str = "INFO",
    log_reads: int = 0,
) -> int:
    """
    Build a matrix with rows = mature miRNAs, columns = samples.

    Counts are performed at mature level using a mature-only GFF.

    `multi="fractional"` splits each read equally across all matching mature miRNAs.
    """
    logger = _make_logger(log_level)

    # Load mature annotations
    mature_index = load_mature_only_gff(gff_path, logger=logger)

    # Expand BAM patterns
    bam_list = _expand_bam_patterns(bam_paths)
    if not bam_list:
        logger.error("No BAMs found.")
        return 1

    logger.info(
        f"{len(bam_list)} BAM(s) to process; mode={mode}, shift={shift}, "
        f"max_nh={max_nh}, multi={multi}"
    )
    sample_names = [Path(b).stem for b in bam_list]

    per_sample_counts: List[Dict[str, List[float]]] = []
    per_sample_aligned: List[int] = []

    for b in bam_list:
        try:
            mc, aligned = _map_one_bam(
                b,
                mature_index,
                shift=shift,
                max_nh=max_nh,
                mode=mode,
                multi=multi,
                logger=logger,
                log_reads=log_reads,
            )
            per_sample_counts.append(mc)
            per_sample_aligned.append(aligned)
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

    # Feature set = all mature names observed
    features = sorted({m for mc in per_sample_counts for m in mc.keys()})
    if not features:
        logger.warning(
            "No mature miRNAs matched any BAM. Check genome build (chr1 vs 1), "
            "strand orientation, or shift/NH/multi settings."
        )

    # Metric selection
    metric = metric.lower()
    cpm = metric.endswith("_cpm")
    base = metric.replace("_cpm", "")
    idx = {"exact": 0, "approx": 1, "nonspecific": 2}.get(base)
    if idx is None:
        logger.error(
            "--metric must be one of: exact, approx, nonspecific, "
            "exact_cpm, approx_cpm, nonspecific_cpm"
        )
        return 2

    # Write matrix
    outp = Path(out_path)
    outp.parent.mkdir(parents=True, exist_ok=True)

    with open(outp, "w", encoding="utf-8") as fh:
        fh.write("feature\t" + "\t".join(sample_names) + "\n")
        for feat in features:
            vals: List[str] = []
            for s_idx, mc in enumerate(per_sample_counts):
                cnt = mc.get(feat, [0.0, 0.0, 0.0])[idx]
                if cpm:
                    aligned = per_sample_aligned[s_idx] or 1
                    val = 1_000_000 * cnt / aligned
                    vals.append(f"{val:.4f}")
                else:
                    vals.append(f"{cnt:.4f}")
            fh.write(feat + "\t" + "\t".join(vals) + "\n")

    logger.info(f"Wrote matrix to {outp} with {len(features)} features and {len(sample_names)} samples")
    return 0

