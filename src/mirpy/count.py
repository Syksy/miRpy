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


@dataclass
class Mature:
    name: str
    start: int  # 1-based inclusive (from GFF)
    end: int    # 1-based inclusive
    precursor_id: str  # ID=MI... it derives from
    precursor_name: Optional[str] = None


@dataclass
class Precursor:
    id: str
    name: Optional[str]
    chr: str
    strand: str
    start: int
    end: int
    matures: List[Mature]  # length 1 or 2 typically

def _make_logger(level: str) -> logging.Logger:
    lvl = getattr(logging, level.upper(), logging.INFO)
    logger = logging.getLogger("mirpy.count")
    # Configure once
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


def load_mirbase_like_gff(
    gff_path: str | Path,
    logger: logging.Logger | None = None,
) -> Tuple[Dict[str, Dict[str, List[Precursor]]], Dict[str, Mature], Dict[str, str]]:
    """Load miRBase-like GFF3 (primary_transcript + miRNA)."""
    prec_by_id: Dict[str, Precursor] = {}
    mature_index: Dict[str, Mature] = {}

    with _open_text_auto(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start_s, end_s, _score, strand, _phase, attrs = cols
            A = _parse_attrs(attrs)
            try:
                start = int(start_s); end = int(end_s)
            except ValueError:
                continue

            if feature == "miRNA_primary_transcript":
                pid = A.get("ID"); pname = A.get("Name")
                if not pid:
                    continue
                if pid not in prec_by_id:
                    prec_by_id[pid] = Precursor(
                        id=pid, name=pname, chr=chrom, strand=strand, start=start, end=end, matures=[]
                    )
            elif feature == "miRNA":
                name = A.get("Name"); derives = A.get("Derives_from")
                if not name or not derives:
                    continue
                m = Mature(name=name, start=start, end=end, precursor_id=derives)
                mature_index[name] = m

    for m in mature_index.values():
        prec = prec_by_id.get(m.precursor_id)
        if prec is not None:
            m.precursor_name = prec.name or prec.id
            prec.matures.append(m)

    premiR_index: Dict[str, Dict[str, List[Precursor]]] = {}
    for prec in prec_by_id.values():
        premiR_index.setdefault(prec.chr, {}).setdefault(prec.strand, []).append(prec)
    for chr_ in premiR_index:
        for st in premiR_index[chr_]:
            premiR_index[chr_][st].sort(key=lambda p: p.start)

    mature_to_precursor_name = {m.name: (m.precursor_name or m.precursor_id) for m in mature_index.values()}

    if logger:
        logger.info(f"GFF loaded: {len(prec_by_id)} precursors; {len(mature_index)} mature miRNAs")
        # Show chromosomes/strands for context
        if logger.isEnabledFor(logging.DEBUG):
            for i, (chr_, strands) in enumerate(premiR_index.items()):
                #if i >= 5: break
                for st, precs in strands.items():
                    logger.debug(f"  chr={chr_} strand={st}: {len(precs)} precursors")

    return premiR_index, mature_index, mature_to_precursor_name

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


def _distance_to_mature(start_1b: int, end_1b: int, m: Mature) -> int:
    return abs(m.start - start_1b) + abs(m.end - end_1b)


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
        matches = glob.glob(pat) if any(ch in pat for ch in "*?[]") else ([pat] if os.path.exists(pat) else [])
        for m in sorted(matches):
            if m not in seen:
                seen.add(m)
                out.append(m)
        if not matches:
            print(f"[WARNING] No BAMs matched: {pat}")
    return out


def _count_one_bam(
    bam_path: str | Path,
    premiR_index: Dict[str, Dict[str, List[Precursor]]],
    *,
    shift: int,
    max_nh: int,
    mode: str = "qname",  # "qname", "nh-bucket", or "auto"
    logger: logging.Logger | None = None,
    log_reads: int = 0
) -> Tuple[Dict[str, List[int]], int]:
    """
    Count per mature for a single BAM.
    Returns (mature_counts, aligned_reads).
    mature_counts[name] -> [exact, approx, nonspec]
    """
    # Open a short-lived handle to inspect sortedness (and that the BAM is readable)
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
            # Sorted, use the fast path
            return _count_qname_sorted(bam_path, premiR_index, shift=shift, max_nh=max_nh, logger=logger, log_reads=log_reads)

        elif mode == "nh-bucket":
            # Stream unsorted, using NH-bucket logic
            return _count_unsorted_nh_bucket(bam_path, premiR_index, shift=shift, max_nh=max_nh, logger=logger, log_reads=log_reads)

        elif mode == "auto":
            if _is_qname_sorted(test):
                return _count_qname_sorted(bam_path, premiR_index, shift=shift, max_nh=max_nh, logger=logger, log_reads=log_reads)
            return _count_unsorted_nh_bucket(bam_path, premiR_index, shift=shift, max_nh=max_nh, logger=logger, log_reads=log_reads)

        else:
            raise RuntimeError(f"Unknown mode: {mode}")

    finally:
        # Always close the probe handle
        try:
            test.close()
        except Exception:
            pass

def _count_qname_sorted(
    bam_path: str | Path,
    premiR_index: Dict[str, Dict[str, List[Precursor]]],
    *,
    shift: int,
    max_nh: int,
    logger: logging.Logger | None = None,
    log_reads: int = 0,
) -> Tuple[Dict[str, List[int]], int]:
    """
    Stream a QNAME-sorted BAM and count matches at the mature level.

    Returns
    -------
    (mature_counts, aligned_reads)
      mature_counts[name] -> [exact, approx, nonspecific]
      aligned_reads: number of mapped alignments observed (CPM denominator upstream)
    """
    def ensure(name: str, table: Dict[str, List[int]]):
        if name not in table:
            table[name] = [0, 0, 0]
    if logger:
        logger.info(f"Counting (QNAME) in {bam_path}")
        if logger.isEnabledFor(logging.DEBUG):
            try:
                with bn.AlignmentFile(str(bam_path), "rb") as btmp:
                    logger.debug(f"BAM references: {list(getattr(btmp, 'references', []))[:10]}")
            except Exception:
                pass

    mature_counts: Dict[str, List[int]] = {}
    aligned_reads = 0

    # Per-read bucket while names are contiguous
    current_qname: Optional[str] = None
    bucket: List = []
    reads_logged = 0  # limit DEBUG detail to the first N reads

    def process_bucket():
        """Apply Perl-like logic to a single read's alignments."""
        nonlocal reads_logged
        nonlocal mature_counts
        if not bucket:
            return

        # NH: genomic matches per read; use first alignment carrying NH if present
        nh = None
        for aln in bucket:
            try:
                nh = aln.opt("NH")
                break
            except Exception:
                pass
        genomic_matches = nh if isinstance(nh, int) else 1
        if genomic_matches > max_nh:
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                logger.debug(f"skip read: NH={genomic_matches} > max_nh={max_nh}")
            return

        sum_matches = 0
        exact_tmp: Dict[str, None] = {}
        approx_tmp: Dict[str, None] = {}
        matched_names: Dict[str, None] = {}

        for aln in bucket:
            if getattr(aln, "is_unmapped", False):
                continue

            chr_ = getattr(aln, "reference_name", None)
            if chr_ is None:
                continue
            strand = "-" if getattr(aln, "is_reverse", False) else "+"
            start_1b, end_1b = _aln_locus_1b(aln)

            # Exact genomic match via NM
            nm = None
            try:
                nm = aln.opt("NM")
            except Exception:
                pass
            non_pm = 0 if (nm == 0) else 1

            pre_list = premiR_index.get(chr_, {}).get(strand, [])
            if not pre_list:
                continue

            # Find overlapping precursor(s), then compare to its mature(s)
            for prec in pre_list:
                if start_1b <= prec.end and end_1b >= prec.start:
                    for m in prec.matures:
                        d = _distance_to_mature(start_1b, end_1b, m)
                        if non_pm == 0 and d == 0:
                            exact_tmp[m.name] = None
                            matched_names[m.name] = None
                            sum_matches += 1
                        elif d <= shift:
                            approx_tmp[m.name] = None
                            matched_names[m.name] = None
                            sum_matches += 1

        if sum_matches == 0:
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                # Inspect why: first aln chrom/strand presence?
                a0 = next((a for a in bucket if not getattr(a, "is_unmapped", False)), None)
                if a0 is not None:
                    chr_ = getattr(a0, "reference_name", None)
                    st  = "-" if getattr(a0, "is_reverse", False) else "+"
                    has_chr = chr_ in premiR_index
                    has_str = has_chr and (st in premiR_index[chr_])
                    logger.debug(f"no match: chr={chr_} in_gff={has_chr}, strand_ok={has_str}, "
                                 f"start/end={_aln_locus_1b(a0)}")
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                reads_logged += 1
            return

        unique = list(matched_names.keys())
        single = (len(unique) == 1)

        # If all genomic hits are mature and agree on a single mature, award exact/approx (+ nonspecific per Perl)
        if sum_matches == genomic_matches and single:
            for k in exact_tmp.keys():
                ensure(k, mature_counts)
                mature_counts[k][0] += 1  # exact
                mature_counts[k][1] += 1  # approx (also incremented for exact)
                mature_counts[k][2] += 1  # nonspecific (Perl increments all three)
            for k in approx_tmp.keys():
                ensure(k, mature_counts)
                mature_counts[k][1] += 1  # approx
                mature_counts[k][2] += 1  # nonspecific
        else:
            # Not all hits mature OR ambiguous w.r.t mature name:
            # If single mature name, count as nonspecific only (as in Perl).
            if single:
                for k in exact_tmp.keys():
                    ensure(k, mature_counts)
                    mature_counts[k][2] += 1
                for k in approx_tmp.keys():
                    ensure(k, mature_counts)
                    mature_counts[k][2] += 1
            # If multiple mature names, then ignore the read entirely (Perl behavior).

    progress_every = 100000
    # Single streaming pass (names are contiguous in QNAME-sorted BAM)
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

    # Flush last bucket
    process_bucket()

    if logger:
        logger.info(f"Done {bam_path}: aligned={aligned_reads}, features_matched={len(mature_counts)}")

    return mature_counts, aligned_reads

def _count_unsorted_nh_bucket(
    bam_path: str | Path,
    premiR_index: Dict[str, Dict[str, List[Precursor]]],
    *,
    shift: int,
    max_nh: int,
    logger: logging.Logger | None = None,
    log_reads: int = 0,
) -> Tuple[Dict[str, List[int]], int]:
    """Stream an unsorted BAM, buffering per read until we've seen NH alignments."""
    if logger:
        logger.info(f"Counting (NH-bucket) in {bam_path}")
        if logger.isEnabledFor(logging.DEBUG):
            try:
                with bn.AlignmentFile(str(bam_path), "rb") as btmp:
                    logger.debug(f"BAM references: {list(getattr(btmp, 'references', []) )[:10]}")
            except Exception:
                pass

    mature_counts: Dict[str, List[int]] = {}
    aligned_reads = 0

    def ensure(name: str):
        if name not in mature_counts:
            mature_counts[name] = [0, 0, 0]

    # per-read buckets
    buckets: Dict[str, List] = defaultdict(list)
    nh_by_read: Dict[str, Optional[int]] = {}
    reads_logged = 0  # limit DEBUG detail to the first N reads

    def process_bucket(qn: str):
        nonlocal reads_logged
        nonlocal mature_counts
        bucket = buckets.pop(qn, [])
        nh = nh_by_read.pop(qn, None)
        genomic_matches = nh if isinstance(nh, int) else len(bucket)
        if genomic_matches > max_nh:
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                logger.debug(f"{qn}: skip NH={genomic_matches} > max_nh={max_nh}")
            return

        sum_matches = 0
        exact_tmp: Dict[str, None] = {}
        approx_tmp: Dict[str, None] = {}
        matched_names: Dict[str, None] = {}

        for aln in bucket:
            if getattr(aln, "is_unmapped", False):
                continue
            chr_ = getattr(aln, "reference_name", None)
            if chr_ is None:
                continue
            strand = "-" if getattr(aln, "is_reverse", False) else "+"
            start_1b, end_1b = _aln_locus_1b(aln)
            nm = None
            try:
                nm = aln.opt("NM")
            except Exception:
                pass
            non_pm = 0 if (nm == 0) else 1

            pre_list = premiR_index.get(chr_, {}).get(strand, [])
            for prec in pre_list:
                if start_1b <= prec.end and end_1b >= prec.start:
                    for m in prec.matures:
                        d = _distance_to_mature(start_1b, end_1b, m)
                        if non_pm == 0 and d == 0:
                            exact_tmp[m.name] = None; matched_names[m.name] = None; sum_matches += 1
                        elif d <= shift:
                            approx_tmp[m.name] = None; matched_names[m.name] = None; sum_matches += 1

        if sum_matches == 0:
            if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
                reads_logged += 1
                a0 = next((a for a in bucket if not getattr(a, "is_unmapped", False)), None)
                if a0 is not None:
                    chr_ = getattr(a0, "reference_name", None)
                    st = "-" if getattr(a0, "is_reverse", False) else "+"
                    has_chr = chr_ in premiR_index
                    has_str = has_chr and (st in premiR_index.get(chr_, {}))
                    logger.debug(f"{qn}: no match: chr={chr_} in_gff={has_chr}, strand_ok={has_str}, "
                                 f"start/end={_aln_locus_1b(a0)}")
            return

        unique = list(matched_names.keys())
        single = (len(unique) == 1)

        if sum_matches == genomic_matches and single:
            for k in exact_tmp.keys():
                ensure(k); mature_counts[k][0]+=1; mature_counts[k][1]+=1; mature_counts[k][2]+=1
            for k in approx_tmp.keys():
                ensure(k); mature_counts[k][1]+=1; mature_counts[k][2]+=1
        else:
            if single:
                for k in exact_tmp.keys():
                    ensure(k); mature_counts[k][2]+=1
                for k in approx_tmp.keys():
                    ensure(k); mature_counts[k][2]+=1
        # if multiple mature names, then ignore (as before)
        if logger and logger.isEnabledFor(logging.DEBUG) and reads_logged < log_reads:
            reads_logged += 1

    progress_every = 100000
    # Stream once
    for aln in bn.AlignmentFile(str(bam_path), "rb"):
        if getattr(aln, "is_unmapped", False):
            continue
        aligned_reads += 1

        if logger and aligned_reads % progress_every == 0:
            logger.info(f"Processed {aligned_reads:,} alignments...")

        qn = _get_read_name(aln)
        buckets[qn].append(aln)

        # learn NH for this read if present
        if qn not in nh_by_read or nh_by_read[qn] is None:
            try:
                nh_by_read[qn] = aln.opt("NH")
            except Exception:
                nh_by_read.setdefault(qn, None)

        # if we know NH and have seen that many alignments for this read, process now
        nh = nh_by_read.get(qn)
        if isinstance(nh, int) and len(buckets[qn]) >= nh:
            process_bucket(qn)

    # Flush leftovers (reads without NH or where we didn’t hit nh alignments)
    for qn in list(buckets.keys()):
        process_bucket(qn)

    if logger:
        logger.info(f"Done {bam_path}: aligned={aligned_reads}, features_matched={len(mature_counts)}")

    return mature_counts, aligned_reads

def count_matrix(
    bam_paths: List[str],
    gff_path: str | Path,
    out_path: str | Path,
    *,
    level: str = "mature",  # "mature" or "precursor"
    metric: str = "exact",  # exact | approx | nonspecific | exact_cpm | approx_cpm | nonspecific_cpm
    shift: int = 4,
    max_nh: int = 50,
    mode: str = "auto",
    log_level: str = "INFO",
    log_reads: int = 0,
) -> int:
    """
    Build a matrix with rows = features (mature or precursor), columns = samples.
    Metric selects which value is reported (counts or CPM).
    """
    # Load logger
    logger = _make_logger(log_level)
    # Load annotations
    premiR_index, mature_index, mature_to_prec = load_mirbase_like_gff(gff_path, logger=logger)

    # Expand patterns and collect sample names
    bam_list = _expand_bam_patterns(bam_paths)
    if not bam_list:
        logger.error("No BAMs found.")
        return 1

    logger.info(f"{len(bam_list)} BAM(s) to process; mode={mode}, shift={shift}, max_nh={max_nh}")
    sample_names = [Path(b).stem for b in bam_list]

    # Count per sample (at mature level)
    per_sample_counts: List[Dict[str, List[int]]] = []
    per_sample_aligned: List[int] = []

    for b in bam_list:
        try:
            mc, aligned = _count_one_bam(
                b,
                premiR_index,
                shift=shift,
                max_nh=max_nh,
                mode=mode,
                logger=logger,
                log_reads=log_reads,
            )
            per_sample_counts.append(mc)
            per_sample_aligned.append(aligned)
        except Exception as e:
            logger.error(f"{b}: {e}")
            return 1

    # Feature set (rows)
    if level == "mature":
        features = sorted({m for mc in per_sample_counts for m in mc.keys()})
        if not features:
            logger.warning(
                "No features matched any BAM. Common causes: contig prefix mismatch (chr1 vs 1), "
                "wrong strand, or too-strict shift/NH filters."
            )
    elif level == "precursor":
        # derive precursor names from observed matures
        precs = set()
        for mc in per_sample_counts:
            for mname in mc.keys():
                precs.add(mature_to_prec.get(mname, None))
        features = sorted(x for x in precs if x)
    else:
        logger.error("--level must be 'mature' or 'precursor'")
        return 2

    # Decide metric slot and CPM flag
    metric = metric.lower()
    cpm = metric.endswith("_cpm")
    base = metric.replace("_cpm", "")
    idx = {"exact": 0, "approx": 1, "nonspecific": 2}.get(base)
    if idx is None:
        logger.error("--metric must be one of: exact, approx, nonspecific, exact_cpm, approx_cpm, nonspecific_cpm")
        return 2

    # Assemble matrix
    outp = Path(out_path)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, "w", encoding="utf-8") as fh:
        fh.write("feature\t" + "\t".join(sample_names) + "\n")
        for feat in features:
            row_vals: List[str] = []
            for s_idx, mc in enumerate(per_sample_counts):
                if level == "mature":
                    cnt = mc.get(feat, [0, 0, 0])[idx]
                else:  # precursor
                    # sum over matures that derive from this precursor
                    cnt = 0
                    for mname, triple in mc.items():
                        if mature_to_prec.get(mname, None) == feat:
                            cnt += triple[idx]
                if cpm:
                    aligned = per_sample_aligned[s_idx] or 1
                    val = 1_000_000 * cnt / aligned
                    row_vals.append(f"{val:.2f}")
                else:
                    row_vals.append(str(cnt))

            fh.write(feat + "\t" + "\t".join(row_vals) + "\n")

    logger.info(f"Wrote matrix to {outp} with {len(features)} features and {len(sample_names)} samples")

    return 0
