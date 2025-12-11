from __future__ import annotations

import gzip
import io
import os
import urllib.request
from pathlib import Path

MIRBASE_CURRENT_GFF_URL = "https://mirbase.org/download/hsa.gff3"

def _is_gzip_file(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

def _first_line_text(path: Path) -> str:
    # read first line whether gzip or plain
    if _is_gzip_file(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.readline()
    else:
        with open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.readline()

def download_mirbase_gff(
        dest: str | os.PathLike,
        *,
        force: str | None = True,
        url: str | None = None,
) -> Path:
    """
    Download a GFF3 (gz or plain) to 'dest'. If dest ends with .gff3 and the
    download is gzipped, auto-decompress; otherwise keep as-is.
    """
    url = url or MIRBASE_CURRENT_GFF_URL
    out = Path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and not force:
        raise FileExistsError(f"Destination exists: {out}")

    tmp = out.with_suffix(out.suffix + ".partial")
    try:
        # stream download
        with urllib.request.urlopen(url) as resp, open(tmp, "wb") as fh:
            while True:
                chunk = resp.read(65536)
                if not chunk:
                    break
                fh.write(chunk)

        is_gz = _is_gzip_file(tmp)

        # Decide how to place the file based on BOTH dest suffix and actual content
        if out.suffix.lower() == ".gff3":
            if is_gz:
                # auto-decompress into final .gff3
                with gzip.open(tmp, "rb") as gz, open(out, "wb") as outfh:
                    outfh.write(gz.read())
                tmp.unlink(missing_ok=True)
            else:
                # already plain text -> just move/rename
                if out.exists():
                    out.unlink()
                tmp.rename(out)
        else:
            # any other suffix (including .gz): move verbatim, regardless of gzip-ness
            if out.exists():
                out.unlink()
            tmp.rename(out)

        # Sanity-check content: first line should include GFF3 header
        first = _first_line_text(out).strip()
        if not first.startswith("##gff-version 3"):
            # common failure: server returned HTML page (starts with "<")
            if first.startswith("<"):
                raise RuntimeError("Server likely returned HTML instead of GFF3 (wrong URL or needs a direct file link).")
            raise RuntimeError(f"Downloaded file does not look like GFF3 (first line: {first!r}).")

        return out

    except Exception as e:
        try:
            tmp.unlink(missing_ok=True)
        except Exception:
            pass
        raise RuntimeError(f"Download failed: {e}") from e
