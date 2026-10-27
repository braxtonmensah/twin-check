"""twin-check: near-twin leakage between a training split and a test split.

Rule (unchanged from Silicyte's published leakage record,
https://silicyte.com/data/leakage-measurements.json):

    twin rate  = fraction of test rows whose best training match has a Levenshtein ratio
                 of 0.80 or above, where ratio = (len(a) + len(b) - indel(a, b)) / (len(a) + len(b))
                 (rapidfuzz.fuzz.ratio / 100, the python-Levenshtein ratio)
    exact rate = fraction of test rows whose sequence appears verbatim in the training split

VENDORED METRIC CODE
    source repo   https://github.com/braxtonmensah/silicyte-design-loops
    source path   thirdparty/verify_more.py, function twin_rate
    source commit 50a6eb704730607abd02bf381c811c68ebd4d378
    sha256 of the original function source (ast.get_source_segment, UTF-8):
                  a01b9ee7f3b5c034bc6f69dc38c7a17c4d068cd21144b1574ebb8a89c56d210f
    The function `twin_rate` below is a byte-for-byte copy. CI (ci/check_vendored.py) downloads
    the original at that commit, re-hashes it, and checks both copies against the value above.
    The library is rapidfuzz, as imported by thirdparty/verify_thirdparty.py and verify_more.py
    at that commit. Those scripts do not pin a version; this action pins rapidfuzz 3.13.0.

Built by Silicyte (https://silicyte.com) with AI assistance. MIT license.
"""
import argparse
import csv
import hashlib
import json
import os
import sys

import numpy as np  # noqa: F401  (cdist returns a numpy array; imported so a missing install fails early)
from rapidfuzz import fuzz, process

TOOL_VERSION = "1.0.2"
VENDORED_SHA256 = "a01b9ee7f3b5c034bc6f69dc38c7a17c4d068cd21144b1574ebb8a89c56d210f"
VENDORED_COMMIT = "50a6eb704730607abd02bf381c811c68ebd4d378"
RULE = ("test row counts as a twin if any training sequence has Levenshtein ratio "
        "(rapidfuzz.fuzz.ratio/100) >= {thr:.2f}; exact = verbatim sequence match")

THRESH = 80.0          # rapidfuzz.fuzz.ratio is 100x the Levenshtein ratio; set from --threshold


class InputError(Exception):
    pass


# ----------------------------------------------------------------- vendored, do not edit
def twin_rate(train, test, workers=-1, block=256):
    """Fraction of test rows with any train sequence at ratio >= THRESH."""
    tr = list(dict.fromkeys(train))
    te = list(test)
    hits = 0
    for i in range(0, len(te), block):
        chunk = te[i:i + block]
        sim = process.cdist(chunk, tr, scorer=fuzz.ratio,
                            score_cutoff=THRESH, workers=workers)
        hits += int((sim >= THRESH).any(axis=1).sum())
    return hits, len(te)
# ----------------------------------------------------------------- end vendored


def exact_overlap(train, test):
    """Exact matches through a hash set."""
    trset = set(train)
    return sum(1 for t in test if t in trset), len(test)


def length_bound_ok(la, lb, thr):
    """Lossless prefilter. indel(a, b) >= |la - lb|, so ratio <= 2*min(la, lb)/(la + lb).
    If that upper bound is below thr, the pair cannot be a twin."""
    if la + lb == 0:
        return True
    return 2 * min(la, lb) / (la + lb) >= thr


def twin_rate_early_exit(train, test, thr):
    """Independent all-pairs search with early exit on the first hit at or above thr and the
    lossless length prefilter above. Slower; used by the tests to cross-check the vendored path."""
    tr = list(dict.fromkeys(train))
    cutoff = thr * 100.0
    hits = 0
    for t in test:
        lt = len(t)
        for s in tr:
            if not length_bound_ok(lt, len(s), thr):
                continue
            if fuzz.ratio(t, s, score_cutoff=cutoff) >= cutoff:
                hits += 1
                break
    return hits, len(test)


# ----------------------------------------------------------------- input parsing
def read_fasta(path):
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        text = fh.read()
    if not text.strip():
        raise InputError(f"{path}: file is empty")
    seqs, headers, cur, header = [], set(), [], None
    for lineno, line in enumerate(text.splitlines(), 1):
        if line.startswith(">"):
            if header is not None:
                if not cur:
                    raise InputError(f"{path}: record '{header}' has no sequence (line {lineno})")
                seqs.append("".join(cur))
            header = line[1:].strip()
            if header in headers:
                raise InputError(f"{path}: duplicate FASTA header '{header}' (line {lineno})")
            headers.add(header)
            cur = []
        elif line.strip():
            if header is None:
                raise InputError(f"{path}: sequence data before the first '>' header "
                                 f"(line {lineno}); is this really FASTA?")
            cur.append(line.strip())
    if header is None:
        raise InputError(f"{path}: no FASTA records found")
    if not cur:
        raise InputError(f"{path}: record '{header}' has no sequence")
    seqs.append("".join(cur))
    return seqs


def read_csv(path, column):
    if not column:
        raise InputError("format is csv but no column was given")
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        text = fh.read()
    if not text.strip():
        raise InputError(f"{path}: file is empty")
    reader = csv.reader(text.splitlines())
    header = next(reader)
    dupes = sorted({h for h in header if header.count(h) > 1})
    if dupes:
        raise InputError(f"{path}: duplicate CSV header(s) {dupes}")
    if column not in header:
        raise InputError(f"{path}: column '{column}' not found; columns are {header}")
    k = header.index(column)
    seqs = []
    for rowno, row in enumerate(reader, 2):
        if not row:
            continue
        if k >= len(row) or row[k] == "":
            raise InputError(f"{path}: empty '{column}' value on row {rowno}")
        seqs.append(row[k])
    if not seqs:
        raise InputError(f"{path}: header only, no data rows")
    return seqs


def load(path, fmt, column):
    if not path or not os.path.isfile(path):
        raise InputError(f"input file not found: '{path}'")
    if fmt == "fasta":
        return read_fasta(path)
    if fmt == "csv":
        return read_csv(path, column)
    raise InputError(f"unknown format '{fmt}' (use fasta or csv)")


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def color_for(pct):
    """Presentation only, not a rating. Bands are documented in the README."""
    if pct < 1.0:
        return "brightgreen"
    if pct < 10.0:
        return "yellow"
    if pct < 30.0:
        return "orange"
    return "red"


def measure(train, test, thr):
    global THRESH
    THRESH = round(thr * 100.0, 10)
    w, n = twin_rate(train, test)
    e, _ = exact_overlap(train, test)
    return w, e, n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--train", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--format", default="fasta")
    ap.add_argument("--column", default="")
    ap.add_argument("--threshold", default="0.80")
    ap.add_argument("--fail-above", default="")
    ap.add_argument("--badge-path", default=".twin-check/badge.json")
    ap.add_argument("--record-path", default=".twin-check/record.json")
    ap.add_argument("--summary-path", default=".twin-check/summary.md")
    ap.add_argument("--label", default="test twins @0.80")
    a = ap.parse_args(argv)

    try:
        thr = float(a.threshold)
        if not 0.0 < thr <= 1.0:
            raise InputError(f"threshold must be in (0, 1], got {a.threshold}")
        fail_above = float(a.fail_above) if str(a.fail_above).strip() else None
        fmt = a.format.strip().lower()
        train = load(a.train, fmt, a.column)
        test = load(a.test, fmt, a.column)
    except (InputError, ValueError) as exc:
        print(f"::error title=twin-check input error::{exc}")
        print(f"twin-check: input error: {exc}", file=sys.stderr)
        return 2

    w, e, n = measure(train, test, thr)
    twin_pct = 100.0 * w / n
    exact_pct = 100.0 * e / n
    twin_s, exact_s = f"{twin_pct:.2f}", f"{exact_pct:.2f}"

    record = {
        "schema_version": "1.0",
        "tool": "twin-check",
        "tool_version": TOOL_VERSION,
        "rule": RULE.format(thr=thr),
        "threshold_levenshtein_ratio": thr,
        "twin_at_080": round(twin_pct, 2),
        "exact_overlap": round(exact_pct, 2),
        "twin_rows": w,
        "exact_rows": e,
        "n_test": n,
        "n_train": len(train),
        "n_train_unique": len(set(train)),
        "md5_train": md5_file(a.train),
        "md5_test": md5_file(a.test),
        "library": "rapidfuzz 3.13.0",
        "vendored_from": f"braxtonmensah/silicyte-design-loops@{VENDORED_COMMIT}:thirdparty/verify_more.py#twin_rate",
        "vendored_sha256": VENDORED_SHA256,
        "note": "A twin rate is not a score inflation; nothing is retrained.",
    }
    if abs(thr - 0.80) > 1e-12:
        record["twin_at_080"] = None
        record["twin_at_threshold"] = round(twin_pct, 2)
    label = str(a.label)
    if label == "test twins @0.80" and abs(thr - 0.80) > 1e-12:
        label = f"test twins @{thr:.2f}"   # do not show the default 0.80 label for another threshold
    badge = {"schemaVersion": 1, "label": label, "message": f"{twin_s}%",
             "color": color_for(twin_pct)}

    summary = "\n".join([
        "### twin-check",
        "",
        "| metric | value |",
        "|---|---|",
        f"| test rows with a training row at Levenshtein ratio >= {thr:.2f} | {w} of {n} ({twin_s}%) |",
        f"| test rows with an exact training match | {e} of {n} ({exact_s}%) |",
        f"| training rows (unique) | {len(train)} ({len(set(train))}) |",
        f"| md5 train / test | `{record['md5_train']}` / `{record['md5_test']}` |",
        "",
        "A twin rate is not a score inflation. Sequence identity is not activity.",
        "",
    ])

    for p, payload in ((a.badge_path, json.dumps(badge)),
                       (a.record_path, json.dumps(record, indent=1)),
                       (a.summary_path, summary)):
        if p:
            d = os.path.dirname(p)
            if d:
                os.makedirs(d, exist_ok=True)
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(payload + ("\n" if not payload.endswith("\n") else ""))

    gs = os.environ.get("GITHUB_STEP_SUMMARY")
    if gs:
        with open(gs, "a", encoding="utf-8") as fh:
            fh.write(summary + "\n")
    go = os.environ.get("GITHUB_OUTPUT")
    if go:
        with open(go, "a", encoding="utf-8") as fh:
            fh.write(f"twin_rate={twin_s}\nexact_rate={exact_s}\n"
                     f"n_train={len(train)}\nn_test={n}\n")
    print(json.dumps(record))

    if fail_above is not None and twin_pct > fail_above:
        print(f"::error title=twin-check::twin rate {twin_s}% is above fail-above {fail_above}%")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
