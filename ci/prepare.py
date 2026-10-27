"""Download one dataset's public upstream split files, check each MD5 against ci/sources.json,
and write train/test inputs for the action. Runs only on GitHub Actions runners.

usage: python ci/prepare.py <dataset> <outdir>
"""
import csv
import gzip
import hashlib
import io
import json
import os
import pickle
import sys
import time
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))


def fetch(url, tries=4):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "twin-check-ci/1.0"})
            with urllib.request.urlopen(req, timeout=600) as r:
                return r.read()
        except Exception as exc:  # network flake on a shared runner
            last = exc
            time.sleep(10 * (k + 1))
    raise SystemExit(f"download failed after {tries} tries: {url}: {last}")


def get(spec, name):
    data = fetch(spec["url"])
    md5 = hashlib.md5(data).hexdigest()
    ok = md5 == spec["md5"]
    print(f"MD5 {name}: got {md5} expected {spec['md5']} {'OK' if ok else 'MISMATCH'}")
    if not ok:
        raise SystemExit(f"MD5 mismatch for {name} ({spec['url']})")
    return data


def write_fasta(path, seqs):
    with open(path, "w", encoding="utf-8") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">r{i}\n{s}\n")


def write_csv(path, seqs):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sequence"])
        for s in seqs:
            w.writerow([s])


def lrgb_rows(files):
    raw = get(files["data"], "lrgb data")
    rows = list(csv.DictReader(io.StringIO(gzip.decompress(raw).decode("utf-8", "replace"))))
    return [r["peptide_seq"] for r in rows]


def main():
    name, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    src = json.load(open(os.path.join(HERE, "sources.json")))["datasets"][name]
    files = src["files"]
    fmt = "fasta"
    if name == "lrgb-peptides-func":
        import numpy as np
        seqs = lrgb_rows(files)
        sp = pickle.loads(get(files["split"], "lrgb split"))
        tr = [seqs[i] for i in np.asarray(sp["train"]).ravel()]
        te = [seqs[i] for i in np.asarray(sp["test"]).ravel()]
    elif name == "lrgb-corrected-v1":
        seqs = lrgb_rows(files)
        sp = json.loads(get(files["split"], "lrgb corrected split"))["splits"]
        tr = [seqs[int(i)] for i in sp["train"]]
        te = [seqs[int(i)] for i in sp["test"]]
    elif name == "amplify-positives":
        # Upstream FASTA is passed to the action unchanged.
        for k in ("train", "test"):
            with open(os.path.join(out, f"{k}.fasta"), "wb") as fh:
                fh.write(get(files[k], f"amplify {k}"))
        print(f"rows: upstream FASTA written unchanged")
        _meta(out, name, src, "fasta")
        return
    elif name == "flip-meltome-mixed":
        z = zipfile.ZipFile(io.BytesIO(get(files["zip"], "flip meltome splits.zip")))
        names = [n for n in z.namelist() if n.endswith("mixed_split.csv")]
        assert len(names) == 1, names
        rows = list(csv.DictReader(io.StringIO(z.read(names[0]).decode("utf-8", "replace"))))
        tr = [r["sequence"] for r in rows if r["set"] == "train"]
        te = [r["sequence"] for r in rows if r["set"] == "test"]
        fmt = "csv"
    else:
        raise SystemExit(f"unknown dataset {name}")
    print(f"rows: train {len(tr)} test {len(te)}")
    if fmt == "csv":
        write_csv(os.path.join(out, "train.csv"), tr)
        write_csv(os.path.join(out, "test.csv"), te)
    else:
        write_fasta(os.path.join(out, "train.fasta"), tr)
        write_fasta(os.path.join(out, "test.fasta"), te)
    _meta(out, name, src, fmt)


def _meta(out, name, src, fmt):
    gh = os.environ.get("GITHUB_OUTPUT")
    ext = "csv" if fmt == "csv" else "fasta"
    lines = [f"format={fmt}", f"train={out}/train.{ext}", f"test={out}/test.{ext}",
             f"label={src['label']}", f"record_id={src['record_id'] if src['record_id'] is not None else ''}"]
    if gh:
        with open(gh, "a") as fh:
            fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
