import json
import os
import random
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import twin_check as tc  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
FIX = os.path.join(ROOT, "tests", "fixtures")


@pytest.fixture(autouse=True)
def _no_github_files(monkeypatch):
    monkeypatch.delenv("GITHUB_OUTPUT", raising=False)
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)


def rand_seqs(n, length, seed):
    rng = random.Random(seed)
    return ["".join(rng.choice(AA) for _ in range(length)) for _ in range(n)]


def write_fasta(path, seqs, wrap=None):
    with open(path, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(f">r{i}\n")
            if wrap:
                for k in range(0, len(s), wrap):
                    fh.write(s[k:k + wrap] + "\n")
            else:
                fh.write(s + "\n")


def run(tmp_path, train, test, *extra):
    tr, te = tmp_path / "train.fa", tmp_path / "test.fa"
    write_fasta(tr, train)
    write_fasta(te, test)
    rec = tmp_path / "rec.json"
    code = tc.main(["--train", str(tr), "--test", str(te), "--record-path", str(rec),
                    "--badge-path", str(tmp_path / "badge.json"),
                    "--summary-path", str(tmp_path / "summary.md"), *extra])
    return code, json.loads(rec.read_text()) if rec.exists() else None


# ---------------------------------------------------------------- must-fail controls
def test_identical_sets_are_all_twins_and_all_exact(tmp_path):
    s = rand_seqs(100, 20, 1)
    code, rec = run(tmp_path, s, s)
    assert code == 0
    assert rec["twin_at_080"] == 100.00 and rec["exact_overlap"] == 100.00


def test_disjoint_random_sets_give_zero(tmp_path):
    code, rec = run(tmp_path, rand_seqs(300, 20, 2), rand_seqs(300, 20, 3))
    assert rec["twin_at_080"] == 0.00 and rec["exact_overlap"] == 0.00


def test_one_substitution_mutants_are_twins_not_exact(tmp_path):
    train = rand_seqs(200, 20, 4)
    rng = random.Random(5)
    test = []
    for s in train:
        p = rng.randrange(20)
        c = rng.choice([a for a in AA if a != s[p]])
        m = s[:p] + c + s[p + 1:]
        assert tc.fuzz.ratio(s, m) == 95.0  # ratio 0.95
        test.append(m)
    code, rec = run(tmp_path, train, test)
    assert rec["twin_at_080"] == 100.00 and rec["exact_overlap"] == 0.00


def test_boundary_079_not_counted_080_counted(tmp_path):
    base = "A" * 100
    at079 = "A" * 79 + "C" * 21   # LCS 79, indel 42, ratio 158/200 = 0.79
    at080 = "A" * 80 + "C" * 20   # LCS 80, indel 40, ratio 160/200 = 0.80
    assert tc.fuzz.ratio(base, at079) == pytest.approx(79.0)
    assert tc.fuzz.ratio(base, at080) == pytest.approx(80.0)
    _, rec = run(tmp_path, [base], [at079])
    assert rec["twin_at_080"] == 0.00
    _, rec = run(tmp_path, [base], [at080])
    assert rec["twin_at_080"] == 100.00


def test_early_exit_reference_agrees_with_vendored_path():
    rng = random.Random(9)
    train = rand_seqs(150, 18, 10) + rand_seqs(50, 30, 11)
    test = []
    for _ in range(120):
        s = rng.choice(train)
        k = rng.randrange(0, 6)
        s = list(s)
        for _ in range(k):
            s[rng.randrange(len(s))] = rng.choice(AA)
        if rng.random() < 0.3:
            s = s[: rng.randrange(10, len(s) + 1)]
        test.append("".join(s))
    for thr in (0.70, 0.80, 0.90):
        w, e, n = tc.measure(train, test, thr)
        ref, _ = tc.twin_rate_early_exit(train, test, thr)
        assert w == ref, (thr, w, ref)
    tc.measure(train, test, 0.80)


def test_length_prefilter_is_lossless():
    rng = random.Random(12)
    for _ in range(3000):
        a = "".join(rng.choice("ACG") for _ in range(rng.randrange(1, 25)))
        b = "".join(rng.choice("ACG") for _ in range(rng.randrange(1, 25)))
        if tc.fuzz.ratio(a, b) >= 80.0:
            assert tc.length_bound_ok(len(a), len(b), 0.80)


# ---------------------------------------------------------------- outputs
def test_badge_matches_shields_endpoint_schema(tmp_path):
    run(tmp_path, rand_seqs(20, 20, 1), rand_seqs(20, 20, 1))
    b = json.loads((tmp_path / "badge.json").read_text())
    assert b["schemaVersion"] == 1
    for k in ("label", "message", "color"):
        assert isinstance(b[k], str) and b[k]
    assert b["message"] == "100.00%"
    assert b["color"] == "red"


def test_color_bands():
    assert tc.color_for(0.0) == "brightgreen"
    assert tc.color_for(0.99) == "brightgreen"
    assert tc.color_for(1.0) == "yellow"
    assert tc.color_for(9.99) == "yellow"
    assert tc.color_for(10.0) == "orange"
    assert tc.color_for(30.0) == "red"


def test_record_fields(tmp_path):
    _, rec = run(tmp_path, rand_seqs(20, 20, 1), rand_seqs(20, 20, 2))
    for k in ("schema_version", "twin_at_080", "exact_overlap", "n_test", "n_train",
              "md5_train", "md5_test", "rule", "tool_version"):
        assert k in rec
    assert rec["schema_version"] == "1.0"


def test_summary_has_table_header(tmp_path):
    run(tmp_path, rand_seqs(20, 20, 1), rand_seqs(20, 20, 2))
    assert "| metric | value |" in (tmp_path / "summary.md").read_text()


def test_fail_above(tmp_path):
    s = rand_seqs(30, 20, 1)
    code, _ = run(tmp_path, s, s, "--fail-above", "1")
    assert code == 1
    code, _ = run(tmp_path, s, s, "--fail-above", "100")
    assert code == 0


def test_example_fixture_is_40_percent(tmp_path):
    rec = tmp_path / "r.json"
    code = tc.main(["--train", os.path.join(FIX, "example_train.fasta"),
                    "--test", os.path.join(FIX, "example_test.fasta"),
                    "--record-path", str(rec), "--badge-path", "", "--summary-path", ""])
    r = json.loads(rec.read_text())
    assert code == 0 and r["n_test"] == 50 and r["twin_at_080"] == 40.00 and r["exact_overlap"] == 0.00


# ---------------------------------------------------------------- input parsing
def test_wrapped_fasta_is_joined(tmp_path):
    s = rand_seqs(10, 130, 7)
    tr, te = tmp_path / "tr.fa", tmp_path / "te.fa"
    write_fasta(tr, s, wrap=60)
    write_fasta(te, s, wrap=None)
    assert tc.read_fasta(str(tr)) == s
    rec = tmp_path / "r.json"
    tc.main(["--train", str(tr), "--test", str(te), "--record-path", str(rec),
             "--badge-path", "", "--summary-path", ""])
    assert json.loads(rec.read_text())["exact_overlap"] == 100.00


def test_csv_named_column(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("id,seq,label\n1,ACDE,0\n2,FGHI,1\n")
    assert tc.read_csv(str(p), "seq") == ["ACDE", "FGHI"]


@pytest.mark.parametrize("content,fmt,col,needle", [
    ("", "fasta", "", "empty"),
    ("   \n\n", "fasta", "", "empty"),
    ("ACDEFG\n", "fasta", "", "before the first"),
    (">a\nACD\n>a\nEFG\n", "fasta", "", "duplicate FASTA header"),
    (">a\n>b\nEFG\n", "fasta", "", "has no sequence"),
    (">a\nACD\n>b\n", "fasta", "", "has no sequence"),
    ("", "csv", "seq", "empty"),
    ("seq,seq\nA,B\n", "csv", "seq", "duplicate CSV header"),
    ("id,seq\n", "csv", "seq", "no data rows"),
    ("id,sequence\n1,AC\n", "csv", "seq", "not found"),
    ("id,seq\n1,\n", "csv", "seq", "empty 'seq'"),
])
def test_bad_inputs_raise_clear_errors(tmp_path, content, fmt, col, needle):
    p = tmp_path / "bad"
    p.write_text(content)
    with pytest.raises(tc.InputError, match=needle):
        tc.load(str(p), fmt, col)


def test_bad_input_exits_nonzero_not_zero_rate(tmp_path, capsys):
    p = tmp_path / "empty.fa"
    p.write_text("")
    good = tmp_path / "g.fa"
    write_fasta(good, ["ACDEFGHIK"])
    code = tc.main(["--train", str(p), "--test", str(good), "--record-path", str(tmp_path / "r.json")])
    assert code == 2
    assert not (tmp_path / "r.json").exists()
    assert "input error" in capsys.readouterr().err
