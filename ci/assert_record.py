"""Assert the action's outputs equal the published record to two decimals.

usage: python ci/assert_record.py <dataset> <twin_rate> <exact_rate> <n_train> <n_test>
"""
import json
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    name, twin, exact, n_train, n_test = sys.argv[1:6]
    src = json.load(open(os.path.join(HERE, "sources.json")))
    ds = src["datasets"][name]
    failures = []
    rows = [("twin rate", twin, None), ("exact rate", exact, None)]
    if ds["record_id"] is None:
        want_twin = f"{ds['expected_twin_pct']:.2f}"
        want = {"twin": want_twin, "exact": "0.00"}
        where = "corrected-split control (expected 0.00, the file's own rebuilt_twin_pct)"
        want_n_test = None
        want_n_train = None
    else:
        req = urllib.request.Request(src["record_url"], headers={"User-Agent": "twin-check-ci/1.0"})
        record = json.load(urllib.request.urlopen(req, timeout=120))
        m = [x for x in record["measurements"] if x["id"] == ds["record_id"]]
        assert len(m) == 1, f"record id {ds['record_id']} not found"
        m = m[0]
        want = {"twin": f"{m['twin_at_080_pct']:.2f}", "exact": f"{m['exact_overlap_pct']:.2f}"}
        want_n_test = m.get("test_rows")
        want_n_train = m.get("train_rows")
        where = f"{src['record_url']} id {ds['record_id']} ({m['benchmark']}, {m['split']})"

    print(f"comparing against {where}")
    lines = ["| check | action | record | |", "|---|---:|---:|---|"]

    def cmp(label, got, exp):
        ok = str(got) == str(exp)
        lines.append(f"| {label} | {got} | {exp} | {'PASS' if ok else 'FAIL'} |")
        if not ok:
            failures.append(f"{label}: action {got} != record {exp}")

    cmp("twin_at_080 %", twin, want["twin"])
    cmp("exact_overlap %", exact, want["exact"])
    if want_n_test is not None:
        cmp("n_test", n_test, want_n_test)
    if want_n_train is not None:
        cmp("n_train", n_train, want_n_train)

    table = f"#### {name}: {where}\n\n" + "\n".join(lines) + "\n"
    print(table)
    gs = os.environ.get("GITHUB_STEP_SUMMARY")
    if gs:
        with open(gs, "a") as fh:
            fh.write(table + "\n")
    if failures:
        raise SystemExit("; ".join(failures))
    print("PASS")


if __name__ == "__main__":
    main()
