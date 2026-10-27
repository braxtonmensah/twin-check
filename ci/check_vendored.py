"""Vendored-code control.

Downloads thirdparty/verify_more.py from braxtonmensah/silicyte-design-loops at the pinned commit,
extracts function twin_rate with ast.get_source_segment, hashes it, and checks that
(1) it equals the sha256 pinned in twin_check.py and
(2) the copy inside twin_check.py is byte-identical.
"""
import ast
import hashlib
import os
import re
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def segment(src, name):
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(src, node)
    raise SystemExit(f"function {name} not found")


def main():
    mine = open(os.path.join(ROOT, "twin_check.py"), encoding="utf-8").read()
    pinned = re.search(r'^VENDORED_SHA256 = "([0-9a-f]{64})"', mine, re.M).group(1)
    commit = re.search(r'^VENDORED_COMMIT = "([0-9a-f]{40})"', mine, re.M).group(1)
    url = (f"https://raw.githubusercontent.com/braxtonmensah/silicyte-design-loops/"
           f"{commit}/thirdparty/verify_more.py")
    orig = urllib.request.urlopen(url, timeout=120).read().decode("utf-8")
    h_orig = hashlib.sha256(segment(orig, "twin_rate").encode("utf-8")).hexdigest()
    h_mine = hashlib.sha256(segment(mine, "twin_rate").encode("utf-8")).hexdigest()
    print(f"pinned   {pinned}\noriginal {h_orig}  ({url})\nvendored {h_mine}  (twin_check.py)")
    if not (pinned == h_orig == h_mine):
        raise SystemExit("FAIL: vendored metric function does not match the pinned original")
    gs = os.environ.get("GITHUB_STEP_SUMMARY")
    if gs:
        with open(gs, "a") as fh:
            fh.write(f"Vendored `twin_rate` sha256 `{h_orig}` matches the original at "
                     f"silicyte-design-loops@{commit[:7]} and the pinned value. PASS\n")
    print("PASS")


if __name__ == "__main__":
    sys.exit(main())
