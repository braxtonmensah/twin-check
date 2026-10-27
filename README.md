# twin-check

[![reproduce-record](https://github.com/braxtonmensah/twin-check/actions/workflows/reproduce-record.yml/badge.svg?branch=main)](https://github.com/braxtonmensah/twin-check/actions/workflows/reproduce-record.yml)
[![test](https://github.com/braxtonmensah/twin-check/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/braxtonmensah/twin-check/actions/workflows/test.yml)

![LRGB Peptides-func](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/braxtonmensah/twin-check/badges/badges/lrgb-peptides-func.json)
![AMPlify positives](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/braxtonmensah/twin-check/badges/badges/amplify-positives.json)
![FLIP meltome mixed](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/braxtonmensah/twin-check/badges/badges/flip-meltome-mixed.json)
![LRGB corrected split](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/braxtonmensah/twin-check/badges/badges/lrgb-corrected-v1.json)

A GitHub Action applying the same 0.80 rule as Silicyte's published leakage record; this repo's CI reproduces that record.

## What it measures

For a training split and a test split of sequences (peptides, proteins, any strings):

- **twin rate**: the percent of test rows whose best training match has a Levenshtein ratio of 0.80 or above.
  The ratio is `(len(a) + len(b) - indel(a, b)) / (len(a) + len(b))`, which is `rapidfuzz.fuzz.ratio / 100`
  and the python-Levenshtein ratio.
- **exact rate**: the percent of test rows whose sequence appears verbatim in the training split.

This is the rule behind [silicyte.com/data/leakage-measurements.json](https://silicyte.com/data/leakage-measurements.json)
and the page [silicyte.com/benchmark-leakage](https://silicyte.com/benchmark-leakage.html).

## Use it in one line

```yaml
      - uses: braxtonmensah/twin-check@v1
        with: { train: data/train.fasta, test: data/test.fasta }
```

A complete workflow is in [examples/workflow.yml](examples/workflow.yml). CSV works too:
`with: { train: train.csv, test: test.csv, format: csv, column: sequence }`.

| input | default | meaning |
|---|---|---|
| `train`, `test` | required | paths to the two splits |
| `format` | `fasta` | `fasta` or `csv` |
| `column` | | sequence column for CSV |
| `threshold` | `0.80` | Levenshtein ratio that counts as a twin |
| `fail-above` | empty | fail the step if the twin rate (percent) is above this; empty never fails |
| `badge-path` | `.twin-check/badge.json` | shields.io endpoint JSON |
| `record-path` | `.twin-check/record.json` | JSON record (leakage-measurements schema 1.0 fields, input MD5s, rule, tool version) |
| `summary-path` | `.twin-check/summary.md` | copy of the markdown table written to the job summary |
| `label` | `test twins @0.80` | badge label |

Outputs: `twin_rate`, `exact_rate` (percent, two decimals), `n_train`, `n_test`.

Inputs that would otherwise read as a silent zero stop the step with an error instead: an empty file, sequence
lines before any FASTA header, a header with no sequence, duplicate FASTA headers, duplicate CSV column
names, a missing CSV column, an empty sequence cell, or a header-only CSV. Wrapped FASTA lines are joined
and a leading UTF-8 byte order mark is ignored. Sequences are otherwise compared exactly as given: case
matters and nothing is normalized, so uppercase both splits first if your files mix case.

This repo's CI runs the action on `ubuntu-latest` only. The install step uses the venv's `bin/python`, so it
expects a Linux or macOS runner; Windows runners are not supported.

### Badge colors

For presentation, not a rating: below 1% `brightgreen`, 1% to below 10% `yellow`, 10% to below 30% `orange`,
30% and above `red`. To show your own badge, commit the badge JSON somewhere public and point
`https://img.shields.io/endpoint?url=<raw URL>` at it, as this repo does on its `badges` branch.

## The published record, rerun by CI

[`reproduce-record.yml`](.github/workflows/reproduce-record.yml) runs on every push, pull request, weekly, and on
demand. It downloads each benchmark's own public split files from the URLs the record's scripts use
([ci/sources.json](ci/sources.json)), checks every MD5, runs this action through `uses: ./`, and asserts the
values in the live record to two decimals:

| benchmark | split | record | n_test |
|---|---|---:|---:|
| LRGB Peptides-func | `splits_random_stratified_peptide.pickle` | 65.29% | 2,331 |
| AMPlify | positive class, `train_common` to `test_common` | 63.47% | 835 |
| FLIP meltome | `mixed` | 0.29% | 3,134 |
| LRGB Peptides-func, Silicyte corrected split v1 | single linkage at 0.80 | 0.00% | 2,331 |

The same job runs the action on LRGB with `fail-above: 1` and asserts that the step fails. The LRGB and FLIP MD5s
are the ones in the record's write-up (the LRGB pair also matches LRGB's own loader). The record does not list
MD5s for the AMPlify files; this repo pinned them on 2026-10-27, and both files last changed upstream on 2022-07-03.
This repo redistributes none of the benchmark files: CI downloads them from their own hosts at run time, under
their own licenses. The `badges` branch holds only the computed numbers and the input MD5s.

Other controls in [`test.yml`](.github/workflows/test.yml): test equal to train gives 100.00% twins and 100.00%
exact; seeded disjoint random sequences give 0.00%; one-substitution mutants of 20-mers (ratio 0.95) give
100.00% twins and 0.00% exact; a pair at ratio 0.79 is not counted and one at 0.80 is; an independent early-exit
search with a length prefilter agrees with the main path; and the 50-row example fixture reads 40.00%.

## Where the metric code comes from

`twin_rate` in [twin_check.py](twin_check.py) is a byte-for-byte copy of the function in
[`thirdparty/verify_more.py`](https://github.com/braxtonmensah/silicyte-design-loops/blob/50a6eb704730607abd02bf381c811c68ebd4d378/thirdparty/verify_more.py)
at commit `50a6eb7`, sha256 of the function source
`a01b9ee7f3b5c034bc6f69dc38c7a17c4d068cd21144b1574ebb8a89c56d210f`. CI downloads the original, re-hashes it, and
checks both copies against that value. The library is rapidfuzz, as in the record's scripts. Those scripts do not
pin a version; this action installs rapidfuzz 3.13.0 and numpy 2.2.6 with pip `--require-hashes`, and
`actions/setup-python` is pinned by commit SHA. The LRGB value in the record was originally produced by the sibling
script `verify_thirdparty.py` (rapidfuzz `extractOne`, same rule); the CI run shows this copy gives the same 65.29.

The search is all pairs: rapidfuzz `cdist` with a score cutoff, in blocks of 256 test rows. Exact matches go
through a hash set. A length prefilter appears in the early-exit reference path (not the main path), and it is lossless:
since `indel(a, b) >= |len(a) - len(b)|`, the ratio is at most `2 * min(len) / (len(a) + len(b))`.

## What it does not tell you

- **Sequence identity is not activity.** A near-twin can have a different label, and a distant sequence can
  behave the same.
- **A twin rate is not a score inflation.** Nothing is retrained, so the number says nothing about how much any
  reported accuracy is overstated. In Silicyte's AMPlify follow-up, 517 of the 530 test positives with a training
  twin had no opposite-label twin
  ([RESULT_AMPLIFY_CROSS_LABEL.md](https://github.com/braxtonmensah/silicyte-design-loops/blob/50a6eb704730607abd02bf381c811c68ebd4d378/thirdparty/RESULT_AMPLIFY_CROSS_LABEL.md));
  that describes what the twins are, and it does not say how much a model's score moves because of them.
- **0.80 is one threshold**, the one the record uses. It is an input.
- **Edit distance on strings is not alignment.** It does not see homology that an aligner such as MMseqs2 would.
- Cost grows with `n_train * n_test * length`. Peptide sets take seconds; FLIP meltome (3,134 full-length test
  proteins against 20,206 unique training proteins) is the slow case in this repo's own CI.

## Prior art

Checking and fixing train/test similarity is established work. Use these when you need more than this check:

- [DataSAIL](https://github.com/kalininalab/DataSAIL): similarity-aware data splitting.
- [GraphPart](https://github.com/graph-part/graph-part): homology partitioning for protein sequences.
- [CD-HIT-2D](https://github.com/weizhongli/cdhit): compares two sequence sets at an identity threshold.
- [hashFrag](https://github.com/de-Boer-Lab/hashFrag): homology-aware splits for DNA sequence datasets.
- [deepchecks](https://github.com/deepchecks/deepchecks): train/test validation checks for tabular, vision and NLP data.

What this repo adds is narrow: the record's own rule as a one-line Action with a badge, and a CI that reruns the
published numbers from the upstream files.

## License

MIT. See [LICENSE](LICENSE).

Built by Silicyte ([silicyte.com](https://silicyte.com)) with AI assistance.
