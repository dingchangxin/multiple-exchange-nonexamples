# Multiple-exchange nonexamples

Exact, reproducible computations over the triangular hyperfield T1.
All matroids considered here are **strong** matroids.

## Run in your browser

[Open the rank-4 example](https://dingchangxin.github.io/multiple-exchange-nonexamples/lab/index.html?path=rank4-me1-me2.ipynb)

Choose **Run → Run All Cells**. The site runs Python in your browser using
JupyterLite; no account or local Python installation is needed by readers.

## Contents

- [`content/t1_checks.py`](content/t1_checks.py): reusable exact functions
  `check_strong_t1`, `check_ME1`, `check_ME2`, and `exchange_witness`.
- [`content/rank4-me1-me2.ipynb`](content/rank4-me1-me2.ipynb): the rank-4
  example on eight elements with values 1, 3, 4, 7, all strong GP checks,
  exhaustive ME checks, and both exchange tables from the paper.
- [`tests/test_t1_checks.py`](tests/test_t1_checks.py): regression tests,
  including non-uniform inputs, exact boundaries, a weak-but-not-strong
  example, and independent small-case cross-checks.

## Reusable functions

Place new notebooks alongside `t1_checks.py` in `content/`:

```python
from t1_checks import check_strong_t1, check_ME1, check_ME2

# Define E, r and chi(B), with B passed as a frozenset.
strong = check_strong_t1(E, r, chi)
assert strong['is_strong'], strong
me1 = check_ME1(E, r, chi)
me2 = check_ME2(E, r, chi)
print(me1['holds'], me2['holds'])
```

The checkers enumerate all relevant subsets, not random samples. They return
report dictionaries containing verdicts, counts and a first failure when one
exists. Read `result['holds']`, not `bool(result)`.

Coordinates may be integers, `fractions.Fraction`, `decimal.Decimal`, or rational
strings such as `'1/3'`. Floating-point and symbolic irrational coordinates are
not supported. Zero coordinates are allowed; complete mappings must explicitly
include them. The functions use only the Python standard library.

### Exchange domain

Both ME checkers always use `I1 ⊆ B1 - B2`, as in the corrected paper, and
`I2 ⊆ B2 - B1` with `|I2| = |I1|`. There is no alternative-domain option.
This does **not** require `B1` and `B2` themselves to be disjoint.
Empty `I1` is included and has exactly one choice, `I2 = ∅`; every exchanged
set has size `r`. The helper `exchange_witness` raises `ValueError` when
`I1` is not a subset of `B1 - B2`.

## Verified results

For the supplied example, all 3,136 indexed strong GP relations pass:
560 with 2 terms, 1,680 with 3, 840 with 4, and 56 with 5.
Each property is checked on 22,470 ordered triples;
there are four violations of ME1 and four of ME2, all of exchange size 2.
The specific witnesses give `49 not in [0, 46]` and `9 not in [12, 86]`.

## Test locally

From the repository root, run:

```bash
python -m unittest discover -s tests -v
```

The deployment workflow runs these tests, including all notebook code cells,
before building and publishing the site. Notebook output cells record a local
CPython execution; running all cells recomputes them in the browser.

## Reference

Matthew Baker and Nathan Bowler, *Matroids over hyperfields*,
[arXiv:1601.01204v5](https://arxiv.org/abs/1601.01204v5), Definition 3.11,
Theorem 3.13, and Example 3.30 (used as a negative regression test).
