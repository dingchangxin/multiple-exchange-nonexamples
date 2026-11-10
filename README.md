# Multiple-exchange nonexamples

Exact, reproducible computations over the triangular hyperfield T1 and a pasture P.
The T1 examples check **strong** matroids. The pasture example checks matroidal
support and three-term Plucker relations in the sense of Baker--Lorscheid; it
does not claim higher-term relations or test multiple-exchange properties.

## Run in your browser

[Open the rank-4 triangular-hyperfield example (ME1 and ME2)](https://dingchangxin.github.io/multiple-exchange-nonexamples/lab/index.html?path=rank4-me1-me2.ipynb)

[Open the rank-6 triangular-hyperfield example (ME3)](https://dingchangxin.github.io/multiple-exchange-nonexamples/lab/index.html?path=rank6-me3.ipynb)

[Open the pasture example](https://dingchangxin.github.io/multiple-exchange-nonexamples/lab/index.html?path=rank4-pasture-p7.ipynb)

Choose **Run → Run All Cells**. The site runs Python in your browser using
JupyterLite; no account or local Python installation is needed by readers.

## Contents

- [`content/t1_checks.py`](content/t1_checks.py): reusable exact functions
  `check_strong_t1`, `check_ME1`, `check_ME2`, `check_ME3`, and `exchange_witness`.
- [`content/rank4-me1-me2.ipynb`](content/rank4-me1-me2.ipynb): the rank-4
  example on eight elements with values 1, 3, 4, 7, all strong GP checks,
  exhaustive ME checks, and both exchange tables from the paper.
- [`tests/test_t1_checks.py`](tests/test_t1_checks.py): regression tests,
  including non-uniform inputs, exact boundaries, a weak-but-not-strong
  example, and independent small-case cross-checks.

- [`content/rank6-me3.ipynb`](content/rank6-me3.ipynb): a self-contained
  computational supplement presenting the full ME3 counterexample referred
  to in the paper, including all auxiliary values and the twenty-term witness.
- [`content/rank6_me3_data.py`](content/rank6_me3_data.py): all 924 exact
  coordinates, adapted without numerical changes from counterexampleII.py
  and cross-checked against counterexampleII.tex.
- [`tests/test_me3.py`](tests/test_me3.py): rank-six regression tests and
  independent checks of the ME3 enumeration and exact-integer arithmetic.

- [`content/p7_checks.py`](content/p7_checks.py): exact multiplication and null
  triples for the specified pasture P, and the reusable `check_p_matroid`,
  `check_support_matroid`, and `check_three_term_plucker` functions.
- [`content/rank4-pasture-p7.ipynb`](content/rank4-pasture-p7.ipynb): the new
  rank-4 example on A, B, C, D, X1, X2, Y1, Y2; support and three-term checks only.
- [`tests/test_p7_checks.py`](tests/test_p7_checks.py): pasture regression tests,
  including exhaustive null-table invariance, zero terms, nonmatroidal support,
  and independent ordered-tuple cross-checks.

## Triangular-hyperfield functions

Place new notebooks alongside `t1_checks.py` in `content/`:

```python
from t1_checks import check_strong_t1, check_ME1, check_ME2, check_ME3

# Define E, r and chi(B), with B passed as a frozenset.
strong = check_strong_t1(E, r, chi)
assert strong['is_strong'], strong
me1 = check_ME1(E, r, chi)
me2 = check_ME2(E, r, chi)
me3 = check_ME3(E, r, chi)
print(me1['holds'], me2['holds'], me3['holds'])
```

The checkers enumerate all relevant subsets, not random samples. They return
report dictionaries containing verdicts, counts and a first failure when one
exists. Read `result['holds']`, not `bool(result)`.

The strong checker skips the automatic relations with `J ⊆ I`: their two
products are identical, including when zero. All other strong GP relations
are explicitly evaluated. The report separates `relations_total`,
`relations_skipped_trivial`, and `relations_checked`; `by_terms` contains only
evaluated relations, not a two-term row.

Coordinates may be integers, `fractions.Fraction`, `decimal.Decimal`, or rational
strings such as `'1/3'`. Floating-point and symbolic irrational coordinates are
not supported. Zero coordinates are allowed; complete mappings must explicitly
include them. The functions use only the Python standard library.

For larger rational examples, `check_strong_t1(..., arithmetic="integer")`
clears denominators exactly. Multiplying every coordinate by a common positive
integer D multiplies every product and slack by D squared; the verdict is
unchanged. Returned slacks and failure witnesses use the original scale.
The default `arithmetic="fraction"` retains the original exact arithmetic.

### Exchange domain

All ME checkers always use `I1 ⊆ B1 - B2`, as in the corrected paper, and
`I2 ⊆ B2 - B1` with `|I2| = |I1|`. There is no alternative-domain option.
This does **not** require `B1` and `B2` themselves to be disjoint.
Empty `I1` is included and has exactly one choice, `I2 = ∅`; every exchanged
set has size `r`. The helper `exchange_witness` raises `ValueError` when
`I1` is not a subset of `B1 - B2`.

ME3 tests only triples with zero target product. Pairs with two nonzero
coordinates are omitted before forming exchange products because their
hypothesis is false. Its report distinguishes `triples_total`,
`triples_checked`, and `triples_skipped_nonzero`. The helper also reports
`ME3_hypothesis` and `violates_ME3` for a specified exchange.

## Triangular-hyperfield results

### Rank four: ME1 and ME2

The first three cases defining the example are mutually exclusive; their
order does not matter, with the final case remaining the default.

For the supplied example, all 3,136 indexed strong GP relations are satisfied:
560 automatic two-term relations are skipped, and the remaining 2,576 are
explicitly checked (1,680 with 3 terms, 840 with 4, and 56 with 5).
Each property is checked on 22,470 ordered triples;
there are four violations of ME1 and four of ME2, all of exchange size 2.
The specific witnesses give `49 not in [0, 46]` and `9 not in [12, 86]`.

### Rank six: ME3

The new example has 924 coordinates and exactly one nonbasis, `{7,8,9,10,11,12}`.
All 627,264 indexed strong GP relations hold: 16,632 automatic relations are
skipped and 610,632 are explicitly evaluated, including seven-term relations.
The ME3 checker evaluates all 17,977 triples with zero target and finds two
violations, both of exchange size three. It omits 8,287,859 vacuous triples
with nonzero target. The displayed twenty exchanged products have largest
term 1 and sum of the other nineteen 77/80. Their hypersum is `[3/80, 157/80]`,
which does not contain the target 0.

The complete coordinate definition is in the notebook, so the paper can
refer readers directly to it without reproducing the example's large table.

## Pasture functions and results

For the specific pasture with unit group Z/7Z and -1=1, use:

```python
from p7_checks import P_ZERO, P_ONE, x_power, check_p_matroid

# Define E, r and chi(B) for this pasture.
# P_ZERO=None is pasture zero; P_ONE=0 is the identity; x_power(k) encodes x^k.
report = check_p_matroid(E, r, chi)
print(report['is_p_matroid'])
print(report['support'])
print(report['plucker'])
```

The integer 0 encodes x^0=1, **not** pasture zero. All nonzero values are encoded
by exponents modulo 7. For three units, normalize by the third and exclude the
pairs (3,3), (0,4), (4,0). With zero terms, 0+a+b is null exactly when a=b.

The combined checker separately verifies nonempty support and basis exchange,
and then all three-term relations, including zero terms. The standalone
`check_three_term_plucker` does not imply that the support is a matroid.

The pasture example has 70 coordinates and 66 bases. All 8,712 basis-exchange
conditions and all 420 three-term relations pass. Among the latter, 142 have
two nonzero products and 278 have three. No higher-term Plucker relations or
multiple-exchange properties are tested for this example.

## Test locally

From the repository root, run:

```bash
python -m unittest discover -s tests -v
```

The deployment workflow runs these tests, including all notebook code cells,
before building and publishing the site. Notebook output cells record a local
CPython execution; running all cells recomputes them in the browser.

## References

Matthew Baker and Nathan Bowler, *Matroids over hyperfields*,
[arXiv:1601.01204v5](https://arxiv.org/abs/1601.01204v5), Definition 3.11,
Theorem 3.13, and Example 3.30 (used as a negative regression test).

Matthew Baker and Oliver Lorscheid, *Foundations of matroids, Part 1: Matroids
without large uniform minors*, [arXiv:2008.00014v1](https://arxiv.org/abs/2008.00014v1),
p. 5 for pasture representations and Definition 2.1 for pastures.
