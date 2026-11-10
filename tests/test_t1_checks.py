"""Regression tests; run with python -m unittest discover -s tests -v."""
import contextlib
import io
import json
from pathlib import Path
import random
import sys
import unittest
from decimal import Decimal
from fractions import Fraction
from itertools import combinations, product

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'content'))
from t1_checks import (check_strong_t1, check_ME1, check_ME2,
                       exchange_witness, t1_hypersum_interval, t1_contains)

E8 = tuple(range(1, 9))
L = frozenset({1, 2, 5, 6})
H = frozenset({3, 4, 5, 6})


def paper_chi(B):
    if B in (L, frozenset(E8) - L):
        return 1
    if B in (H, frozenset(E8) - H):
        return 7
    if len(B & L) == 2 and len(B & H) in (1, 3):
        return 4
    return 3


def independent_me_counts(E, r, values, which, restricted):
    """Independent ordered-tuple implementation, not using checker helpers."""
    applicable = failures = checked = 0
    bases = list(combinations(E, r))
    for B1 in bases:
        for B2 in bases:
            target = values[frozenset(B1)] * values[frozenset(B2)]
            for bits in product((0, 1), repeat=r):
                I1 = tuple(B1[i] for i in range(r) if bits[i])
                if restricted and set(I1) & set(B2):
                    continue
                terms = []
                for I2 in combinations([x for x in B2 if x not in B1], len(I1)):
                    left = tuple(x for x in B1 if x not in I1) + I2
                    right = tuple(x for x in B2 if x not in I2) + I1
                    term = 0 if len(set(right)) != r else (
                        values[frozenset(left)] * values[frozenset(right)])
                    terms.append(term)
                checked += 1
                m = max(terms, default=0)
                if (target >= m) == (which == 1):
                    applicable += 1
                    polygon = [target] + terms
                    if any(x > sum(polygon) - x for x in polygon):
                        failures += 1
    return checked, applicable, failures


class TestT1(unittest.TestCase):
    def test_interval_and_boundaries(self):
        self.assertEqual(t1_hypersum_interval([]), (0, 0))
        self.assertEqual(t1_hypersum_interval([3]), (3, 3))
        self.assertEqual(t1_hypersum_interval([0, 0]), (0, 0))
        self.assertEqual(t1_hypersum_interval([49, 9, 9, 9, 9, 1]), (12, 86))
        self.assertTrue(t1_contains(0, [1, 1]))
        self.assertTrue(t1_contains(3, [1, 2]))
        self.assertFalse(t1_contains(0, [1, 2]))
        self.assertEqual(t1_hypersum_interval(['1/3', Decimal('0.5')]),
                         (Fraction(1, 6), Fraction(5, 6)))

    def test_invalid_inputs(self):
        for x in (0.1, float('nan'), float('inf'), 1j):
            with self.assertRaises(TypeError):
                check_strong_t1([1, 2], 1, lambda B: x)
        for x in (-1, 'NaN', Decimal('Infinity')):
            with self.assertRaises(ValueError):
                check_strong_t1([1, 2], 1, lambda B: x)
        for E, r in (([1, 1], 1), ([1], 2), ([1], -1), ([1], True)):
            with self.assertRaises(ValueError):
                check_strong_t1(E, r, lambda B: 1)
        with self.assertRaises(ValueError):
            check_strong_t1([1, 2, 3], 2, {(1, 2): 1})
        with self.assertRaises(ValueError):
            exchange_witness([1, 2, 3], 2, lambda B: 1, [1, 2], [2, 3], [3])

    def test_zero_and_extreme_ranks(self):
        self.assertFalse(check_strong_t1([1, 2], 1, lambda B: 0)['is_strong'])
        for E, r in (([], 0), ([1, 2], 0), ([1, 2], 2)):
            self.assertTrue(check_strong_t1(E, r, lambda B: 1)['is_strong'])
        for f in (check_ME1, check_ME2):
            self.assertTrue(f([], 0, lambda B: 1)['holds'])

    def test_nonuniform_and_bad_support(self):
        # Absolute 2x2 minors of columns (1,0), (0,1), (1,1), (0,0).
        good = {frozenset(B): int(4 not in B) for B in combinations(range(1, 5), 2)}
        self.assertTrue(check_strong_t1(range(1, 5), 2, good)['is_strong'])
        for f in (check_ME1, check_ME2):
            self.assertTrue(f(range(1, 5), 2, good, disjoint_only=True)['holds'])
        bad = lambda B: int(B in (frozenset({1, 2}), frozenset({3, 4})))
        self.assertFalse(check_strong_t1(range(1, 5), 2, bad)['is_strong'])

    def test_weak_is_not_strong(self):
        # Baker--Bowler arXiv:1601.01204v5, Example 3.30.
        def weak(B):
            if B == frozenset({1, 5, 6}):
                return 4
            if 1 in B and len(B & {2, 3, 4}) == 1 and len(B & {5, 6}) == 1:
                return 2
            return 1
        result = check_strong_t1(range(1, 7), 3, weak)
        self.assertFalse(result['is_strong'])
        self.assertEqual(result['by_terms'][3]['failures'], 0)
        self.assertGreater(result['by_terms'][4]['failures'], 0)

    def test_paper_strong(self):
        result = check_strong_t1(E8, 4, paper_chi)
        self.assertTrue(result['is_strong'])
        self.assertEqual(result['relations_checked'], 3136)
        self.assertEqual({k: v['checked'] for k, v in result['by_terms'].items()},
                         {2: 560, 3: 1680, 4: 840, 5: 56})

    def test_paper_witnesses(self):
        w1 = exchange_witness(E8, 4, paper_chi, H, frozenset(E8)-H, {5, 6})
        w2 = exchange_witness(E8, 4, paper_chi, {1,2,3,4}, {5,6,7,8}, {1,2})
        self.assertEqual(w1['products'], (9, 9, 9, 9, 9, 1))
        self.assertEqual(w2['products'], (49, 9, 9, 9, 9, 1))
        self.assertEqual(w1['hypersum_interval'], (0, 46))
        self.assertEqual(w2['hypersum_interval'], (12, 86))
        self.assertTrue(w1['violates_ME1'])
        self.assertTrue(w2['violates_ME2'])

    def test_paper_exhaustive_counts(self):
        expected = {(False,1):(78400,71408,55934), (False,2):(78400,6992,4),
                    (True,1):(22470,15478,4), (True,2):(22470,6992,4)}
        for (restricted, which), counts in expected.items():
            f = check_ME1 if which == 1 else check_ME2
            result = f(E8, 4, paper_chi, disjoint_only=restricted)
            self.assertEqual(tuple(result[k] for k in
                ('triples_checked','hypothesis_satisfied','violations')), counts)
            self.assertFalse(result['holds'])

    def test_literal_domain_and_equality(self):
        one = lambda B: 1
        literal = check_ME1([1, 2, 3], 2, one)
        self.assertFalse(literal['holds'])
        self.assertEqual(literal['first_failure']['products'], ())
        self.assertTrue(check_ME1([1, 2, 3], 2, one, disjoint_only=True)['holds'])
        w = exchange_witness([1, 2], 1, one, [1], [2], [1])
        self.assertTrue(w['ME1_hypothesis'])
        self.assertFalse(w['ME2_hypothesis'])
        self.assertTrue(w['membership'])

    def test_independent_random_crosscheck(self):
        rng = random.Random(20261110)
        E = tuple(range(4))
        for _ in range(20):
            values = {frozenset(B): rng.randrange(5) for B in combinations(E, 2)}
            for restricted in (False, True):
                for which, f in ((1, check_ME1), (2, check_ME2)):
                    result = f(E, 2, values, disjoint_only=restricted)
                    expected = independent_me_counts(E, 2, values, which, restricted)
                    self.assertEqual(tuple(result[k] for k in
                        ('triples_checked','hypothesis_satisfied','violations')), expected)
            # Cross-check ALL ordered tuples, including repeated indices.
            def value(t):
                return values[frozenset(t)] if len(set(t)) == 2 else 0
            passes = any(values.values())
            for I in product(E, repeat=3):
                for j in E:
                    terms = [value(I[:i]+I[i+1:]) * value((I[i], j)) for i in range(3)]
                    passes = passes and all(t <= sum(terms)-t for t in terms)
            self.assertEqual(check_strong_t1(E, 2, values)['is_strong'], passes)

    def test_rational_scaling_and_labels(self):
        scaled = lambda B: Fraction(paper_chi(B), 3)
        self.assertTrue(check_strong_t1(E8, 4, scaled)['is_strong'])
        self.assertEqual(check_ME2(E8, 4, scaled, disjoint_only=True)['violations'], 4)
        self.assertTrue(check_strong_t1(['a', ('b', 2)], 1, lambda B: '1/7')['is_strong'])

    def test_notebooks_execute(self):
        for path in sorted((ROOT/'content').glob('*.ipynb')):
            notebook = json.loads(path.read_text())
            scope = {'__name__': '__main__'}
            with contextlib.redirect_stdout(io.StringIO()):
                for cell in notebook['cells']:
                    if cell['cell_type'] == 'code':
                        source = cell['source']
                        exec(compile(''.join(source), str(path), 'exec'), scope)



if __name__ == '__main__':
    unittest.main()
