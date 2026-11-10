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
from math import comb

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


def independent_me_counts(E, r, values, which):
    """Independent ordered-tuple implementation, not using checker helpers."""
    applicable = failures = checked = 0
    bases = list(combinations(E, r))
    for B1 in bases:
        for B2 in bases:
            target = values[frozenset(B1)] * values[frozenset(B2)]
            for bits in product((0, 1), repeat=r):
                I1 = tuple(B1[i] for i in range(r) if bits[i])
                if set(I1) & set(B2):
                    continue
                terms = []
                for I2 in combinations([x for x in B2 if x not in B1], len(I1)):
                    left = tuple(x for x in B1 if x not in I1) + I2
                    right = tuple(x for x in B2 if x not in I2) + I1
                    assert len(set(left)) == len(set(right)) == r
                    term = values[frozenset(left)] * values[frozenset(right)]
                    terms.append(term)
                checked += 1
                m = max(terms)
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
            self.assertTrue(f(range(1, 5), 2, good)['holds'])
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
        self.assertEqual(result['relations_checked'], 2576)
        self.assertEqual(result['relations_skipped_trivial'], 560)
        self.assertEqual(result['relations_total'], 3136)
        self.assertEqual({k: v['checked'] for k, v in result['by_terms'].items()},
                         {3: 1680, 4: 840, 5: 56})

    def test_paper_cases_mutually_exclusive(self):
        ground = frozenset(E8)
        for B in map(frozenset, combinations(E8, 4)):
            cases = (B in (L, ground - L), B in (H, ground - H),
                     len(B & L) == 2 and len(B & H) in (1, 3))
            self.assertLessEqual(sum(cases), 1)
            matching = [v for applies, v in zip(cases, (1, 7, 4)) if applies]
            self.assertEqual(paper_chi(B), matching[0] if matching else 3)

    def test_all_automatic_gp_relations(self):
        # Rank 1 and corank 1: every indexed relation has J subset I.
        # Include zero coordinates and reject the all-zero function.
        for n, r in ((4, 1), (4, 3), (2, 1)):
            E = tuple(range(n))
            special = frozenset(E[:r])
            for nonzero in (False, True):
                values = {frozenset(B): int(nonzero and frozenset(B) == special)
                          for B in combinations(E, r)}
                result = check_strong_t1(E, r, values)
                self.assertEqual(result['is_strong'], nonzero)
                self.assertEqual(result['relations_checked'], 0)
                self.assertEqual(result['relations_skipped_trivial'],
                                 comb(n, r + 1) * comb(n, r - 1))
                self.assertEqual(result['relations_total'],
                                 result['relations_skipped_trivial'])
                self.assertEqual(result['by_terms'], {})
                self.assertEqual(result['violations'], 0)
        for E, r in (([], 0), ([1, 2], 0), ([1, 2], 2)):
            result = check_strong_t1(E, r, lambda B: 1)
            self.assertEqual(result['relations_total'], 0)
            self.assertEqual(result['relations_skipped_trivial'], 0)

    def test_skip_matches_unoptimized_gp_reference(self):
        # Independently evaluate EVERY pair, including the skipped identities.
        rng = random.Random(20261111)
        for n, r in ((4, 2), (5, 2), (5, 3), (6, 3)):
            E = tuple(range(n))
            for _ in range(10):
                values = {frozenset(B): Fraction(rng.randrange(8), 3)
                          for B in combinations(E, r)}
                stats = {}
                for I in map(frozenset, combinations(E, r + 1)):
                    for J in map(frozenset, combinations(E, r - 1)):
                        terms = [values[I - {e}] * values[J | {e}] for e in I - J]
                        slack = sum(terms) - 2 * max(terms)
                        k = len(terms)
                        entry = stats.setdefault(k, {'checked': 0, 'failures': 0,
                                                     'minimum_slack': None})
                        entry['checked'] += 1
                        entry['failures'] += int(slack < 0)
                        if entry['minimum_slack'] is None or slack < entry['minimum_slack']:
                            entry['minimum_slack'] = slack
                        if J <= I:
                            self.assertEqual(terms[0], terms[1])
                            self.assertEqual(slack, 0)
                result = check_strong_t1(E, r, values)
                failures = sum(v['failures'] for v in stats.values())
                self.assertEqual(result['is_strong'], any(values.values()) and failures == 0)
                self.assertEqual(result['violations'], failures)
                self.assertEqual(result['relations_total'], comb(n, r + 1) * comb(n, r - 1))
                self.assertEqual(result['relations_skipped_trivial'], stats[2]['checked'])
                self.assertEqual(result['by_terms'], {k: v for k, v in stats.items() if k != 2})
                self.assertEqual(result['relations_checked'],
                                 result['relations_total'] - result['relations_skipped_trivial'])

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
        expected = {1: (22470, 15478, 4), 2: (22470, 6992, 4)}
        for which, counts in expected.items():
            f = check_ME1 if which == 1 else check_ME2
            result = f(E8, 4, paper_chi)
            self.assertEqual(tuple(result[k] for k in
                ('triples_checked','hypothesis_satisfied','violations')), counts)
            self.assertFalse(result['holds'])
            self.assertEqual(result['violations_by_exchange_size'], {2: 4})
            w = result['first_failure']
            self.assertTrue(set(w['I1']) <= set(w['B1']) - set(w['B2']))

    def test_exchange_domain_and_equality(self):
        one = lambda B: 1
        for f in (check_ME1, check_ME2):
            result = f([1, 2, 3], 2, one)
            self.assertTrue(result['holds'])
            self.assertEqual(result['triples_checked'], 15)
        self.assertEqual(check_ME1([1, 2, 3], 2, one)['hypothesis_satisfied'], 15)
        self.assertEqual(check_ME2([1, 2, 3], 2, one)['hypothesis_satisfied'], 0)
        w = exchange_witness([1, 2], 1, one, [1], [2], [1])
        self.assertTrue(w['ME1_hypothesis'])
        self.assertFalse(w['ME2_hypothesis'])
        self.assertTrue(w['membership'])

    def test_witness_rejects_forbidden_exchange_elements(self):
        for B1, B2, I1 in (([1, 2], [2, 3], [2]),
                           ([1, 2], [2, 3], [1, 2]),
                           ([1, 2], [1, 2], [1])):
            with self.assertRaisesRegex(ValueError, 'subset of B1 - B2'):
                exchange_witness([1, 2, 3], 2, lambda B: 1, B1, B2, I1)

    def test_overlapping_bases_and_empty_exchange(self):
        w = exchange_witness([1, 2, 3], 2, lambda B: 1, [1, 2], [2, 3], [1])
        self.assertEqual(len(w['rows']), 1)
        self.assertEqual(w['rows'][0]['I2'], (3,))
        self.assertEqual(w['rows'][0]['left'], (2, 3))
        self.assertEqual(w['rows'][0]['right'], (1, 2))
        self.assertTrue(w['membership'])
        for E, r, B in (([], 0, []), ([1, 2], 0, []), ([1, 2], 2, [1, 2])):
            for f in (check_ME1, check_ME2):
                result = f(E, r, lambda S: 1)
                self.assertTrue(result['holds'])
                self.assertEqual(result['triples_checked'], 1)
            w = exchange_witness(E, r, lambda S: 1, B, B, [])
            self.assertEqual(w['products'], (1,))
            self.assertEqual(w['rows'][0]['I2'], ())
            self.assertTrue(w['membership'])

    def test_independent_random_crosscheck(self):
        rng = random.Random(20261110)
        E = tuple(range(4))
        for _ in range(20):
            values = {frozenset(B): rng.randrange(5) for B in combinations(E, 2)}
            for which, f in ((1, check_ME1), (2, check_ME2)):
                result = f(E, 2, values)
                expected = independent_me_counts(E, 2, values, which)
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
        self.assertEqual(check_ME2(E8, 4, scaled)['violations'], 4)
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
