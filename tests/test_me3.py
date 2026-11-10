"""ME3 and rank-six regression tests; all calculations are exact."""
import hashlib
import random
import sys
import unittest
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'content'))
from t1_checks import check_strong_t1, check_ME3, exchange_witness
from rank6_me3_data import E, r, G, A, B, eps, extra, chi, coordinates


def independent_me3(E, r, values):
    """Literal reference enumeration; no optimized checker helpers."""
    total = checked = violations = 0
    subsets = list(combinations(E, r))
    for B1 in map(set, subsets):
        for B2 in map(set, subsets):
            pool = [e for e in E if e in B1 and e not in B2]
            other = [e for e in E if e in B2 and e not in B1]
            for k in range(len(pool) + 1):
                for I1 in map(set, combinations(pool, k)):
                    total += 1
                    if values[frozenset(B1)] * values[frozenset(B2)] != 0:
                        continue
                    checked += 1
                    terms = [values[frozenset((B1-I1) | set(I2))]
                             * values[frozenset((B2-set(I2)) | I1)]
                             for I2 in combinations(other, k)]
                    violations += any(t > sum(terms)-t for t in terms)
    return total, checked, violations


class TestME3(unittest.TestCase):
    def test_integer_mode_matches_fraction_mode(self):
        rng = random.Random(20261110)
        for n, rank in ((0, 0), (4, 0), (4, 1), (4, 2), (5, 3), (6, 3), (4, 4)):
            ground = tuple(range(n))
            for _ in range(8):
                values = {frozenset(S): Fraction(rng.randrange(9), rng.choice((2, 3, 5, 7, 11)))
                          for S in combinations(ground, rank)}
                self.assertEqual(
                    check_strong_t1(ground, rank, values, arithmetic='fraction'),
                    check_strong_t1(ground, rank, values, arithmetic='integer'))
        with self.assertRaises(ValueError):
            check_strong_t1([1, 2], 1, lambda S: 1, arithmetic='float')

    def test_integer_failure_report_original_units(self):
        values = {frozenset(S): Fraction(1, 3) for S in combinations(range(4), 2)}
        values[frozenset((0, 1))] = Fraction(10, 7)
        a = check_strong_t1(range(4), 2, values)
        b = check_strong_t1(range(4), 2, values, arithmetic='integer')
        self.assertFalse(b['is_strong'])
        self.assertEqual(a, b)
        self.assertTrue(any(t.denominator != 1 for t in b['first_failure']['terms']))

    def test_me3_reference_random_inputs(self):
        rng = random.Random(20261112)
        for n, rank in ((4, 2), (5, 2), (5, 3)):
            ground = tuple(range(n))
            for _ in range(10):
                values = {frozenset(S): Fraction(rng.randrange(4), 3)
                          for S in combinations(ground, rank)}
                result = check_ME3(ground, rank, values)
                expected = independent_me3(ground, rank, values)
                self.assertEqual(tuple(result[k] for k in
                    ('triples_total', 'triples_checked', 'violations')), expected)
                self.assertEqual(result['triples_skipped_nonzero'], expected[0]-expected[1])
                self.assertEqual(result['holds'], expected[2] == 0)
                if result['first_failure']:
                    w = result['first_failure']
                    self.assertEqual(w['target'], 0)
                    self.assertTrue(w['violates_ME3'])
                    self.assertTrue(set(w['I1']) <= set(w['B1']) - set(w['B2']))

    def test_vacuity_zero_and_boundary_ranks(self):
        for ground, rank in (((), 0), ((1, 2), 0), ((1, 2), 2), (tuple(range(4)), 2)):
            for value in (0, Fraction(2, 7)):
                result = check_ME3(ground, rank, lambda S: value)
                self.assertTrue(result['holds'])
                self.assertEqual(result['violations'], 0)
                self.assertEqual(result['triples_checked'], result['triples_total'] if value == 0 else 0)
        # Nonuniform realizable rank-two data: a loop at 4.
        values = {frozenset(S): int(4 not in S) for S in combinations(range(1, 5), 2)}
        self.assertTrue(check_ME3(range(1, 5), 2, values)['holds'])
        with self.assertRaises(ValueError):
            check_ME3([1, 2, 3], 2, {(1, 2): 1})


class TestRankSixData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.strong = check_strong_t1(E, r, coordinates, arithmetic='integer')
        cls.me3 = check_ME3(E, r, coordinates)

    def test_exact_dataset_and_domain(self):
        self.assertEqual(len(extra), 155)
        self.assertEqual(set(coordinates), {frozenset(S) for S in combinations(E, r)})
        self.assertEqual([S for S, v in coordinates.items() if v == 0], [A | B])
        self.assertTrue(all(isinstance(v, Fraction) and v >= 0 for v in coordinates.values()))
        self.assertEqual(coordinates[G], Fraction(1511, 10000))
        self.assertEqual(eps, Fraction(1, 40))
        canonical = '\n'.join(','.join(map(str, sorted(S))) + ':' + str(v)
                              for S, v in coordinates.items())
        # Fingerprint matched against all 924 coordinates of counterexampleII.py.
        self.assertEqual(hashlib.sha256(canonical.encode()).hexdigest(),
            '12c2c31dcd05f7b3a66b835dcb15d16d4ee09273309e22ab3086e3a82b32321d')
        for S in ((1, 2, 3), (1, 2, 3, 4, 5, 13), (1, 1, 2, 3, 4, 5)):
            with self.assertRaises(ValueError):
                chi(S)

    def test_all_strong_relations(self):
        result = self.strong
        self.assertTrue(result['is_strong'])
        self.assertEqual(result['support_size'], 923)
        self.assertEqual(result['relations_total'], 627264)
        self.assertEqual(result['relations_skipped_trivial'], 16632)
        self.assertEqual(result['relations_checked'], 610632)
        self.assertEqual(result['violations'], 0)
        self.assertEqual({k: v['checked'] for k, v in result['by_terms'].items()},
                         {3: 138600, 4: 277200, 5: 166320, 6: 27720, 7: 792})
        self.assertEqual(result['by_terms'][7]['minimum_slack'],
                         Fraction(19865607066213170573, 312500000000000000000))

    def test_explicit_twenty_term_witnesses(self):
        for I1 in (A, B):
            w = exchange_witness(E, r, coordinates, A | B, G, I1)
            self.assertEqual(w['target'], 0)
            self.assertTrue(w['ME3_hypothesis'] and w['violates_ME3'])
            self.assertEqual(len(w['products']), 20)
            self.assertEqual(w['products'].count(1), 1)
            self.assertEqual(sum(w['products']), Fraction(157, 80))
            self.assertEqual(w['hypersum_interval'], (Fraction(3, 80), Fraction(157, 80)))
            self.assertEqual(sum(w['products']) - 2*max(w['products']), Fraction(-3, 80))

    def test_exhaustive_me3(self):
        result = self.me3
        self.assertFalse(result['holds'])
        self.assertEqual(result['pairs_with_zero_target'], 1847)
        self.assertEqual(result['triples_total'], 8305836)
        self.assertEqual(result['triples_checked'], 17977)
        self.assertEqual(result['triples_skipped_nonzero'], 8287859)
        self.assertEqual(result['violations'], 2)
        self.assertEqual(result['violations_by_exchange_size'], {3: 2})

    def test_notebook_displays_complete_auxiliary_table(self):
        import json, re
        notebook = json.loads((ROOT/'content/rank6-me3.ipynb').read_text())
        cell = next(c for c in notebook['cells'] if c['id'] == 'me3-auxiliary')
        table = ''.join(cell['source'])
        parsed = {}
        for number, keys in re.findall(r'\| `([^`]+)` \| (.*?) \|', table):
            for a, c, H in re.findall(r'v_\{(\d),(\d),([0-9]+)\}', keys):
                key = (int(a), int(c), tuple(map(int, H)))
                self.assertNotIn(key, parsed)
                parsed[key] = Fraction(number)
        self.assertEqual(parsed, extra)


if __name__ == '__main__':
    unittest.main()
