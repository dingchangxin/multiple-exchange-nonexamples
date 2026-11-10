"""Tests of the pasture encoding, matroid support and three-term checker."""
import json
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import random
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'content'))
from p7_checks import (P_ZERO, P_ONE, x_power, p_label, p_multiply, p_is_null,
                       check_support_matroid, check_three_term_plucker, check_p_matroid)


def example_scope():
    notebook = json.loads((ROOT/'content/rank4-pasture-p7.ipynb').read_text())
    scope = dict(P_ZERO=P_ZERO, P_ONE=P_ONE, x_power=x_power)
    cell = next(c for c in notebook['cells'] if c['id'] == 'p7-example')
    exec(''.join(cell['source']), scope)
    return scope


def independent_null(terms):
    """Use the forbidden scaling/permutation orbit, not normalization."""
    nonzero = [t for t in terms if t is not None]
    if len(nonzero) < 3:
        return len(nonzero) == 0 or (len(nonzero) == 2 and nonzero[0] == nonzero[1])
    forbidden = {tuple(sorted((k, (k+3) % 7, (k+3) % 7))) for k in range(7)}
    return tuple(sorted(nonzero)) not in forbidden


class TestP7(unittest.TestCase):
    def test_group_and_labels(self):
        self.assertEqual(p_label(P_ZERO), '0')
        self.assertEqual(p_label(P_ONE), '1')
        self.assertEqual(p_label(x_power(1)), 'x')
        self.assertEqual(p_label(x_power(-1)), 'x^6')
        for a in range(7):
            self.assertIsNone(p_multiply(P_ZERO, a))
            self.assertEqual(p_multiply(a, P_ONE), a)
            self.assertEqual(p_multiply(a, x_power(-a)), P_ONE)
            for b in range(7):
                self.assertEqual(p_multiply(a, b), (a+b) % 7)

    def test_normalized_rule(self):
        for m, n in product(range(7), repeat=2):
            self.assertEqual(p_is_null(m, n, P_ONE),
                             (m, n) not in {(3, 3), (0, 4), (4, 0)})

    def test_zero_patterns(self):
        self.assertTrue(p_is_null(P_ZERO, P_ZERO, P_ZERO))
        for a in range(7):
            self.assertFalse(p_is_null(a, P_ZERO, P_ZERO))
            for b in range(7):
                self.assertEqual(p_is_null(P_ZERO, a, b), a == b)

    def test_full_null_table_and_pasture_invariances(self):
        elements = (P_ZERO, *range(7))
        total = 0
        for terms in product(elements, repeat=3):
            expected = independent_null(terms)
            self.assertEqual(p_is_null(*terms), expected)
            total += expected
            for perm in permutations(terms):
                self.assertEqual(p_is_null(*perm), expected)
            for scalar in range(7):
                scaled = tuple(p_multiply(scalar, a) for a in terms)
                self.assertEqual(p_is_null(*scaled), expected)
        self.assertEqual(total, 344)

    def test_invalid_inputs(self):
        for value in (True, 1.0, 'x', complex(1)):
            with self.assertRaises(TypeError):
                p_is_null(value, P_ONE, P_ONE)
        for value in (-1, 7):
            with self.assertRaises(ValueError):
                check_p_matroid([1, 2], 1, lambda B: value)
        for value in (True, 1.0, '1'):
            with self.assertRaises(TypeError):
                x_power(value)
        for E, r in (([1, 1], 1), ([1], -1), ([1], 2), ([1], True)):
            with self.assertRaises(ValueError):
                check_p_matroid(E, r, lambda B: P_ONE)
        with self.assertRaises(TypeError):
            check_p_matroid([[1]], 1, lambda B: P_ONE)
        with self.assertRaises(TypeError):
            check_p_matroid([1], 1, 0)
        for mapping in ({(1, 2): P_ONE},
                        {(1, 1): P_ONE},
                        {(1, 2): P_ONE, (2, 1): P_ONE}):
            with self.assertRaises(ValueError):
                check_p_matroid([1, 2, 3], 2, mapping)

    def test_case_partition_and_coordinate_counts(self):
        s = example_scope()
        values = {}
        for B in combinations(s['E'], 4):
            S = frozenset(B)
            mixed = len(S & s['X']) == len(S & s['Y']) == 1
            cases = (S in s['zero_sets'], S in s['x_sets'], S in s['x2_sets'],
                     mixed and S & s['core'] in s['six_pairs'],
                     mixed and S & s['core'] in s['five_pairs'])
            self.assertLessEqual(sum(cases), 1)
            values[S] = s['chi'](S)
        self.assertEqual(Counter(values.values()),
                         Counter({None: 4, 0: 42, 1: 4, 2: 4, 5: 8, 6: 8}))
        self.assertEqual({B for B, v in values.items() if v is None}, s['zero_sets'])

    def test_example_support(self):
        s = example_scope()
        result = check_support_matroid(s['E'], 4, s['chi'])
        self.assertTrue(result['is_matroid'])
        self.assertEqual(result['bases'], 66)
        self.assertEqual(result['basis_exchange_checks'], 8712)
        self.assertEqual(result['violations'], 0)

    def test_example_relations(self):
        s = example_scope()
        result = check_three_term_plucker(s['E'], 4, s['chi'])
        self.assertTrue(result['holds'])
        self.assertEqual(result['relations_checked'], 420)
        self.assertEqual(result['violations'], 0)
        self.assertEqual({k: v['checked'] for k, v in result['by_nonzero_terms'].items()},
                         {0: 0, 1: 0, 2: 142, 3: 278})

    def test_combined_callable_mapping_and_single_evaluation(self):
        s = example_scope()
        calls = Counter()
        def chi(B):
            calls[B] += 1
            return s['chi'](B)
        result = check_p_matroid(s['E'], 4, chi)
        self.assertTrue(result['is_p_matroid'])
        self.assertEqual(len(calls), 70)
        self.assertEqual(set(calls.values()), {1})
        values = {B: s['chi'](B) for B in calls}
        self.assertEqual(check_p_matroid(s['E'], 4, values), result)

    def test_support_is_an_independent_requirement(self):
        E = tuple(range(6))
        bases = {frozenset({0, 1, 2}), frozenset({3, 4, 5})}
        chi = lambda B: P_ONE if B in bases else P_ZERO
        result = check_p_matroid(E, 3, chi)
        # Every three-term product is zero, yet basis exchange fails.
        self.assertTrue(result['plucker']['holds'])
        self.assertFalse(result['support']['is_matroid'])
        self.assertEqual(result['support']['violations'], 6)
        self.assertIsNotNone(result['support']['first_failure'])
        self.assertFalse(result['is_p_matroid'])
        zero = check_p_matroid(E, 3, lambda B: P_ZERO)
        self.assertTrue(zero['plucker']['holds'])
        self.assertFalse(zero['support']['nonempty'])
        self.assertFalse(zero['is_p_matroid'])

    def test_nonnull_three_unit_relation(self):
        E = tuple(range(1, 5))
        values = {frozenset(B): P_ONE for B in combinations(E, 2)}
        values[frozenset({1, 2})] = x_power(3)
        values[frozenset({1, 3})] = x_power(3)
        result = check_p_matroid(E, 2, values)
        self.assertTrue(result['support']['is_matroid'])
        self.assertFalse(result['is_p_matroid'])
        self.assertEqual(result['plucker']['violations'], 1)
        self.assertEqual(result['plucker']['first_failure']['terms'], (3, 3, 0))

    def test_nonnull_two_unit_relation(self):
        E = tuple(range(1, 5))
        values = {frozenset(B): P_ONE for B in combinations(E, 2)}
        values[frozenset({1, 2})] = P_ZERO
        values[frozenset({1, 4})] = x_power(1)
        result = check_p_matroid(E, 2, values)
        self.assertTrue(result['support']['is_matroid'])
        self.assertFalse(result['is_p_matroid'])
        self.assertEqual(result['plucker']['first_failure']['terms'], (None, 0, 1))

    def test_extreme_ranks_and_mixed_labels(self):
        for E, r in (([], 0), ([1, 2], 0), ([1, 2], 1), ([1, 2], 2)):
            result = check_p_matroid(E, r, lambda B: P_ONE)
            self.assertTrue(result['is_p_matroid'])
            self.assertEqual(result['plucker']['relations_checked'], 0)
        self.assertTrue(check_p_matroid(['a', ('b', 2), 3], 1,
                                        lambda B: P_ONE)['is_p_matroid'])

    def test_scaling_and_relabeling(self):
        s = example_scope()
        baseline = check_p_matroid(s['E'], 4, s['chi'])
        for scalar in range(7):
            result = check_p_matroid(s['E'], 4,
                                    lambda B: p_multiply(scalar, s['chi'](B)))
            self.assertEqual(result, baseline)
        labels = tuple(range(8))
        result = check_p_matroid(labels, 4,
                                lambda B: s['chi'](frozenset(s['E'][i] for i in B)))
        self.assertEqual(result, baseline)

    def test_independent_all_ordered_tuples(self):
        rng = random.Random(7)
        for E, r in ((tuple(range(4)), 2), (tuple(range(5)), 3)):
            for trial in range(5):
                values = {frozenset(B): (P_ONE if trial == 0 else
                                        rng.choice((P_ZERO, *range(7))))
                          for B in combinations(E, r)}
                def value(t):
                    return values[frozenset(t)] if len(set(t)) == r else P_ZERO
                failures = 0
                # Include repeated J and repeated a,b,c,d, not only distinct sets.
                for J in product(E, repeat=r-2):
                    for a,b,c,d in product(E, repeat=4):
                        terms = []
                        for i,j,k,l in ((a,b,c,d),(a,c,b,d),(a,d,b,c)):
                            v, w = value(J+(i,j)), value(J+(k,l))
                            terms.append(None if v is None or w is None else (v+w)%7)
                        failures += not independent_null(terms)
                result = check_three_term_plucker(E, r, values)
                self.assertEqual(result['holds'], failures == 0)
                self.assertEqual(failures, 24 * result['violations'])


if __name__ == '__main__':
    unittest.main()
