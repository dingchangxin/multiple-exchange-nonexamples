"""Exact checks for the specific pasture P with unit group Z/7Z and -1=1.

Encoding: P_ZERO=None; an integer e in 0,...,6 represents x^e.
In particular, the integer 0 is P_ONE, NOT the zero of the pasture.
Use x_power(k), P_ONE and P_ZERO when defining coordinates.

For three units, divide by the third to apply the supplied normalized rule.
For triples containing zero, use the pasture axiom 0+a+b null iff a=b.
Reference: Baker--Lorscheid, Foundations of matroids, Part 1,
arXiv:2008.00014v1, p. 5 (representations), Definition 2.1 (pastures).

The combined checker verifies nonempty matroidal support and all three-term
Plucker relations. It does not assert higher-term/strong relations or check
any multiple-exchange property. Only the Python standard library is needed.
"""
from collections.abc import Mapping
from itertools import combinations
from typing import Optional

PElement = Optional[int]
P_ZERO = None
P_ONE = 0
FORBIDDEN_PAIRS = frozenset({(3, 3), (0, 4), (4, 0)})

__all__ = ["P_ZERO", "P_ONE", "x_power", "p_label", "p_multiply", "p_is_null",
           "check_support_matroid", "check_three_term_plucker", "check_p_matroid"]


def _validate(value: PElement) -> PElement:
    if value is None:
        return value
    if type(value) is not int:
        raise TypeError("Use P_ZERO=None or integer exponents 0,...,6; "
                        "the integer 0 represents P_ONE.")
    if not 0 <= value < 7:
        raise ValueError("Exponents must be in 0,...,6; use x_power(k) to reduce k.")
    return value


def x_power(k: int) -> int:
    """Encode x^k, reducing its integer exponent modulo 7."""
    if type(k) is not int:
        raise TypeError("The exponent must be an integer (not a bool or float).")
    return k % 7


def p_label(a: PElement) -> str:
    """Human-readable label, distinguishing pasture zero from exponent zero."""
    _validate(a)
    if a is None:
        return "0"
    if a == 0:
        return "1"
    return "x" if a == 1 else "x^%s" % a


def _mul(a, b):
    return P_ZERO if a is None or b is None else (a + b) % 7


def p_multiply(a: PElement, b: PElement) -> PElement:
    """Multiply two encoded pasture elements; zero is absorbing."""
    return _mul(_validate(a), _validate(b))


def _null(a, b, c):
    nonzero = tuple(t for t in (a, b, c) if t is not None)
    if len(nonzero) == 0:
        return True
    if len(nonzero) == 1:
        return False
    if len(nonzero) == 2:
        return nonzero[0] == nonzero[1]
    return ((a - c) % 7, (b - c) % 7) not in FORBIDDEN_PAIRS


def p_is_null(a: PElement, b: PElement, c: PElement) -> bool:
    """Test a+b+c in N_P, including every possible zero pattern."""
    return _null(_validate(a), _validate(b), _validate(c))


class _Data:
    """Validate input and evaluate chi once on each r-subset."""
    def __init__(self, E, r, chi):
        self.E = tuple(E)
        try:
            ground = frozenset(self.E)
        except TypeError as exc:
            raise TypeError("Ground-set labels must be hashable.") from exc
        if len(ground) != len(self.E):
            raise ValueError("Ground-set labels must be distinct.")
        if type(r) is not int or not 0 <= r <= len(self.E):
            raise ValueError("r must be an integer between 0 and len(E).")
        self.r = r
        self.subsets = tuple(frozenset(B) for B in combinations(self.E, r))
        if isinstance(chi, Mapping):
            normalized = {}
            for key, value in chi.items():
                labels = tuple(key)
                B = frozenset(labels)
                if len(labels) != r or len(B) != r or not B <= ground:
                    raise ValueError("Every mapping key must be an r-subset of E.")
                if B in normalized:
                    raise ValueError("Duplicate subset keys after normalization.")
                normalized[B] = value
            if set(normalized) != set(self.subsets):
                raise ValueError("Supply every r-subset, explicitly including zeros.")
            self.values = {B: _validate(normalized[B]) for B in self.subsets}
        elif callable(chi):
            self.values = {B: _validate(chi(B)) for B in self.subsets}
        else:
            raise TypeError("chi must be callable or a complete mapping.")
        self.bases = tuple(B for B in self.subsets if self.values[B] is not None)

    def ordered(self, S):
        return tuple(e for e in self.E if e in S)


def _support(data):
    bases = set(data.bases)
    checked = failures = 0
    first = None
    for B1 in data.bases:
        for B2 in data.bases:
            candidates = data.ordered(B2 - B1)
            for e in data.ordered(B1 - B2):
                checked += 1
                if not any((B1 - {e}) | {f} in bases for f in candidates):
                    failures += 1
                    if first is None:
                        first = {"B1": data.ordered(B1), "B2": data.ordered(B2),
                                 "e": e, "candidate_f": candidates}
    return {"is_matroid": bool(bases) and failures == 0,
            "nonempty": bool(bases), "bases": len(bases),
            "basis_exchange_checks": checked, "violations": failures,
            "first_failure": first}


def _plucker(data):
    checked = failures = 0
    first = None
    stats = {k: {"checked": 0, "failures": 0} for k in range(4)}
    if data.r >= 2:
        for common in combinations(data.E, data.r - 2):
            J = frozenset(common)
            rest = tuple(e for e in data.E if e not in J)
            for a, b, c, d in combinations(rest, 4):
                terms = tuple(_mul(data.values[J | {i, j}], data.values[J | {k, l}])
                              for i, j, k, l in
                              ((a, b, c, d), (a, c, b, d), (a, d, b, c)))
                n = sum(t is not None for t in terms)
                checked += 1
                stats[n]["checked"] += 1
                if not _null(*terms):
                    failures += 1
                    stats[n]["failures"] += 1
                    if first is None:
                        first = {"J": common, "abcd": (a, b, c, d),
                                 "terms": terms,
                                 "term_labels": tuple(p_label(t) for t in terms)}
    return {"holds": failures == 0, "relations_checked": checked,
            "violations": failures, "by_nonzero_terms": stats,
            "first_failure": first}


def check_support_matroid(E, r, chi):
    """Check nonempty support and the basis-exchange axiom exhaustively.

    chi uses this module's P encoding; in particular P_ONE=0 is nonzero.
    Every ordered pair B1,B2 of supported r-subsets and e in B1-B2 is tested.
    """
    return _support(_Data(E, r, chi))


def check_three_term_plucker(E, r, chi):
    """Check all three-term relations, WITHOUT assuming matroidal support.

    For each (r-2)-subset J and unordered 4-subset {a,b,c,d} disjoint from J,
    test the products chi(Jab)chi(Jcd), chi(Jac)chi(Jbd), chi(Jad)chi(Jbc).
    Signs do not matter because -1=1. Permutations only permute the products.
    Repeated-index cases reduce to (0,0,0) or (0,t,t), hence are automatic.
    Include relations with zero products. For r<2 there are none to test.
    Read report['holds']; an all-zero chi passes here but is NOT a P-matroid.
    """
    return _plucker(_Data(E, r, chi))


def check_p_matroid(E, r, chi):
    """Check a P-representation in the sense of Baker--Lorscheid.

    Return a report with the verdict in ['is_p_matroid'] and separate
    ['support'] / ['plucker'] subreports. Do not use bool(report) as a verdict.
    chi may be callable on frozensets or a complete mapping on r-subsets.
    Alternation is supplied by symmetry on distinct tuples (-1=1) and by
    assigning pasture zero to tuples with repetitions. No higher-term or
    multiple-exchange relation is checked. Each coordinate is evaluated once.
    """
    data = _Data(E, r, chi)
    support = _support(data)
    plucker = _plucker(data)
    return {"is_p_matroid": support["is_matroid"] and plucker["holds"],
            "exact": True, "rank": r, "coordinates": len(data.subsets),
            "support_size": len(data.bases), "support": support, "plucker": plucker}
