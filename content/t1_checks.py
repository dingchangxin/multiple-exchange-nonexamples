"""Exact exhaustive checks for strong T1-matroids and ME1 / ME2.

T1 is the triangular hyperfield: a (+) b = [abs(a-b), a+b].
Only the Python standard library is used. Integers, Fraction, Decimal and
rational/decimal strings are accepted; floats are deliberately rejected.

chi is a callable on frozensets, or a COMPLETE mapping on r-element subsets.
Zero coordinates are allowed. Missing mapping entries are errors, not zeros.
All r-subsets are checked, including nonbases. Each function snapshots chi
once, so the supplied callable should be deterministic.

ME functions default to the paper's I1 subset B1, with repeated-element
coordinates interpreted as zero, max(empty)=0 and empty hypersum={0}.
Use disjoint_only=True to impose I1 subset B1 - B2 instead. Without this
restriction ME1 fails trivially whenever rank>0 and chi is nonzero.

Reference for strong GP axioms: Baker--Bowler, Matroids over hyperfields,
arXiv:1601.01204v5, Definition 3.11 and Theorem 3.13.
"""
from collections import Counter
from collections.abc import Mapping
from decimal import Decimal
from fractions import Fraction
from itertools import combinations
from numbers import Rational

__all__ = ["check_strong_t1", "check_ME1", "check_ME2",
           "exchange_witness", "t1_hypersum_interval", "t1_contains"]


def _exact(x):
    """Validate and normalize a nonnegative exact rational number."""
    if not isinstance(x, (Rational, Decimal, str)):
        raise TypeError("Use int, Fraction, Decimal or a rational string; "
                        "floats and symbolic irrationals are not supported.")
    try:
        q = Fraction(x)
    except (ValueError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError("Coordinates must be finite rational numbers.") from exc
    if q < 0:
        raise ValueError("T1 coordinates must be nonnegative.")
    return q.numerator if q.denominator == 1 else q


def _interval(terms):
    total = sum(terms)
    largest = max(terms, default=0)
    return max(0, 2 * largest - total), total


def t1_hypersum_interval(terms):
    """Endpoints of the CLOSED T1 hypersum; the empty hypersum is {0}."""
    return _interval(tuple(_exact(x) for x in terms))


def t1_contains(target, terms):
    """Whether target belongs to the hypersum of terms (exactly)."""
    target = _exact(target)
    lo, hi = t1_hypersum_interval(terms)
    return lo <= target <= hi


class _Data:
    def __init__(self, E, r, chi):
        self.E = tuple(E)
        try:
            self.ground = frozenset(self.E)
        except TypeError as exc:
            raise TypeError("Ground-set elements must be hashable.") from exc
        if len(self.ground) != len(self.E):
            raise ValueError("E must not contain repeated elements.")
        if isinstance(r, bool) or not isinstance(r, int) or not 0 <= r <= len(self.E):
            raise ValueError("r must be an integer with 0 <= r <= len(E).")
        self.r = r
        self.subsets = tuple(frozenset(B) for B in combinations(self.E, r))
        if isinstance(chi, Mapping):
            normalized = {}
            for key, value in chi.items():
                B = tuple(key)
                S = frozenset(B)
                if len(B) != r or len(S) != r or not S <= self.ground:
                    raise ValueError("Every mapping key must be an r-subset of E.")
                if S in normalized:
                    raise ValueError("Duplicate subset keys after normalization.")
                normalized[S] = value
            if set(normalized) != set(self.subsets):
                raise ValueError("Provide every r-subset explicitly, including zeros.")
            self.values = {B: _exact(normalized[B]) for B in self.subsets}
        elif callable(chi):
            self.values = {B: _exact(chi(B)) for B in self.subsets}
        else:
            raise TypeError("chi must be callable or a complete mapping.")

    def ordered(self, S):
        return tuple(e for e in self.E if e in S)

    def subset(self, B, size=None):
        B = tuple(B)
        S = frozenset(B)
        if len(S) != len(B) or not S <= self.ground:
            raise ValueError("Expected a subset of E without repetitions.")
        if size is not None and len(S) != size:
            raise ValueError("Expected a subset of size %s." % size)
        return S


def check_strong_t1(E, r, chi):
    """Check nonzeroness and ALL strong GP relations; return a report dict.

    For each I of size r+1 and J of size r-1, set
    t_e = chi(I-{e}) * chi(J+{e}) for e in I-J.
    The relation holds iff 2*max(t_e) <= sum(t_e).
    Terms with e in J are zero by alternation, not extra conditions.
    No three-term-only shortcut is used; equality is allowed.
    Rank zero and full rank have no GP relations to enumerate.

    report['is_strong'] is the verdict. report['by_terms'] groups checks
    by |I-J| (formal terms after repeated-element zeros are omitted).
    Additional zero terms from nonbases remain in that formal count.
    """
    data = _Data(E, r, chi)
    nonzero = any(data.values.values())
    checked = failures = 0
    first = None
    stats = {}
    if 0 < r < len(data.E):
        large = tuple(frozenset(I) for I in combinations(data.E, r + 1))
        small = tuple(frozenset(J) for J in combinations(data.E, r - 1))
        for I in large:
            for J in small:
                active = data.ordered(I - J)
                terms = tuple(data.values[I - {e}] * data.values[J | {e}]
                              for e in active)
                slack = sum(terms) - 2 * max(terms, default=0)
                k = len(terms)
                entry = stats.setdefault(k, {"checked": 0, "failures": 0,
                                              "minimum_slack": None})
                entry["checked"] += 1
                if entry["minimum_slack"] is None or slack < entry["minimum_slack"]:
                    entry["minimum_slack"] = slack
                checked += 1
                if slack < 0:
                    failures += 1
                    entry["failures"] += 1
                    if first is None:
                        first = {"I": data.ordered(I), "J": data.ordered(J),
                                 "elements": active, "terms": terms, "slack": slack}
    return {"is_strong": nonzero and failures == 0, "nonzero": nonzero,
            "exact": True, "coordinates": len(data.subsets),
            "support_size": sum(x != 0 for x in data.values.values()),
            "relations_checked": checked, "violations": failures,
            "by_terms": dict(sorted(stats.items())), "first_failure": first}


def _rows(data, B1, B2, I1):
    for chosen in combinations(data.ordered(B2 - B1), len(I1)):
        I2 = frozenset(chosen)
        left = (B1 - I1) | I2
        right = (B2 - I2) | I1
        # If I1 meets B2, the second alternating tuple has a repetition.
        # Its set has size < r, so its coordinate is interpreted as zero.
        a = data.values[left]
        b = data.values.get(right, 0)
        yield I2, left, right, a, b


def _witness(data, B1, B2, I1):
    rows = tuple(_rows(data, B1, B2, I1))
    products = tuple(a * b for _, _, _, a, b in rows)
    target = data.values[B1] * data.values[B2]
    largest = max(products, default=0)
    lo, hi = _interval(products)
    holds = lo <= target <= hi
    return {"B1": data.ordered(B1), "B2": data.ordered(B2),
            "I1": data.ordered(I1), "target": target,
            "products": products, "max_product": largest,
            "hypersum_interval": (lo, hi), "membership": holds,
            "ME1_hypothesis": target >= largest,
            "ME2_hypothesis": target < largest,
            "violates_ME1": target >= largest and not holds,
            "violates_ME2": target < largest and not holds,
            "rows": tuple({"I2": data.ordered(I2), "left": data.ordered(left),
                           "right": data.ordered(right), "left_value": a,
                           "right_value": b, "product": a * b}
                          for I2, left, right, a, b in rows)}


def exchange_witness(E, r, chi, B1, B2, I1):
    """Evaluate one specified exchange and return its full table and verdicts."""
    data = _Data(E, r, chi)
    B1, B2 = data.subset(B1, r), data.subset(B2, r)
    I1 = data.subset(I1)
    if not I1 <= B1:
        raise ValueError("I1 must be a subset of B1.")
    return _witness(data, B1, B2, I1)


def _check_me(E, r, chi, which, disjoint_only):
    if not isinstance(disjoint_only, bool):
        raise TypeError("disjoint_only must be True or False.")
    data = _Data(E, r, chi)
    checked = applicable = failures = 0
    first = None
    by_size = Counter()
    for B1 in data.subsets:
        for B2 in data.subsets:
            target = data.values[B1] * data.values[B2]
            pool = data.ordered(B1 - B2 if disjoint_only else B1)
            for k in range(len(pool) + 1):
                for chosen in combinations(pool, k):
                    I1 = frozenset(chosen)
                    products = tuple(a * b for _, _, _, a, b in _rows(data, B1, B2, I1))
                    largest = max(products, default=0)
                    checked += 1
                    hypothesis = (target >= largest if which == 1 else target < largest)
                    if not hypothesis:
                        continue
                    applicable += 1
                    lo, hi = _interval(products)
                    if not lo <= target <= hi:
                        failures += 1
                        by_size[k] += 1
                        if first is None:
                            first = _witness(data, B1, B2, I1)
    return {"property": "ME%s" % which, "holds": failures == 0, "exact": True,
            "disjoint_only": disjoint_only, "triples_checked": checked,
            "hypothesis_satisfied": applicable, "violations": failures,
            "violations_by_exchange_size": dict(sorted(by_size.items())),
            "first_failure": first}


def check_ME1(E, r, chi, *, disjoint_only=False):
    """Exhaustively check ME1: target >= max(products) implies membership.

    Default: I1 subset B1, as written in the paper, with the zero-extension
    conventions documented above. disjoint_only=True restricts I1 to B1-B2.
    Returns a report; read report['holds'], NOT bool(report).
    The strong GP axioms are a separate check: call check_strong_t1 first.
    """
    return _check_me(E, r, chi, 1, disjoint_only)


def check_ME2(E, r, chi, *, disjoint_only=False):
    """Exhaustively check ME2: target < max(products) implies membership.

    Same input conventions and report fields as check_ME1. All triples are
    checked (including empty I1); enumeration does not stop at a first failure.
    """
    return _check_me(E, r, chi, 2, disjoint_only)
