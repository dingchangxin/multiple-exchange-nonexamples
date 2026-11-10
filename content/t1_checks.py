"""Exact exhaustive checks for strong T1-matroids and ME1 / ME2 / ME3.

T1 is the triangular hyperfield: a (+) b = [abs(a-b), a+b].
Only the Python standard library is used. Integers, Fraction, Decimal and
rational/decimal strings are accepted; floats are deliberately rejected.

chi is a callable on frozensets, or a COMPLETE mapping on r-element subsets.
Zero coordinates are allowed. Missing mapping entries are errors, not zeros.
All r-subsets are checked, including nonbases. Each function snapshots chi
once, so the supplied callable should be deterministic.

ME functions always require I1 subset B1 - B2 and enumerate I2 subset
B2 - B1 of the same size. The bases B1 and B2 need not be disjoint.
Empty I1 is allowed; its exchange family contains the single choice I2=empty.

Reference for strong GP axioms: Baker--Bowler, Matroids over hyperfields,
arXiv:1601.01204v5, Definition 3.11 and Theorem 3.13.
"""
from collections import Counter
from collections.abc import Mapping
from decimal import Decimal
from fractions import Fraction
from itertools import combinations
from numbers import Rational
from math import comb, lcm

__all__ = ["check_strong_t1", "check_ME1", "check_ME2", "check_ME3",
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


def check_strong_t1(E, r, chi, *, arithmetic="fraction"):
    """Check nonzeroness and ALL strong GP relations; return a report dict.

    For each I of size r+1 and J of size r-1, set
    t_e = chi(I-{e}) * chi(J+{e}) for e in I-J.
    The relation holds iff 2*max(t_e) <= sum(t_e).
    Terms with e in J are zero by alternation, not extra conditions.
    When J is a subset of I, I-J={a,b} and both products equal
    chi(J+{a}) * chi(J+{b}). These automatic relations are skipped before
    any products are computed. All other relations are checked, not just
    three-term relations; equality is allowed.
    Rank zero and full rank have no GP relations to enumerate.

    report['is_strong'] is the verdict. report['relations_checked'] counts
    only explicitly evaluated relations; report['relations_skipped_trivial']
    counts automatic J subset I relations. Their sum is relations_total.
    report['by_terms'] groups evaluated relations by |I-J| (formal terms
    after repeated-element zeros are omitted), so it has no two-term entry.
    Additional zero terms from nonbases remain in that formal count.

    arithmetic="fraction" (default) uses the original exact coordinates.
    arithmetic="integer" clears their denominators by a common positive
    integer D. Each product and slack is multiplied by D**2, so the verdict
    is unchanged. Returned products and slacks are converted back to the
    original scale. No floating-point arithmetic or tolerance is used.
    """
    if arithmetic not in ("fraction", "integer"):
        raise ValueError("arithmetic must be 'fraction' or 'integer'.")
    data = _Data(E, r, chi)
    values = data.values
    scale = 1
    if arithmetic == "integer":
        for value in values.values():
            scale = lcm(scale, value.denominator)
        values = {B: value.numerator * (scale // value.denominator)
                  for B, value in values.items()}
    nonzero = any(values.values())
    checked = failures = skipped = 0
    first = None
    stats = {}
    if 0 < r < len(data.E):
        large = tuple(frozenset(I) for I in combinations(data.E, r + 1))
        small = tuple(frozenset(J) for J in combinations(data.E, r - 1))
        for I in large:
            for J in small:
                # I-J has two elements; the two products are identical.
                # This holds even when either coordinate is zero.
                if J <= I:
                    skipped += 1
                    continue
                active = data.ordered(I - J)
                terms = tuple(values[I - {e}] * values[J | {e}]
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
    # Restore original units in the report, including failure diagnostics.
    if scale != 1:
        square = scale * scale
        for entry in stats.values():
            entry["minimum_slack"] = Fraction(entry["minimum_slack"], square)
        if first is not None:
            first["terms"] = tuple(Fraction(t, square) for t in first["terms"])
            first["slack"] = Fraction(first["slack"], square)
    return {"is_strong": nonzero and failures == 0, "nonzero": nonzero,
            "exact": True, "coordinates": len(data.subsets),
            "support_size": sum(x != 0 for x in data.values.values()),
            "relations_checked": checked, "relations_skipped_trivial": skipped,
            "relations_total": checked + skipped, "violations": failures,
            "by_terms": dict(sorted(stats.items())), "first_failure": first}


def _rows(data, B1, B2, I1):
    for chosen in combinations(data.ordered(B2 - B1), len(I1)):
        I2 = frozenset(chosen)
        left = (B1 - I1) | I2
        right = (B2 - I2) | I1
        # I1 subset B1-B2 and I2 subset B2-B1 preserve both cardinalities.
        a = data.values[left]
        b = data.values[right]
        yield I2, left, right, a, b


def _witness(data, B1, B2, I1):
    rows = tuple(_rows(data, B1, B2, I1))
    products = tuple(a * b for _, _, _, a, b in rows)
    target = data.values[B1] * data.values[B2]
    largest = max(products)
    lo, hi = _interval(products)
    holds = lo <= target <= hi
    return {"B1": data.ordered(B1), "B2": data.ordered(B2),
            "I1": data.ordered(I1), "target": target,
            "products": products, "max_product": largest,
            "hypersum_interval": (lo, hi), "membership": holds,
            "ME1_hypothesis": target >= largest,
            "ME2_hypothesis": target < largest,
            "ME3_hypothesis": target == 0,
            "violates_ME1": target >= largest and not holds,
            "violates_ME2": target < largest and not holds,
            "violates_ME3": target == 0 and not holds,
            "rows": tuple({"I2": data.ordered(I2), "left": data.ordered(left),
                           "right": data.ordered(right), "left_value": a,
                           "right_value": b, "product": a * b}
                          for I2, left, right, a, b in rows)}


def exchange_witness(E, r, chi, B1, B2, I1):
    """Evaluate one exchange with I1 subset B1-B2; return its table and verdicts.

    B1 and B2 must be r-subsets of E. Reject I1 outside B1-B2, including
    elements of B1 intersection B2. Empty I1 is allowed.
    """
    data = _Data(E, r, chi)
    B1, B2 = data.subset(B1, r), data.subset(B2, r)
    I1 = data.subset(I1)
    if not I1 <= B1 - B2:
        raise ValueError("I1 must be a subset of B1 - B2.")
    return _witness(data, B1, B2, I1)


def _check_me(E, r, chi, which):
    data = _Data(E, r, chi)
    checked = applicable = failures = 0
    first = None
    by_size = Counter()
    for B1 in data.subsets:
        for B2 in data.subsets:
            target = data.values[B1] * data.values[B2]
            pool = data.ordered(B1 - B2)
            for k in range(len(pool) + 1):
                for chosen in combinations(pool, k):
                    I1 = frozenset(chosen)
                    products = tuple(a * b for _, _, _, a, b in _rows(data, B1, B2, I1))
                    largest = max(products)
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
            "triples_checked": checked,
            "hypothesis_satisfied": applicable, "violations": failures,
            "violations_by_exchange_size": dict(sorted(by_size.items())),
            "first_failure": first}


def check_ME1(E, r, chi):
    """Exhaustively check ME1: target >= max(products) implies membership.

    Always enumerate I1 subset B1-B2 and I2 subset B2-B1 with |I2|=|I1|.
    B1 and B2 range over all r-subsets, including nonbases and overlapping
    subsets. Empty I1 is included. Every exchanged coordinate has size r.
    Returns a report; read report['holds'], NOT bool(report).
    The strong GP axioms are a separate check: call check_strong_t1 first.
    """
    return _check_me(E, r, chi, 1)


def check_ME2(E, r, chi):
    """Exhaustively check ME2: target < max(products) implies membership.

    Always require I1 subset B1-B2, with the same conventions as check_ME1.
    All allowed triples are checked (including empty I1); enumeration does
    not stop at a first failure.
    """
    return _check_me(E, r, chi, 2)


def check_ME3(E, r, chi):
    """Check ME3: target == 0 implies membership in the exchange hypersum.

    Always use I1 subset B1-B2, including empty I1. Every allowed triple
    satisfying the hypothesis is checked; do not stop at the first failure.
    Pairs with two nonzero coordinates are omitted before forming products.
    They cannot violate this implication. The report distinguishes:
      triples_total: all allowed (B1, B2, I1), including vacuous cases;
      triples_checked: cases with chi(B1)*chi(B2) == 0;
      triples_skipped_nonzero: cases with a nonzero target.
    A nonzero function with no zero coordinates satisfies ME3 vacuously.
    Matroid validity remains a separate check: call check_strong_t1 first.
    """
    data = _Data(E, r, chi)
    zeros = tuple(B for B in data.subsets if data.values[B] == 0)
    n = len(data.E)
    total = len(data.subsets) * sum(
        comb(r, d) * comb(n - r, d) * (1 << d)
        for d in range(min(r, n - r) + 1))
    checked = failures = pairs = 0
    first = None
    by_size = Counter()
    for B1 in data.subsets:
        # This enumerates every ordered pair with at least one zero once.
        candidates = data.subsets if data.values[B1] == 0 else zeros
        for B2 in candidates:
            pairs += 1
            pool = data.ordered(B1 - B2)
            for k in range(len(pool) + 1):
                for chosen in combinations(pool, k):
                    I1 = frozenset(chosen)
                    products = tuple(a * b for _, _, _, a, b
                                     in _rows(data, B1, B2, I1))
                    checked += 1
                    if sum(products) < 2 * max(products):
                        failures += 1
                        by_size[k] += 1
                        if first is None:
                            first = _witness(data, B1, B2, I1)
    return {"property": "ME3", "holds": failures == 0, "exact": True,
            "triples_total": total, "triples_checked": checked,
            "triples_skipped_nonzero": total - checked,
            "pairs_with_zero_target": pairs, "hypothesis_satisfied": checked,
            "violations": failures,
            "violations_by_exchange_size": dict(sorted(by_size.items())),
            "first_failure": first}
