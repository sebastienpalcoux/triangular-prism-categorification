"""Count marked tetrahedral and ordered associator entries in the local F210 prisms.

The five drawings in the simplified TPE use the edge order
LR, RT, TL, LO, TO, RO (left, right, top and interior vertices).
A coefficient key consists of its drawing type and its ordered edge labels.
The four marked corners in a drawing are fixed: the two middle-cut factors
have one common corner pattern; the three side-cut factors have another.
When the arrows are retained the two middle-cut factors have different
arrow patterns. All three side-cut factors use the same arrows.

We use coherent bases for the one-dimensional vertex spaces here. Every
admissible triple encountered has a repeated label. Over C, normalize a
vector in the palindromic order to be its own dual and use cyclic transports
for the other orders as in the proof of Theorem 6.1. The fixed corner of a
triple of equal labels is retained. No move of a marked corner, rotation
or reflection of a tetrahedron is used to match coefficient keys. Unit-edge
coefficients are counted. This routine is for these multiplicity-one,
self-dual F210 instances, not for general fusion rings.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

UNIT = 0


def prism_edges(center, selected):
    """The X1,...,X9 edge labels of the three specializations in Corollary 6.2."""
    k = center
    for a, b in product(selected, repeat=2):
        yield (a, b, UNIT, k, k, k, k, k, k)
        yield (a, b, b, k, k, k, k, k, k)
        yield (k, k, a, b, k, k, b, k, k)


# Entries refer to the manuscript's X_1,...,X_9; 0 denotes the summation
# label X. The order in each row is LR,RT,TL,LO,TO,RO in the printed drawing.
MIDDLE_CUT = ((6, 4, 5, 3, 2, 1), (9, 8, 7, 1, 2, 3))
SIDE_CUT = ((1, 7, 4, 6, 0, 9), (2, 8, 5, 4, 0, 7), (3, 9, 6, 5, 0, 8))


def _polynomial_key(terms):
    """Collect commuting monomials and identify nonzero scalar multiples."""
    nonzero = sorted((monomial, coeff) for monomial, coeff in terms.items() if coeff)
    if not nonzero:
        return ()
    first = nonzero[0][1]
    return tuple((monomial, Fraction(coeff, first)) for monomial, coeff in nonzero)


def direct_prism_counts(fusion, dimensions, instances, *, arrows=False):
    """Count the fixed nine-edge prism instances without tetrahedral symmetry.

    The iterable 'instances' contains the nine labels in order X_1,...,X_9,
    with tensor unit at label 0. Dimensions must be exact integers/rationals.
    The result includes the nontrivial coincidences between different
    nine-edge labellings, making the polynomial count auditable.
    """
    scalar_keys = set()
    equations = defaultdict(set)
    prism_labels = set()
    indexed = 0
    nonzero_instances = 0
    for labels in instances:
        labels = tuple(labels)
        if len(labels) != 9:
            raise ValueError("each prism must specify nine edge labels")
        prism_labels.add(labels)
        indexed += 1
        x = (None, *labels)

        # Check that label keys need no multiplicity index or arbitrary
        # identification of vectors of higher-dimensional vertex spaces.
        def multiplicity(p, q, r):
            mult = fusion[p][q][r]
            if mult and (mult != 1 or len({p, q, r}) > 2):
                raise ValueError(f"unsupported vertex type {(p, q, r)}")
            return mult

        external = ((4, 1, 6), (5, 2, 4), (6, 3, 5),
                    (9, 1, 7), (7, 2, 8), (8, 3, 9))
        if not all(multiplicity(x[h], x[j], x[k]) for h, j, k in external):
            raise ValueError(f"inadmissible external prism labelling {labels}")

        # Each key records one drawing exactly as printed. Scalar products
        # commute, so sort factors in a monomial, never edges within a key.
        def coefficient(pattern, edge_slots, interior=None):
            key = (pattern, tuple(x[j] if j else interior for j in edge_slots))
            scalar_keys.add(key)
            return key

        terms = Counter()
        if multiplicity(x[3], x[2], x[1]):
            patterns = ("left", "right") if arrows else ("middle", "middle")
            factors = [coefficient(p, t) for p, t in zip(patterns, MIDDLE_CUT)]
            terms[tuple(sorted(factors))] += 1
        for i, dim in enumerate(dimensions):
            if all(multiplicity(x[j], x[j+3], i) for j in (4, 5, 6)):
                factors = [coefficient("side", t, i) for t in SIDE_CUT]
                terms[tuple(sorted(factors))] -= dim
        key = _polynomial_key(terms)
        if key:
            nonzero_instances += 1
            equations[key].add(labels)

    by_pattern = Counter(p for p, edges in scalar_keys)
    nonunit_by_pattern = Counter(p for p, edges in scalar_keys if 0 not in edges)
    return {
        "indexed": indexed,
        "distinct_prisms": len(prism_labels),
        "nonzero_instances": nonzero_instances,
        "entries_by_pattern": dict(sorted(by_pattern.items())),
        "nonunit_entries_by_pattern": dict(sorted(nonunit_by_pattern.items())),
        "entries": len(scalar_keys),
        "polynomials": len(equations),
        "coincident_prisms": sorted(tuple(sorted(labels)) for labels in equations.values()
                                    if len(labels) > 1),
    }


def pentagon_terms(fusion, edges):
    """Expand the two routes between (((ab)c)d) and (a(b(cd))).

    An F entry is (three inputs, output, left channel, right channel).
    All fusion spaces actually used have multiplicity one. The assertions
    check this on the boundary and in each nonempty internal summand.
    """
    # Theorem 5.3: (X1,...,X9) = (i3*,i1,i2,i6,i4,i5*,i9*,i7,i8).
    # Here a=i5, b=i4, c=i7 and d=i9. These simple objects are self-dual.
    total, cd, ab, bcd, b, a, d, c, abc = edges
    boundary = ((a, b, ab), (ab, c, abc), (abc, d, total),
                (c, d, cd), (b, cd, bcd), (a, bcd, total))
    assert all(fusion[x][y][z] == 1 for x, y, z in boundary)

    # The two-move route passes through (ab)(cd).
    middle = fusion[ab][cd][total]
    if middle:
        assert middle == 1
        yield 1, ((ab, c, d, total, abc, cd), (a, b, cd, total, ab, bcd))

    # The three-move route sums over the channel t for bc.
    for t in range(len(fusion)):
        multiplicity = fusion[b][c][t] * fusion[a][t][abc] * fusion[t][d][bcd]
        if multiplicity:
            assert multiplicity == 1
            yield -1, ((a, b, c, abc, ab, t), (a, t, d, total, abc, bcd),
                       (b, c, d, bcd, t, cd))


def polynomial(terms, normalize_unit):
    coefficients = Counter()
    for sign, factors in terms:
        # Only the scalar factors commute. Never permute the labels of an F.
        monomial = tuple(sorted(entry for entry in factors
                                if not normalize_unit or UNIT not in entry[:3]))
        coefficients[monomial] += sign
    return _polynomial_key(coefficients)


def count(fusion, instances, normalize_unit):
    rows = [polynomial(pentagon_terms(fusion, edges), normalize_unit)
            for edges in instances]
    distinct = {row for row in rows if row}
    variables = {entry for row in distinct for term, _ in row for entry in term}
    return len(variables), len(distinct), sum(bool(row) for row in rows)


def verify(fusion):
    assert all(fusion[i][i][UNIT] == 1 for i in range(len(fusion)))
    for center, selected in ((1, (0, 1, 3)), (3, (0, 2, 3))):
        instances = list(prism_edges(center, selected))
        assert (len(instances), len(set(instances))) == (27, 22)
        dimensions = (1, 5, 5, 5, 6, 7, 7)
        for arrows, expected_patterns, expected_rows in (
                (True, {"left": 15, "right": 15, "side": 25}, 19),
                (False, {"middle": 21, "side": 25}, 16)):
            report = direct_prism_counts(fusion, dimensions, instances, arrows=arrows)
            assert report["entries_by_pattern"] == expected_patterns, report
            assert report["polynomials"] == expected_rows, report
            assert report["nonzero_instances"] == 27, report
            convention = "arrowed" if arrows else "unarrowed"
            print(f'F210 center {center}, {convention} TPE: {report["entries"]} marked entries, '
                  f'{report["polynomials"]} distinct nonzero equations')
        raw = count(fusion, instances, normalize_unit=False)
        normalized = count(fusion, instances, normalize_unit=True)
        assert raw == (58, 22, 27), raw
        assert normalized == (45, 19, 24), normalized
        print(f'F210 center {center}, ordered F entries / distinct nonzero pentagons: '
              f'{raw[0]}/{raw[1]} before unit normalization, '
              f'{normalized[0]}/{normalized[1]} after unit normalization '
              f'({normalized[2]} of {len(instances)} indexed instances nonzero)')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run without -O: verification uses assertions.')
    path = Path(__file__).resolve().parent / 'data' / 'f210.json'
    verify(json.loads(path.read_text()))
