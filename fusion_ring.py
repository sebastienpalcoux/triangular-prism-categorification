"""Read and validate finite fusion rules, without assuming commutativity.

The convention is N[i][j][k] = multiplicity of k in i*j; the unit is 0.
JSON input is either the raw rank-by-rank-by-rank tensor, or an object with
``fusion_matrices`` and optional ``name``, ``labels``, ``dimensions``, and
``unit_index`` (which must be 0). Other metadata is retained. If supplied,
``dimensions`` is a positive rational dimension character, represented by
integers or strings such as "3/2". Omit it for irrational dimensions.
It is not the potentially signed ``categorical_dimensions`` of localization.
"""
from fractions import Fraction
from pathlib import Path
import json
import re


class FusionRingError(ValueError):
    """Input does not satisfy the documented exact fusion-ring convention."""


def require(condition, message):
    if not condition:
        raise FusionRingError(message)


def rational(value):
    """An exact rational, rejecting bools, floats, and approximate numbers."""
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    if isinstance(value, str) and re.fullmatch(r'[+-]?\d+(?:/\d+)?', value):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError) as error:
            raise FusionRingError(f'Invalid rational: {value!r}') from error
    raise FusionRingError(f'Expected an integer or exact rational string, got {value!r}')


def validate_fusion_ring(N, dimensions=None):
    """Check all based-ring axioms used here, with explicit exceptions.

    Return rank, the derived dual involution, commutativity, and normalized
    optional dimensions. Validation does not establish categorifiability.
    """
    sequence = lambda x: isinstance(x, (list, tuple))
    require(sequence(N) and len(N) > 0, 'fusion_matrices must be a nonempty tensor')
    n = len(N)
    require(all(sequence(A) and len(A) == n and
                all(sequence(row) and len(row) == n for row in A) for A in N),
            f'Expected a {n} by {n} by {n} tensor')
    require(all(type(c) is int and c >= 0 for A in N for row in A for c in row),
            'Fusion coefficients must be nonnegative integers (not bools or floats)')
    I = range(n)
    require(all(N[0][i][j] == N[i][0][j] == int(i == j) for i in I for j in I),
            'Basis index 0 must be a two-sided unit')
    dual = []
    for i in I:
        candidates = [j for j in I if N[i][j][0]]
        require(len(candidates) == 1 and N[i][candidates[0]][0] == 1,
                f'Basis index {i} must have one dual with unit coefficient 1')
        dual.append(candidates[0])
    require(all(dual[dual[i]] == i for i in I), 'Duality must be an involution')
    for i in I:
        for j in I:
            for k in I:
                require(N[i][j][k] == N[dual[i]][k][j] == N[k][dual[j]][i],
                        f'Frobenius reciprocity fails at ({i}, {j}, {k})')
                require(N[i][j][k] == N[dual[j]][dual[i]][dual[k]],
                        f'Duality on products fails at ({i}, {j}, {k})')
                for ell in I:
                    left = sum(N[i][j][s] * N[s][k][ell] for s in I)
                    right = sum(N[j][k][s] * N[i][s][ell] for s in I)
                    require(left == right,
                            f'Associativity fails at ({i}, {j}, {k}, {ell}): {left} != {right}')
    normalized = None
    if dimensions is not None:
        require(sequence(dimensions) and len(dimensions) == n,
                'dimensions must have one entry per basis element')
        d = [rational(x) for x in dimensions]
        require(d[0] == 1 and all(x > 0 for x in d),
                'dimensions must be positive and have unit dimension 1')
        for i in I:
            for j in I:
                require(sum(N[i][j][k] * d[k] for k in I) == d[i] * d[j],
                        f'Dimension character fails at ({i}, {j})')
        normalized = [x.numerator if x.denominator == 1 else str(x) for x in d]
    return {'rank': n, 'dual': dual,
            'commutative': all(N[i][j] == N[j][i] for i in I for j in I),
            'dimensions': normalized}


def load_fusion_ring(path):
    """Load JSON and return (tensor, metadata); validation is a separate call."""
    def invalid_constant(value):
        raise FusionRingError(f'Non-finite JSON number is not permitted: {value}')
    data = json.loads(Path(path).read_text(), parse_constant=invalid_constant)
    if isinstance(data, dict):
        require('fusion_matrices' in data, 'JSON object requires fusion_matrices')
        N = data['fusion_matrices']
        metadata = {key: value for key, value in data.items() if key != 'fusion_matrices'}
        unit = metadata.get('unit_index', 0)
        require(type(unit) is int and unit == 0,
                'unit_index must be 0; relabel the tensor before using this program')
        if 'labels' in metadata:
            labels = metadata['labels']
            require(isinstance(N, list) and isinstance(labels, list) and
                    len(labels) == len(N) and all(isinstance(x, str) for x in labels) and
                    len(set(labels)) == len(labels),
                    'labels must be distinct strings, one per basis element')
    else:
        N, metadata = data, {}
    return N, metadata
