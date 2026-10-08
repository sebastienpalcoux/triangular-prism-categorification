"""Certify that the displayed F210 local system has exactly 14 complex points.

The certificate gives a squarefree degree-14 polynomial p(t), with t=V6, and
nine rational-polynomial coordinates. Substitution verifies all 12 inputs
modulo p. The separately certified 14 spanning monomials bound the number of
points from above. Neither a root approximation nor Groebner-basis search is
required. Python standard library only.
"""
from fractions import Fraction
from pathlib import Path
import json
from polynomials import parse, standard_monomials

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def remainder(a, b):
    a = list(a)
    require(bool(b), 'Polynomial division by zero')
    while len(a) >= len(b):
        shift, factor = len(a)-len(b), a[-1]/b[-1]
        for i, c in enumerate(b):
            a[i+shift] -= factor*c
        trim(a)
    return a


def multiply(a, b, modulus):
    if not a or not b:
        return []
    c = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return remainder(trim(c), modulus)


def add_into(a, b, scale=1):
    a.extend([Fraction(0)]*max(0, len(b)-len(a)))
    for i, c in enumerate(b):
        a[i] += scale*c
    trim(a)


def univariate(poly, parameter):
    require(all(sum(e) == e[parameter] for e in poly), 'Expected a univariate polynomial')
    degree = max((e[parameter] for e in poly), default=-1)
    result = [Fraction(0)]*(degree+1)
    for e, c in poly.items():
        result[e[parameter]] = Fraction(c)
    return trim(result)


def verify(basis=None):
    if not __debug__:
        raise RuntimeError('Run without -O: the upstream membership checker uses assertions')
    data = json.loads((ROOT/'data/f210-equations.json').read_text())
    names = data['original_variables']
    require(names == ['U0','U1','U2','V0','V1','V2','V3','V4','V5','V6'], 'Unexpected coordinate order')
    if basis is None:
        # These exact ideal-membership identities, not a basis assertion,
        # supply the upper bound on the quotient dimension.
        from verify_f210 import verify_elimination, verify_membership
        basis = verify_membership(data, verify_elimination(data))
    bound = len(standard_monomials(basis))
    require(bound == 14, 'Unexpected upper bound on the quotient dimension')
    cert = json.loads((ROOT/'certificates/f210-solutions.json').read_text())
    require(cert['parameter'] == 'V6', 'Unexpected parameter')
    parameter = names.index(cert['parameter'])
    p = univariate(parse(cert['univariate'], names), parameter)
    degree = len(p)-1
    require(degree == 14, 'The parameter polynomial must have degree 14')
    derivative = trim([i*p[i] for i in range(1, len(p))])
    a, b = p, derivative
    while b:
        a, b = b, remainder(a, b)
    require(len(a) == 1, 'The parameter polynomial is not squarefree')
    coordinates = cert['coordinate_relations']
    require(set(coordinates) == set(names)-{'V6'}, 'Missing or duplicate coordinate relation')
    images = [None]*len(names)
    images[parameter] = [Fraction(0), Fraction(1)]
    for name, text in coordinates.items():
        index = names.index(name)
        relation = parse(text, names)
        linear = tuple(int(i == index) for i in range(len(names)))
        coefficient = relation.pop(linear, 0)
        require(coefficient != 0, f'Missing nonzero linear coefficient of {name}')
        rest = univariate(relation, parameter)
        images[index] = remainder([-c/coefficient for c in rest], p)
    for i, text in enumerate(data['original_equations'], 1):
        value = []
        for exponents, coefficient in parse(text, names).items():
            term = [Fraction(coefficient)]
            for coordinate, exponent in zip(images, exponents):
                for _ in range(exponent):
                    term = multiply(term, coordinate, p)
            add_into(value, term)
        require(not value, f'Original equation {i} does not vanish modulo p')
    result = {'distinct_solutions': degree, 'spanning_upper_bound': bound,
              'squarefree': True, 'input_equations_verified': len(data['original_equations'])}
    print('F210 local system: exactly 14 distinct solutions: PASS', json.dumps(result))
    return result


if __name__ == '__main__':
    verify()
