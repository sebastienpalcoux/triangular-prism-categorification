#!/usr/bin/env python3
"""Check the displayed F210 character table exactly in prime cyclotomic fields."""
import json
from pathlib import Path


def verify(N=None):
    if N is None:
        N = json.loads((Path(__file__).resolve().parent/'data/f210.json').read_text())
    checked = []
    for p in (2, 7, 5):
        # Z[t]/(1+t+...+t^(p-1)); these small prime cyclotomic rings
        # suffice for the explicitly displayed table, without root approximation.
        def const(a):
            return (a,)+(0,)*(p-2)

        def add(a, b):
            return tuple(x+y for x, y in zip(a, b))

        def scale(c, a):
            return tuple(c*x for x in a)

        def reduce(c):
            return tuple(c[i]-c[p-1] for i in range(p-1))

        def mul(a, b):
            c = [0]*p
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    c[(i+j) % p] += x*y
            return reduce(c)

        def zpower(i):
            c = [0]*p
            c[i % p] = 1
            return reduce(c)

        one, zero = const(1), const(0)
        if p == 2:
            columns = [list(map(const, values)) for values in
                       ([1,5,5,5,6,7,7], [1,-1,-1,-1,0,1,1])]
            codegrees = [210, 6]
        elif p == 7:
            a = [scale(-1, add(zpower(i), zpower(-i))) for i in (1,2,3)]
            columns = [[one, a[j], a[(j+1) % 3], a[(j+2) % 3], const(-1), zero, zero]
                       for j in range(3)]
            codegrees = [7]*3
        else:
            a = [add(zpower(i), zpower(-i)) for i in (1,2)]
            columns = [[one,zero,zero,zero,one,a[j],a[1-j]] for j in range(2)]
            codegrees = [5]*2
        assert len({tuple(c) for c in columns}) == len(columns)
        for column, codegree in zip(columns, codegrees):
            assert column[0] == one
            for i in range(7):
                for j in range(7):
                    rhs = zero
                    for k in range(7):
                        rhs = add(rhs, scale(N[i][j][k], column[k]))
                    assert mul(column[i], column[j]) == rhs, (p, i, j)
            actual = zero
            for value in column:  # Every basis element of F210 is self-dual.
                actual = add(actual, mul(value, value))
            assert actual == const(codegree)
            assert (zero not in column) == (codegree == 210)
            checked.append(codegree)
    # Columns in different fields have different zero patterns, except for
    # the two rational columns, which were compared above.
    assert checked == [210,6,7,7,7,5,5]
    print('F210 cyclotomic character table and dual-Burnside property: PASS')
    return checked


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Run without -O: mathematical checks use assertions')
    verify()
