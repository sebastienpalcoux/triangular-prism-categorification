#!/usr/bin/env python3
"""Check the imported rank <= 8 census and 27 exact primary-3 witnesses.

Run from any directory.  Python 3.10+ standard library only.
This checks the listed rings and their exclusions, not exhaustive enumeration.
"""
import json
from math import lcm
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_ring(row):
    M, d = row['fusion_matrices'], row['dimensions']
    n = len(M)
    require(len(d) == n and d[0] == 1, 'dimension-vector length/unit')
    require(all(type(x) is int and x > 0 for x in d), 'positive integer dimensions')
    require(all(len(A) == n and all(len(v) == n for v in A) for A in M), 'matrix sizes')
    require(all(type(x) is int and x >= 0 for A in M for v in A for x in v), 'nonnegative integers')
    require(all(M[0][i][j] == M[i][0][j] == int(i == j) for i in range(n) for j in range(n)), 'unit')
    star = []
    for i in range(n):
        duals = [j for j in range(n) if M[i][j][0]]
        require(len(duals) == 1 and M[i][duals[0]][0] == 1, 'unique dual')
        star.append(duals[0])
    require(all(star[star[i]] == i for i in range(n)), 'duality involution')
    for i in range(n):
        for j in range(n):
            require(sum(M[i][j][k] * d[k] for k in range(n)) == d[i] * d[j], 'dimension character')
            for k in range(n):
                require(M[i][j][k] == M[j][i][k], 'commutativity')
                require(M[i][j][k] == M[star[i]][k][j] == M[k][star[j]][i], 'Frobenius reciprocity')
                require(M[i][j][k] == M[star[j]][star[i]][star[k]], 'duality on products')
                for ell in range(n):
                    require(sum(M[i][j][s] * M[s][k][ell] for s in range(n)) ==
                            sum(M[j][k][s] * M[i][s][ell] for s in range(n)), 'associativity')
    D = sum(x * x for x in d)
    require(D <= 20000 and all(D % x == 0 for x in d), 'dimension bound/1-Frobenius')
    require(all(x > 1 for x in d[1:]), 'nonpointed perfect type')
    # Any nontrivial fusion subring contains a nonunit basis element.  Its
    # fusion-subring closure must therefore be the whole basis in a simple ring.
    for seed in range(1, n):
        S = {0, seed}
        while True:
            T = S | {star[i] for i in S} | {
                k for i in S for j in S for k in range(n) if M[i][j][k]
            }
            if T == S:
                break
            S = T
        require(len(S) == n, 'simplicity')


def tensor_cube_form(A, q):
    """q^T (A tensor A tensor A) q, using Python integer contractions."""
    n = len(A)
    idx = lambda a, b, c: (a * n + b) * n + c
    u = [0] * (n ** 3)
    v = [0] * (n ** 3)
    w = [0] * (n ** 3)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                u[idx(a,b,c)] = sum(A[c][z] * q[idx(a,b,z)] for z in range(n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                v[idx(a,b,c)] = sum(A[b][y] * u[idx(a,y,c)] for y in range(n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                w[idx(a,b,c)] = sum(A[a][x] * v[idx(x,b,c)] for x in range(n))
    return sum(x * y for x, y in zip(q, w))


def verify():
    rows = json.loads((HERE / 'rank_at_most_8.json').read_text())
    certs = json.loads((HERE / 'primary3_negative_witnesses.json').read_text())
    require(len(rows) == 33, '33 imported rings expected')
    by_id = {row['id']: row for row in rows}
    require(len(by_id) == 33, 'unique identifiers')
    for row in rows:
        check_ring(row)
    survivors = {'S60_0', 'S168_0', 'S210_1', 'S360_1', 'S660_5', 'S660_11'}
    require(len(certs) == 27 and {c['id'] for c in certs} == set(by_id) - survivors,
            'exactly the 27 specified exclusions')
    for c in certs:
        row = by_id[c['id']]
        M, d, q = row['fusion_matrices'], row['dimensions'], c['integer_vector']
        L = lcm(*d)
        require(c['clearing_denominator'] == L, 'clearing denominator')
        require(len(q) == len(M) ** 3 and all(type(x) is int for x in q), 'integer witness size')
        value = sum((L // di) * tensor_cube_form(A, q) for A, di in zip(M, d))
        require(value == c['quadratic_form'] and value < 0, 'strictly negative exact quadratic form')
    # Identify the two inputs of the present article without a relabelling.
    for filename, key in [('f210.json', 'S210_1'), ('f660.json', 'S660_5')]:
        original = json.loads((HERE.parent / 'data' / filename).read_text())
        require(original == by_id[key]['fusion_matrices'], f'{filename} agrees with census')
    groups = json.loads((HERE / 'group_character_rings.json').read_text())
    expected_groups = {'S60_0': ('A5', 60), 'S168_0': ('PSL2_7', 168),
                       'S360_1': ('A6', 360), 'S660_11': ('PSL2_11', 660)}
    require(len(groups) == 4 and {g['census_id'] for g in groups} == set(expected_groups),
            'four expected group character rings')
    for g in groups:
        row = by_id[g['census_id']]
        M, C, p = row['fusion_matrices'], g['fusion_matrices'], g['census_to_group_basis']
        n = len(M)
        require((g['group'], g['order']) == expected_groups[g['census_id']], 'group name/order')
        require(sorted(p) == list(range(n)) and p[0] == 0, 'unit-preserving basis bijection')
        require(all(row['dimensions'][i] == g['dimensions'][p[i]] for i in range(n)),
                'group dimensions match')
        require(all(M[i][j][k] == C[p[i]][p[j]][p[k]]
                    for i in range(n) for j in range(n) for k in range(n)),
                'group character tensor matches census')
    print('PASS: 33 imported fusion rings satisfy the checked exact axioms, dimensions, and simplicity.')
    print('PASS: 27 exact negative primary-3 quadratic forms exclude unitary categorification.')
    print('PASS: F210 and F660 agree entry by entry with the census inputs.')
    print('PASS: four remaining inputs match character-ring tensors generated independently with GAP 4.12.1.')
    print('Enumeration completeness is imported from the cited census; group tensor regeneration is optional GAP.')


if __name__ == '__main__':
    verify()
