"""Verify the minimum rank for the zero and one spectrum criteria.

Run alongside spectrum.py and fusion_ring.py with Python 3. The finite
lower-bound check allows arbitrary fusion multiplicities: cap x at 2.
The rank-four example is F_{3,20}, in the basis (1,x,y,z) of line 33531
of the census: https://github.com/sebastienpalcoux/fusion-ring-census/blob/afe3ffa0a371daf5131f733c9459c2eb85e18c83/results/rank4/tables.txt.gz .
"""
from itertools import combinations_with_replacement, permutations, product

from spectrum import evaluate_witness, scan


def selfdual_tensor(rank, entries):
    """Fill the symmetric tensor in lexicographic order of nonunit triples."""
    labels = range(rank)
    rules = [[[0 for _ in labels] for _ in labels] for _ in labels]
    for i in labels:
        rules[0][i][i] = rules[i][0][i] = rules[i][i][0] = 1
    triples = list(combinations_with_replacement(range(1, rank), 3))
    if len(entries) != len(triples):
        raise ValueError('Wrong number of nonunit coefficients')
    for triple, value in zip(triples, entries):
        for a, b, c in permutations(triple):
            rules[a][b][c] = value
    return rules


def associative_after_capping(rules):
    labels = range(len(rules))
    for a, b, c, d in product(labels, repeat=4):
        left = sum(rules[a][b][i] * rules[i][c][d] for i in labels)
        right = sum(rules[b][c][i] * rules[a][i][d] for i in labels)
        if min(left, 2) != min(right, 2):
            return False
    return True


def main():
    for rank, parameters, expected in [(1, 0, 1), (2, 1, 3), (3, 4, 50)]:
        checked = 0
        for entries in product(range(3), repeat=parameters):
            rules = selfdual_tensor(rank, entries)
            if not associative_after_capping(rules):
                continue
            # True associativity over N is not needed for this finite
            # superset of all truncated fusion tensors.
            report = scan(rules, mode='count', validate=False)
            if report['obstruction_found'] or not report['exhaustive']:
                raise RuntimeError((rank, entries, report))
            checked += 1
        if checked != expected:
            raise RuntimeError((rank, checked, expected))
        print(f'Rank {rank}: {checked} associative capped self-dual tables; no witness.')
    cyclic = [[[int((a + b) % 3 == c) for c in range(3)]
               for b in range(3)] for a in range(3)]
    if scan(cyclic)['obstruction_found']:
        raise RuntimeError('Unexpected witness in Z[C3]')
    print('The only non-self-dual rank-three ring is Z[C3]; no witness.')
    rules = selfdual_tensor(4, [1, 0, 1, 1, 1, 2, 1, 0, 2, 3])
    for criterion, indices in [('zero', [2, 2, 1, 1, 2, 3, 2, 2, 2]),
                               ('one', [2, 3, 2, 2, 1, 2, 2, 2, 2, 1])]:
        record = evaluate_witness(rules, indices, criterion)
        if not record['passes']:
            raise RuntimeError(record)
        print(f'Rank 4, {criterion}: {indices}; alternative sums {record["alternative_sums"]}.')
    counts = scan(rules)['witness_count']
    if counts != {0: 12, 1: 12}:
        raise RuntimeError(counts)
    print('F_{3,20}: 12 zero spectrum and 12 one spectrum witnesses.')


if __name__ == '__main__':
    main()
