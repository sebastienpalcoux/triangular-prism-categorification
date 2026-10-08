"""Exact zero/one-spectrum criteria for general, possibly noncommutative rings.

All indices are zero-based and N[i][j][k] is the coefficient of k in i*j.
The six required nonzero coefficients prune the search; no symmetry quotient
or commutativity assumption is used. A failed search is not categorification.
"""
from functools import lru_cache
from itertools import product
from fusion_ring import FusionRingError, validate_fusion_ring


def evaluate_witness(N, indices, criterion, *, dual=None):
    """Evaluate every printed hypothesis directly, independently of pruning.

    indices are [i1,...,i9] for zero, or [i0,i1,...,i9] for one. This function
    assumes already validated fusion rules if dual is provided.
    """
    if criterion not in ('zero', 'one'):
        raise FusionRingError('criterion must be zero or one')
    if dual is None:
        dual = validate_fusion_ring(N)['dual']
    n = len(N)
    expected = 9 if criterion == 'zero' else 10
    if len(indices) != expected or any(type(i) is not int or not 0 <= i < n for i in indices):
        raise FusionRingError(f'{criterion} witness requires {expected} indices in [0, {n})')
    i1,i2,i3,i4,i5,i6,i7,i8,i9 = indices[-9:]
    i0 = indices[0] if criterion == 'one' else None
    dot = lambda a,b,c,d: sum(N[a][b][k]*N[c][d][k] for k in range(n))
    six = [N[i4][i1][i6], N[i5][i4][i2], N[i5][i6][i3],
           N[i7][i9][i1], N[i2][i7][i8], N[i8][i9][i3]]
    terms = [N[i4][i7][k]*N[dual[i5]][i8][k]*N[i6][dual[i9]][k] for k in range(n)]
    if criterion == 'zero':
        alternatives = [
            [dot(i5,i6,i2,i1), dot(i5,i4,i3,dual[i1]), dot(i2,dual[i4],i3,dual[i6])],
            [dot(i2,i1,i8,i9), dot(i2,i7,i3,dual[i9]), dot(i8,dual[i7],i3,dual[i1])]]
    else:
        alternatives = [
            [dot(i5,i0,i2,i7), dot(i5,i4,i8,dual[i7]), dot(i2,dual[i4],i8,dual[i0])],
            [dot(i5,i6,i8,i9), dot(i5,i0,i3,dual[i9]), dot(i8,dual[i0],i3,dual[i6])],
            [dot(i4,i1,i0,i9), dot(i4,i7,i6,dual[i9]), dot(i0,dual[i7],i6,dual[i1])]]
    target = int(criterion == 'one')
    conditions = {'six_coefficients_nonzero': all(six),
                  'spectrum_sum': sum(terms) == target,
                  'opposite_coefficient': N[i2][i1][i3] == 1-target,
                  'each_alternative_group_contains_one': all(1 in group for group in alternatives)}
    if criterion == 'one':
        spectrum_factors = [N[i4][i7][i0], N[dual[i5]][i8][i0], N[i6][dual[i9]][i0]]
        conditions['specified_unique_spectrum_index'] = spectrum_factors == [1,1,1] and sum(terms) == 1
    else:
        spectrum_factors = None
    return {'criterion': criterion, 'indices': list(indices),
            'index_order': ['i'+str(i) for i in range(1-target, 10)],
            'six_required_coefficients': six, 'spectrum_terms': terms,
            'spectrum_factors_at_i0': spectrum_factors,
            'spectrum_sum': sum(terms), 'opposite_coefficient': N[i2][i1][i3],
            'alternative_sums': alternatives,
            'alternative_order': 'The alternative groups in the order printed in the theorem',
            'conditions': conditions, 'passes': all(conditions.values())}


def scan(N, criterion='both', mode='count', *, validate=True):
    """Search either or both criteria; default preserves scan(N)'s old API.

    mode='count' counts every witness and retains the first detailed witness.
    mode='all' additionally retains every detailed witness (potentially large).
    mode='first' stops once each requested criterion has a witness; counts
    then describe only the visited tuples. ``exhaustive`` states whether the
    scan finished. Absence of a witness never establishes categorifiability.
    """
    if criterion not in ('zero', 'one', 'both') or mode not in ('first', 'count', 'all'):
        raise FusionRingError('criterion must be zero/one/both and mode first/count/all')
    if validate:
        checked = validate_fusion_ring(N)
        n, dual = checked['rank'], checked['dual']
    else:
        n = len(N)
    I = range(n)
    if not validate:
        dual = [next(j for j in I if N[i][j][0]) for i in I]
    selected = {0,1} if criterion == 'both' else {int(criterion == 'one')}
    support = [[tuple(k for k in I if N[i][j][k]) for j in I] for i in I]
    @lru_cache(None)
    def dot(a,b,c,d):
        return sum(N[a][b][k]*N[c][d][k] for k in I)
    def one_of(*keys):
        return any(dot(*key) == 1 for key in keys)
    @lru_cache(None)
    def spectrum(i4,i7,i5s,i8,i6,i9s):
        terms = [N[i4][i7][k]*N[i5s][i8][k]*N[i6][i9s][k] for k in I]
        total = sum(terms)
        return total, terms.index(1) if total == 1 else None
    def witnesses():
        for i4,i5,i7,i9 in product(I, repeat=4):
            for i1 in support[i7][i9]:
                for i2 in support[i5][i4]:
                    for i6 in support[i4][i1]:
                        for i8 in support[i2][i7]:
                            total,i0 = spectrum(i4,i7,dual[i5],i8,i6,dual[i9])
                            if total not in selected:
                                continue
                            if total == 1:
                                if not one_of((i5,i0,i2,i7),(i5,i4,i8,dual[i7]),(i2,dual[i4],i8,dual[i0])):
                                    continue
                                if not one_of((i4,i1,i0,i9),(i4,i7,i6,dual[i9]),(i0,dual[i7],i6,dual[i1])):
                                    continue
                            for i3 in support[i5][i6]:
                                if not N[i8][i9][i3] or N[i2][i1][i3] != 1-total:
                                    continue
                                if total == 0:
                                    if not one_of((i5,i6,i2,i1),(i5,i4,i3,dual[i1]),(i2,dual[i4],i3,dual[i6])):
                                        continue
                                    if not one_of((i2,i1,i8,i9),(i2,i7,i3,dual[i9]),(i8,dual[i7],i3,dual[i1])):
                                        continue
                                elif not one_of((i5,i6,i8,i9),(i5,i0,i3,dual[i9]),(i8,dual[i0],i3,dual[i6])):
                                    continue
                                yield total, ([i0] if total else [])+[i1,i2,i3,i4,i5,i6,i7,i8,i9]
    first = {k: None for k in sorted(selected)}
    counts = {k: 0 for k in sorted(selected)}
    details = {k: [] for k in sorted(selected)}
    exhaustive = True
    for kind, indices in witnesses():
        counts[kind] += 1
        if first[kind] is None or mode == 'all':
            record = evaluate_witness(N, indices, 'one' if kind else 'zero', dual=dual)
            if not record['passes']:
                raise RuntimeError('Search produced a witness that failed direct verification')
            details[kind].append(record)
        if first[kind] is None:
            first[kind] = indices
        if mode == 'first' and all(first[k] is not None for k in selected):
            exhaustive = False
            break
    return {'rank': n, 'dual': dual, 'criterion': criterion, 'mode': mode,
            'witness_count': counts, 'first_witness': first, 'witness_details': details,
            'obstruction_found': any(counts.values()),
            'spectrum_cache_entries': spectrum.cache_info().currsize,
            'exhaustive': exhaustive, 'witness_count_is_complete': exhaustive,
            'interpretation': ('A verified witness excludes categorification over every algebraically closed field. '
                               'No witness does not establish categorifiability.')}


if __name__ == '__main__':
    from criteria import main
    main()
