"""Exact necessary localization systems for spherical categorification.

Coefficient field: Q; solutions may lie in any algebraically closed field of
characteristic zero. Dimensions can be rational specializations or variables.
The module constructs equations; it does not decide their solvability.
"""
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd, lcm
import hashlib
import json
import re

from fusion_ring import validate_fusion_ring
from polynomials import add_into, multiply, order


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _rational(value):
    _require(type(value) is int or isinstance(value, (str, Fraction)),
             'Use integers, rational strings, or Fraction values for categorical dimensions; floats are not exact inputs.')
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('Invalid rational categorical dimension: '+repr(value)) from error


def _index(value, rank, label):
    _require(type(value) is int and 0 <= value < rank, label+' must be a zero-based basis index.')
    return value


def _encode(poly):
    return [[str(c), list(e)] for e, c in sorted(poly.items(), key=lambda pair: order(pair[0]), reverse=True) if c]


def _decode(data):
    return {tuple(e): Fraction(c) for c, e in data if Fraction(c)}


def _normalize(poly):
    """A primitive integer representative of the same rational equation."""
    if not poly:
        return {}
    denominator = lcm(*(Fraction(c).denominator for c in poly.values()))
    integral = {e: int(c*denominator) for e, c in poly.items()}
    content = reduce(gcd, integral.values())
    if integral[max(integral, key=order)] < 0:
        content = -content
    return {e: c//content for e, c in integral.items()}


def _text(poly, variables):
    terms = []
    for exponents, coefficient in sorted(poly.items(), key=lambda pair: order(pair[0]), reverse=True):
        monomial = '*'.join(name+('^'+str(power) if power != 1 else '')
                            for name, power in zip(variables, exponents) if power)
        absolute = abs(coefficient)
        factor = (str(absolute)+'*' if absolute != 1 and monomial else
                  str(absolute) if not monomial else '')
        terms.append(('-' if coefficient < 0 else '+')+factor+monomial)
    return ''.join(terms).lstrip('+') or '0'


class _Space:
    def __init__(self, names):
        self.names = names
        self.indices = {name: i for i, name in enumerate(names)}
        self.zero = (0,)*len(names)

    def constant(self, value):
        return {self.zero: Fraction(value)} if value else {}

    def variable(self, name):
        e = [0]*len(self.names)
        e[self.indices[name]] = 1
        return {tuple(e): Fraction(1)}

    @staticmethod
    def sum(terms):
        result = {}
        for term in terms:
            add_into(result, term)
        return result

    @staticmethod
    def difference(left, right):
        result = left.copy()
        add_into(result, right, -1)
        return result


class _Equations:
    def __init__(self, names):
        self.names, self.rows, self.indices, self.zeros = names, [], {}, []
        self.instances = 0

    def add(self, polynomial, origin):
        self.instances += 1
        p = _normalize(polynomial)
        if not p:
            self.zeros.append(origin)
            return
        key = tuple(sorted(p.items()))
        if key in self.indices:
            self.rows[self.indices[key]]['origins'].append(origin)
            return
        self.indices[key] = len(self.rows)
        self.rows.append({'polynomial': _encode(p), 'text': _text(p, self.names), 'origins': [origin]})


def _validate_dimensions(N, dual, values):
    if values is None:
        return None
    _require(isinstance(values, (list, tuple)) and len(values) == len(N),
             'categorical_dimensions must have one entry for each basis element.')
    d = [_rational(x) for x in values]
    _require(d[0] == 1 and all(d), 'Categorical dimensions must be nonzero and the unit dimension must be 1.')
    _require(all(d[i] == d[dual[i]] for i in range(len(N))),
             'Spherical categorical dimensions must agree on dual pairs.')
    for i, j in product(range(len(N)), repeat=2):
        _require(d[i]*d[j] == sum(N[i][j][a]*d[a] for a in range(len(N))),
                 'The supplied categorical dimensions are not a dimension character.')
    return d


def generate_localization(fusion_matrices, *, k, subset=None,
                          categorical_dimensions=None, mode='corollary',
                          odd_witness=None, namespace=None):
    """Generate Theorem 6.1 or Corollary 6.2 for one distinguished object.

    ``categorical_dimensions=None`` includes symbolic, nonzero, dual-invariant
    dimension-character equations. A rational list instead specializes those
    dimensions; inconsistency then excludes that dimension choice only.
    ``subset`` selects S'_k in corollary mode. Full mode uses all of S_k.
    """
    N = fusion_matrices
    validated = validate_fusion_ring(N)
    rank, dual = validated['rank'], validated['dual']
    k = _index(k, rank, 'k')
    _require(mode in ('corollary', 'full'), 'mode must be corollary or full.')
    _require(dual[k] == k, 'The distinguished object k must be self-dual.')
    _require(all(c <= 1 for c in N[k][k]), 'The tensor square of k must be multiplicity-free.')
    support = [i for i, c in enumerate(N[k][k]) if c]
    _require(all(dual[i] == i for i in support), 'Every summand of k squared must be self-dual.')
    # For self-dual k, N[b][c*][k] = N[c][b*][k]. Cross terms in
    # [X X*:k] therefore pair modulo 2. An odd witness X exists exactly
    # when one of its simple summands is an odd witness.
    witnesses = [b for b in range(rank) if N[b][dual[b]][k] % 2 == 1]
    _require(witnesses, 'No odd-multiplicity witness N[b][b*][k] exists.')
    if odd_witness is None:
        odd_witness = witnesses[0]
    else:
        odd_witness = _index(odd_witness, rank, 'odd_witness')
        _require(odd_witness in witnesses, 'The supplied odd-multiplicity witness fails its hypothesis.')
    if subset is None:
        selected = support[:]
    else:
        _require(isinstance(subset, (list, tuple)), 'subset must be a list of zero-based indices.')
        selected = [_index(i, rank, 'subset entry') for i in subset]
        _require(len(set(selected)) == len(selected), 'subset must not contain duplicate indices.')
        _require(set(selected) <= set(support), 'subset must lie in the support of k squared.')
        selected.sort()
    _require(mode != 'full' or selected == support,
             'Full theorem mode requires all of S_k; use corollary mode for a proper selected subset.')
    d = _validate_dimensions(N, dual, categorical_dimensions)
    namespace = namespace or 'k'+str(k)
    _require(isinstance(namespace, str) and re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', namespace),
             'namespace must be an alphanumeric identifier beginning with a letter.')

    dim_representatives = sorted({min(i, dual[i]) for i in range(1, rank)})
    variables, metadata = [], {}
    def register(name, description):
        variables.append(name)
        metadata[name] = description
    if d is None:
        for i in dim_representatives:
            register('d'+str(i), {'kind': 'categorical_dimension', 'indices': [j for j in range(rank) if min(j,dual[j]) == i]})
            register('r'+str(i), {'kind': 'inverse_dimension', 'indices': [j for j in range(rank) if min(j,dual[j]) == i]})

    if mode == 'corollary':
        xkeys = list(product(support, selected))
        free_x = [key for key in xkeys if key[0] != 0 and N[key[1]][key[1]][key[0]] and key[1] != k]
    else:
        xkeys = sorted({tuple(sorted(t)) for t in product(support, repeat=3)})
        free_x = [t for t in xkeys if 0 not in t and N[t[1]][t[2]][t[0]] and not (k in support and t.count(k) >= 2)]
    for key in free_x:
        register(namespace+'_x_'+'_'.join(map(str,key)), {'kind':'x', 'center':k, 'indices':list(key), 'mode':mode})
    ykeys = sorted({tuple(sorted((i,b))) for i,b in product(support,selected)})
    for key in ykeys:
        if 0 not in key:
            register(namespace+'_y_'+'_'.join(map(str,key)), {'kind':'y', 'center':k, 'indices':list(key)})
    space = _Space(variables)
    if d is None:
        dimensions = [space.constant(1) if i == 0 else space.variable('d'+str(min(i,dual[i]))) for i in range(rank)]
        inverses = [space.constant(1) if i == 0 else space.variable('r'+str(min(i,dual[i]))) for i in range(rank)]
    else:
        dimensions, inverses = [space.constant(c) for c in d], [space.constant(1/c) for c in d]

    def y(i,b):
        if i == 0 or b == 0:
            return inverses[k]
        return space.variable(namespace+'_y_'+'_'.join(map(str,sorted((i,b)))))
    def x(*args):
        if mode == 'corollary':
            i,b = args
            if i == 0:
                return multiply(inverses[b],inverses[k])
            if not N[b][b][i]:
                return {}
            if b == k:
                return multiply(y(i,k),y(i,k))
            key = (i,b)
        else:
            key = tuple(sorted(args))
            i,b,c = key
            if i == 0:
                return multiply(inverses[b],inverses[k]) if b == c else {}
            if not N[b][c][i]:
                return {}
            if k in support and key.count(k) >= 2:
                remaining = list(key);remaining.remove(k);remaining.remove(k)
                return multiply(y(remaining[0],k),y(remaining[0],k))
        return space.variable(namespace+'_x_'+'_'.join(map(str,key)))

    eq = _Equations(variables)
    if d is None:
        for i in dim_representatives:
            eq.add(space.difference(multiply(dimensions[i],inverses[i]),space.constant(1)), {'family':'nonzero_dimension','i':i})
        for i,j in product(range(rank),repeat=2):
            rhs = space.sum(multiply(space.constant(N[i][j][a]),dimensions[a]) for a in range(rank) if N[i][j][a])
            eq.add(space.difference(multiply(dimensions[i],dimensions[j]),rhs), {'family':'dimension_character','i':i,'j':j})
    if mode == 'corollary':
        for a,b in product(selected,repeat=2):
            terms = [multiply(dimensions[i],multiply(y(i,a),y(i,b))) for i in support]
            eq.add(space.difference(multiply(dimensions[b],space.sum(terms)),space.constant(int(a==b))), {'family':'Ia','a':a,'b':b,'center':k})
            terms = [multiply(dimensions[i],multiply(y(i,a),multiply(y(i,b),y(i,b)))) for i in support]
            eq.add(space.difference(space.sum(terms),x(a,b)), {'family':'Ib','a':a,'b':b,'center':k})
            terms = [multiply(dimensions[i],multiply(y(i,a),x(i,b))) for i in support]
            eq.add(space.difference(space.sum(terms),multiply(y(a,b),y(a,b))), {'family':'IIa','a':a,'b':b,'center':k})
    else:
        for a,b,c in product(support,repeat=3):
            terms = [multiply(dimensions[i],multiply(y(i,a),multiply(y(i,b),y(i,c)))) for i in support]
            eq.add(space.difference(space.sum(terms),x(a,b,c)), {'family':'I','a':a,'b':b,'c':c,'center':k})
            terms = [multiply(dimensions[i],multiply(y(i,a),x(i,b,c))) for i in support]
            eq.add(space.difference(space.sum(terms),multiply(y(a,b),y(a,c))), {'family':'II','a':a,'b':b,'c':c,'center':k})
    coordinate_x = {','.join(map(str,key)):_encode(x(*key)) for key in xkeys}
    coordinate_y = {','.join(map(str,(i,b))):_encode(y(i,b)) for i,b in product(support,selected)}
    return {
        'schema':'tpe-localization-1','status':'necessary_system_generated',
        'coefficient_field':'QQ','ground_characteristic':0,
        'interpretation':'Necessary equations for spherical categorification; solvability does not prove categorifiability.',
        'dimension_scope': ('All nonzero spherical dimension characters, including irrational values.' if d is None else
                            'Only the supplied categorical-dimension character; it is not inferred from Frobenius-Perron dimensions.'),
        'fusion_tensor_sha256':hashlib.sha256(json.dumps(N,separators=(',',':')).encode()).hexdigest(),
        'variables':variables,'variable_metadata':metadata,
        'dimensions':{'mode':'symbolic' if d is None else 'rational','values':None if d is None else list(map(str,d)),
                      'polynomials':list(map(_encode,dimensions)),'inverse_polynomials':list(map(_encode,inverses))},
        'localizations':[{'center':k,'mode':mode,'namespace':namespace,'support':support,'selected_support':selected,
                         'odd_witness':odd_witness,'odd_multiplicity':N[odd_witness][dual[odd_witness]][k],
                         'hypotheses_checked':['fusion-ring axioms','self-dual center','multiplicity-free square','self-dual support','odd-multiplicity witness','selected support'],
                         'coordinates':{'x':coordinate_x,'y':coordinate_y}}],
        'equations':eq.rows,'identically_zero_instances':eq.zeros,
        'counts':{'variables':len(variables),'equations':len(eq.rows),'instances':eq.instances,'identically_zero_instances':len(eq.zeros)},
        'basis_convention':'One set of vertex bases, dual bases, and cyclic transports for each local system; y is symmetric, corollary x is an ordered-pair coordinate.',
        'couplings':[]}


def _lift(poly, old_names, new_names):
    positions = [new_names.index(name) for name in old_names]
    result = {}
    for e,c in poly.items():
        f = [0]*len(new_names)
        for i,power in enumerate(e):
            f[positions[i]] = power
        result[tuple(f)] = c
    return result


def combine_localizations(systems, couplings=(), *, compatible_bases=False):
    """Combine necessary systems and optionally add ordered Theorem 6.4 equations.

    ``compatible_bases=True`` records use of a simultaneous compatible choice;
    basis compatibility is a categorical convention, not decidable from N.
    No additional cross-system equality of x or y coordinates is imposed.
    """
    _require(isinstance(systems,(list,tuple)) and systems, 'At least one local system is required.')
    first = systems[0]
    signatures = [(s['fusion_tensor_sha256'],s['dimensions']['mode'],s['dimensions']['values']) for s in systems]
    _require(all(v == signatures[0] for v in signatures), 'Local systems must use the same fusion ring and dimension choice.')
    _require(not couplings or compatible_bases is True, 'Coupling requires compatible_bases=True for a simultaneous choice of vertex bases.')
    variables, metadata, centers = [], {}, {}
    for system in systems:
        _require(len(system['localizations']) == 1, 'Pass individual local systems to combine_localizations.')
        local = system['localizations'][0]
        _require(local['center'] not in centers, 'Distinguished centers must be different.')
        centers[local['center']] = (system,local)
        for name in system['variables']:
            if name in metadata:
                _require(metadata[name] == system['variable_metadata'][name] and metadata[name]['kind'] in ('categorical_dimension','inverse_dimension'),
                         'Namespaces collide between local systems.')
            else:
                variables.append(name);metadata[name] = system['variable_metadata'][name]
    space, eq, locals_out = _Space(variables), _Equations(variables), []
    for system in systems:
        old_names = system['variables']
        for row in system['equations']:
            for origin in row['origins']:
                eq.add(_lift(_decode(row['polynomial']),old_names,variables),origin)
        for origin in system['identically_zero_instances']:
            eq.add({},origin)
        local = dict(system['localizations'][0])
        local['coordinates'] = {kind:{key:_encode(_lift(_decode(value),old_names,variables)) for key,value in rows.items()}
                                for kind,rows in local['coordinates'].items()}
        locals_out.append(local)
    dimensions = [_lift(_decode(p),first['variables'],variables) for p in first['dimensions']['polynomials']]
    by_center = {local['center']:local for local in locals_out}
    coupling_records = []
    for pair in couplings:
        _require(isinstance(pair,(list,tuple)) and len(pair) == 2, 'Each ordered coupling is [k,l].')
        k,l = pair
        _require(type(k) is int and type(l) is int and k != l and k in centers and l in centers,
                 'Coupling centers must be distinct generated centers.')
        K,L = by_center[k],by_center[l]
        _require(l in K['selected_support'] and l in L['selected_support'],
                 'The second coupling center l must belong to both selected supports.')
        def diagonal(local,i,b):
            key = (i,b) if local['mode'] == 'corollary' else tuple(sorted((i,b,b)))
            return _decode(local['coordinates']['x'][','.join(map(str,key))])
        terms = [multiply(dimensions[i],multiply(_decode(L['coordinates']['y'][f'{i},{l}']),diagonal(K,i,l)))
                 for i in sorted(set(K['support']) & set(L['support']))]
        equation = space.difference(space.sum(terms),diagonal(K,l,l))
        origin = {'family':'compatibility','k':k,'l':l,'ordered':True,'compatible_bases':True}
        eq.add(equation,origin)
        coupling_records.append(dict(origin,polynomial=_encode(_normalize(equation)),text=_text(_normalize(equation),variables)))
    result = dict(first)
    result.update(variables=variables,variable_metadata=metadata,localizations=locals_out,equations=eq.rows,
                  identically_zero_instances=eq.zeros,couplings=coupling_records,
                  counts={'variables':len(variables),'equations':len(eq.rows),'instances':eq.instances,'identically_zero_instances':len(eq.zeros)},
                  basis_convention=('Coupled systems use compatible vertex bases in one category; this declaration adds only the stated ordered coupling equations.' if couplings else first['basis_convention']))
    result['dimensions'] = dict(first['dimensions'],polynomials=list(map(_encode,dimensions)),
        inverse_polynomials=[_encode(_lift(_decode(p),first['variables'],variables)) for p in first['dimensions']['inverse_polynomials']])
    return result


def generate_from_spec(fusion_matrices, spec):
    """Build one or several systems from a JSON-compatible specification."""
    _require(isinstance(spec,dict), 'localization specification must be an object.')
    allowed = {'mode','categorical_dimensions','centers','couplings','compatible_bases','ground_characteristic'}
    _require(set(spec) <= allowed, 'Unknown localization specification fields: '+', '.join(sorted(set(spec)-allowed)))
    _require(type(spec.get('ground_characteristic',0)) is int and spec.get('ground_characteristic',0) == 0,
             'This implementation generates characteristic-zero equations only.')
    centers = spec.get('centers')
    _require(isinstance(centers,list) and centers, 'specification.centers must be a nonempty list.')
    systems = []
    for center in centers:
        _require(isinstance(center,dict) and 'k' in center, 'Each center must specify k.')
        _require(set(center) <= {'k','subset','odd_witness','namespace'}, 'Unknown center specification field.')
        systems.append(generate_localization(fusion_matrices,k=center['k'],subset=center.get('subset'),
            categorical_dimensions=spec.get('categorical_dimensions'),mode=spec.get('mode','corollary'),
            odd_witness=center.get('odd_witness'),namespace=center.get('namespace')))
    return combine_localizations(systems,spec.get('couplings',()),compatible_bases=spec.get('compatible_bases',False))


def to_singular(system):
    """Export the generated exact ideal for optional computation in Singular."""
    names = system['variables']
    _require(all(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',name) for name in names), 'Invalid exported variable name.')
    equations = [_text(_decode(row['polynomial']),names) for row in system['equations']]
    if not names:
        names, equations = ['dummy'], equations+['dummy']
    return ('// Necessary localization equations over Q; roots may lie in its algebraic closure.\n'
            '// A non-unit ideal does not establish categorifiability.\n'
            'ring r=0,('+','.join(names)+'),dp;\n'
            'ideal I='+(',\n'.join(equations) or '0')+';\n'
            'ideal G=std(I);\nG;\nquit;\n')
