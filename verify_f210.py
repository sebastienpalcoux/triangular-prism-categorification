"""Verify the characteristic-zero obstruction using only Python's standard library.

The large certificate is an integer matrix T with E*T=G. No Groebner-basis
search or trust in the program that generated T is required for verification.
"""
from pathlib import Path
from fractions import Fraction
from functools import reduce
from math import gcd,lcm
import gzip,json,time
from polynomials import (parse,add_into,multiply,substitute,order,
                         standard_monomials,reduce_mod,determinant_mod)
ROOT=Path(__file__).resolve().parent


def verify_elimination(data):
    names=data['variables'];old=data['original_variables']
    assert names==['U1','U2','V1','V3','V4','V5','V6']
    assert old==['U0','U1','U2','V0','V1','V2','V3','V4','V5','V6']
    images=[parse(s,names,True) for s in [
        '4/5-7/5*U1-7/5*U2','U1','U2',
        '-V1-7/5*V3-7/5*V5-1/5','V1',
        '-V1-7/5*V4-7/5*V6-1/5','V3','V4','V5','V6']]
    reduced=[]
    for text in data['original_equations']:
        f=substitute(parse(text,old),images)
        if not f:continue
        denominator=lcm(*(Fraction(c).denominator for c in f.values()))
        reduced.append({e:int(c*denominator) for e,c in f.items()})
    expected=[parse(s,names) for s in data['reduced_equations']]
    assert reduced==expected, 'Linear elimination does not match the equations'
    assert len(reduced)==9
    return reduced


def verify_membership(data,equations,certificate_directory=None):
    certificate_directory=certificate_directory or ROOT/'certificates'
    names=data['variables']
    G=[parse(s,names) for s in (certificate_directory/'f210-basis.txt').read_text().strip().split(',')]
    with gzip.open(certificate_directory/'f210-lift.txt.gz','rt') as stream:
        entries=stream.read().strip().split(',')
    assert len(G)==28 and len(entries)==len(equations)*len(G)
    for j,g in enumerate(G):
        combination={}
        for i,e in enumerate(equations):
            q=parse(entries[i*len(G)+j],names)
            add_into(combination,multiply(q,e))
        assert combination==g, f'Integer ideal-membership identity {j+1} failed'
    # Remove content after the exact integer identity check. Content is nonzero
    # over Q; no claim that these divisions remain valid in characteristic p.
    return [{e:c//reduce(gcd,g.values()) for e,c in g.items()} for g in G]


def verify_unit(basis,prime=101):
    monomials=standard_monomials(basis);assert len(monomials)==14
    n=len(monomials);indices={b:i for i,b in enumerate(monomials)}
    reduced_basis=[]
    for g in basis:
        h=max(g,key=order)
        assert g[h]%prime!=0, 'A leading coefficient is not a unit modulo p'
        reduced={e:c%prime for e,c in g.items() if c%prime}
        assert max(reduced,key=order)==h
        reduced_basis.append(reduced)
    # Polynomial division uses only inverses of these leading coefficients.
    # It therefore takes place over Z_(p), and commutes with reduction modulo p.
    # Each matrix below is the reduction of an exact rational multiplication
    # matrix on the indicated spanning monomials, not a modular infeasibility test.
    matrices={}
    for variable in (0,1,2,3,5):
        matrix=[[0]*n for _ in range(n)]
        for j,b in enumerate(monomials):
            e=list(b);e[variable]+=1
            remainder=reduce_mod({tuple(e):1},reduced_basis,prime)
            for f,c in remainder.items():matrix[indices[f]][j]=c
        matrices[variable]=matrix
    identity=[[int(i==j) for j in range(n)] for i in range(n)]
    def linear_combination(terms):
        return [[sum(c*m[i][j] for c,m in terms)%prime for j in range(n)]for i in range(n)]
    inv5=pow(5,-1,prime)
    U0=linear_combination([(4*inv5,identity),(-7*inv5,matrices[0]),(-7*inv5,matrices[1])])
    V0=linear_combination([(-inv5,identity),(-1,matrices[2]),(-7*inv5,matrices[3]),(-7*inv5,matrices[5])])
    terms=[(5,U0,V0),(7,matrices[0],matrices[3]),(7,matrices[1],matrices[5]),(-5,U0,identity),(1,identity,identity)]
    M=[[sum(c*A[i][k]*B[j][l] for c,A,B in terms)%prime
        for k in range(n) for l in range(n)]for i in range(n)for j in range(n)]
    determinant=determinant_mod(M,prime)
    assert determinant==85, f'Unexpected determinant residue: {determinant}'
    return {'spanning_monomials':n,'matrix_size':n*n,'prime':prime,'determinant_residue':determinant}


def verify():
    start=time.monotonic()
    data=json.loads((ROOT/'data/f210-equations.json').read_text())
    equations=verify_elimination(data)
    basis=verify_membership(data,equations)
    result=verify_unit(basis,data['prime'])
    result.update(membership_identities=len(basis),seconds=round(time.monotonic()-start,3))
    print('F210 characteristic zero: PASS',json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Run without -O: mathematical checks use assertions')
    verify()
