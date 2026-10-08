"""Small exact arithmetic checker, not a computer algebra search engine."""
from fractions import Fraction
import re


def parse(text, names, rational=False):
    """Read the deliberately restricted certificate grammar; never use eval."""
    text=''.join(text.split())
    if text=='0': return {}
    result={}
    for term in re.findall(r'[+-]?[^+-]+',text):
        sign=-1 if term.startswith('-') else 1
        term=term.lstrip('+-'); coefficient=sign; exponents=[0]*len(names)
        for factor in term.split('*'):
            if re.fullmatch(r'\d+(?:/\d+)?',factor):
                coefficient*=Fraction(factor) if rational else int(factor)
            else:
                match=re.fullmatch(r'([UV]\d)(?:\^(\d+))?',factor)
                if match is None or match[1] not in names: raise ValueError('Invalid monomial')
                exponents[names.index(match[1])]+=int(match[2] or 1)
        e=tuple(exponents);result[e]=result.get(e,0)+coefficient
    return {e:c for e,c in result.items() if c}


def add_into(target, polynomial, scalar=1):
    for e,c in polynomial.items():
        value=target.get(e,0)+scalar*c
        if value:target[e]=value
        else:target.pop(e,None)


def multiply(a,b):
    result={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f));result[g]=result.get(g,0)+c*d
    return {e:c for e,c in result.items() if c}


def substitute(polynomial, images):
    n=len(next(iter(images[0])))
    result={}
    for exponents,c in polynomial.items():
        term={(0,)*n:c}
        for i,power in enumerate(exponents):
            for _ in range(power):term=multiply(term,images[i])
        add_into(result,term)
    return result


def order(exponents):
    """Degree reverse lexicographic order in the declared variable order."""
    return (sum(exponents),tuple(-x for x in exponents[::-1]))


def divides(a,b):return all(x<=y for x,y in zip(a,b))


def standard_monomials(basis):
    n=len(next(iter(basis[0])));heads=[max(g,key=order) for g in basis]
    # These explicit bounds guarantee that the following breadth-first search terminates.
    for i in range(n):
        assert any(h[i]>0 and sum(h)==h[i] for h in heads), 'Missing pure-power bound'
    found={(0,)*n};pending=list(found)
    while pending:
        e=pending.pop()
        for i in range(n):
            f=list(e);f[i]+=1;f=tuple(f)
            if f not in found and not any(divides(h,f) for h in heads):
                found.add(f);pending.append(f)
    return sorted(found,key=order)


def reduce_mod(polynomial,basis,p):
    f=polynomial.copy();remainder={}
    heads=[max(g,key=order) for g in basis]
    inverses=[pow(g[h],-1,p) for g,h in zip(basis,heads)]
    while f:
        m=max(f,key=order);c=f[m]
        for g,h,inverse in zip(basis,heads,inverses):
            if divides(h,m):
                shift=tuple(x-y for x,y in zip(m,h));q=c*inverse%p
                for e,t in g.items():
                    a=tuple(x+y for x,y in zip(e,shift));value=(f.get(a,0)-q*t)%p
                    if value:f[a]=value
                    else:f.pop(a,None)
                break
        else:remainder[m]=f.pop(m)
    return remainder


def determinant_mod(matrix,p):
    a=[row[:] for row in matrix];n=len(a);det=1
    for j in range(n):
        k=next((k for k in range(j,n) if a[k][j]),None)
        if k is None:return 0
        if k!=j:a[j],a[k]=a[k],a[j];det=-det
        pivot=a[j][j];det=det*pivot%p;inverse=pow(pivot,-1,p)
        for k in range(j+1,n):
            if a[k][j]:
                factor=a[k][j]*inverse%p
                for l in range(j,n):a[k][l]=(a[k][l]-factor*a[j][l])%p
    return det%p
