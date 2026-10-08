#!/usr/bin/env python3
"""Verify the paper's frozen arithmetic certificates; Python 3.10+, no packages."""
import hashlib, json, subprocess, sys, time
from pathlib import Path
from spectrum import scan
from verify_localization import verify as verify_localization
from verify_f210 import verify as verify_f210
from verify_f210_characters import verify as verify_characters
ROOT=Path(__file__).resolve().parent


def ring_checks(N, dimensions=None):
    n=len(N);I=range(n)
    assert all(len(A)==n and all(len(row)==n for row in A) for A in N)
    assert all(type(c) is int and c>=0 for A in N for row in A for c in row)
    assert all(N[0][i][j]==N[i][0][j]==int(i==j) for i in I for j in I)
    star=[]
    for i in I:
        duals=[j for j in I if N[i][j][0]]
        assert len(duals)==1 and N[i][duals[0]][0]==1
        star.append(duals[0])
    for i in I:
        assert star[star[i]]==i
        for j in I:
            if dimensions is not None:
                assert sum(N[i][j][k]*dimensions[k] for k in I)==dimensions[i]*dimensions[j]
            for k in I:
                assert N[i][j][k]==N[j][i][k]
                assert N[i][j][k]==N[star[i]][k][j]==N[k][star[j]][i]
                assert N[i][j][k]==N[star[j]][star[i]][star[k]]
                for l in I:
                    assert sum(N[i][j][s]*N[s][k][l] for s in I)==sum(N[j][k][s]*N[i][s][l] for s in I)
    return star


def ring_product(a,b,N):
    n=len(N)
    return [sum(a[i]*b[j]*N[i][j][k] for i in range(n)for j in range(n))for k in range(n)]


def matmul(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]


def characteristic_polynomial(A):
    n=len(A);B=[[int(i==j)for j in range(n)]for i in range(n)];coefficients=[1]
    for k in range(1,n+1):
        B=matmul(A,B);trace=sum(B[i][i]for i in range(n));assert trace%k==0
        c=-trace//k;coefficients.append(c)
        for i in range(n):B[i][i]+=c
    assert not any(any(row)for row in B)
    return coefficients


def positive_characteristic(N):
    n=len(N);dimensions=[1,5,5,5,6,7,7]
    C=[[0]*n for _ in range(n)]
    for A in N:
        AAT=matmul(A,list(map(list,zip(*A))))
        for i in range(n):
            for j in range(n):C[i][j]+=AAT[i][j]
    expected=[1]
    for root in [5,5,6,7,7,7,210]:
        next_coefficients=[0]*(len(expected)+1)
        for i,c in enumerate(expected):next_coefficients[i]+=c;next_coefficients[i+1]-=root*c
        expected=next_coefficients
    assert characteristic_polynomial(C)==expected
    product=[1]+[0]*(n-1)
    for i in range(1,n):product=ring_product(product,[int(i==j)for j in range(n)],N)
    assert product==[175*d for d in dimensions]
    assert ring_product(dimensions,dimensions,N)==[210*d for d in dimensions]
    assert ring_product(product,product,N)==[36750*d for d in product]
    assert all(36750%p==0 for p in (2,3,5,7))
    print('F210 codegrees and positive-characteristic identities: PASS')


def f660_witness(N,star):
    i1,i2,i3,i4,i5,i6,i7,i8,i9=(1,3,4,1,1,3,4,2,2);I=range(len(N))
    dot=lambda a,b,c,d:sum(N[a][b][k]*N[c][d][k]for k in I)
    six=[N[i4][i1][i6],N[i5][i4][i2],N[i5][i6][i3],N[i7][i9][i1],N[i2][i7][i8],N[i8][i9][i3]]
    spectrum=[N[i4][i7][k]*N[star[i5]][i8][k]*N[i6][star[i9]][k]for k in I]
    first=[dot(i5,i6,i2,i1),dot(i5,i4,i3,star[i1]),dot(i2,star[i4],i3,star[i6])]
    second=[dot(i2,i1,i8,i9),dot(i2,i7,i3,star[i9]),dot(i8,star[i7],i3,star[i1])]
    assert all(six) and sum(spectrum)==0 and N[i2][i1][i3]==1 and 1 in first and 1 in second
    print('F660 zero-spectrum witness: PASS',json.dumps({'six_coefficients':six,'spectrum':spectrum,'first_sums':first,'second_sums':second}))


def main():
    start=time.monotonic()
    manifest=json.loads((ROOT/'data/manifest-sha256.json').read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest, f'Input checksum mismatch: {name}'
    print('Frozen input and certificate SHA256 hashes: PASS')
    rings={name:json.loads((ROOT/'data'/f'{name}.json').read_text())for name in ('f210','f660','rank6_zero_only','rank6_one_only')}
    dimensions={'f210':[1,5,5,5,6,7,7],'f660':[1,5,5,10,10,11,12,12]}
    stars={name:ring_checks(N,dimensions.get(name))for name,N in rings.items()}
    print('Four fusion tensors and stated dimensions: PASS')
    verify_localization(json.loads((ROOT/'data/f210-equations.json').read_text()),rings['f210'])
    positive_characteristic(rings['f210'])
    verify_characters(rings['f210'])
    f660_witness(rings['f660'],stars['f660'])
    for name,expected in [('rank6_zero_only',{0:12,1:0}),('rank6_one_only',{0:0,1:96})]:
        result=scan(rings[name]);assert result['witness_count']==expected
        print(name+' exhaustive zero/one-spectrum counts: PASS',expected)
    verify_f210()
    subprocess.run([sys.executable,str(ROOT/'census/verify_census.py')],check=True)
    subprocess.run([sys.executable,'-B','-S','-m','unittest','discover','-s',str(ROOT/'tests')],check=True)
    print(f'ALL CHECKS PASS ({time.monotonic()-start:.3f} seconds)')


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Run without -O: mathematical checks use assertions')
    main()
