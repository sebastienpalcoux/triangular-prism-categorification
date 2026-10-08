"""Generate all 27 corollary instances for each localization, using exact fractions."""
from fractions import Fraction as Q
from itertools import product
from polynomials import parse,add_into,multiply,substitute,order


def verify(data,N):
    names=data['original_variables'];constant=lambda c:{(0,)*10:Q(c)} if c else {}
    def scaled(text):return parse(text,names,True)
    def total(terms):
        result={}
        for coefficient,term in terms:add_into(result,term,coefficient)
        return result
    def canonical(f):
        if not f:return ()
        coefficient=f[max(f,key=order)]
        return tuple(sorted((e,Q(c)/coefficient)for e,c in f.items()))
    equations=[parse(s,names)for s in data['original_equations']]
    targets={canonical(p)for p in equations}
    assert len(targets)==12
    Y=[[constant(Q(1,5))for _ in range(3)]for _ in range(5)]
    for i,j,v in [(1,1,0),(1,2,1),(2,1,1),(2,2,2),(3,1,3),(3,2,4),(4,1,5),(4,2,6)]:Y[i][j]=scaled(f'1/5*V{v}')
    def instances(distinguished):
        other=3-distinguished;X=[[{}for _ in range(3)]for _ in range(5)]
        X[0][0]=constant(Q(1,5));X[0][other]=constant(Q(1,25))
        for i in range(5):X[i][distinguished]=multiply(Y[i][distinguished],Y[i][distinguished])
        for i,u in [(other,0),(3,1),(4,2)]:X[i][other]=scaled(f'1/25*U{u}')
        d=[1,5,5,7,7];result=[]
        for a,b in product(range(3),repeat=2):
            first=total((d[b]*d[i],multiply(Y[i][a],Y[i][b]))for i in range(5));add_into(first,constant(int(a==b)),-1)
            second=total((d[i],multiply(Y[i][a],multiply(Y[i][b],Y[i][b])))for i in range(5));add_into(second,X[a][b],-1)
            third=total((d[i],multiply(Y[i][a],X[i][b]))for i in range(5));add_into(third,multiply(Y[a][b],Y[a][b]),-1)
            result.extend([canonical(first),canonical(second),canonical(third)])
        return set(result)-{()}
    assert instances(1)==targets
    permutation=[scaled(v)for v in ['U0','U1','U2','V2','V1','V0','V4','V3','V6','V5']]
    assert instances(2)=={canonical(substitute(p,permutation))for p in equations}
    assert [i for i,c in enumerate(N[1][1])if c]==[0,1,3,5,6]
    assert [i for i,c in enumerate(N[3][3])if c]==[0,2,3,5,6]
    assert N[3][3][1]==N[2][2][3]==0
    print('Localization: 27 instances -> 12 equations, second-system permutation, supports: PASS')
