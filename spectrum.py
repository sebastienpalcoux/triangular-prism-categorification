"""Exhaustive, exact, standard-library scan of the printed spectrum criteria.

This independent audit follows the six required nonzero coefficients to
enumerate admissible labels. Every omitted tuple fails an explicit hypothesis.
All output witness indices are zero-based, in order i0,i1,...,i9 (one) or
i1,...,i9 (zero). No symmetry quotient is used.
"""
import json, time
from itertools import product
from functools import lru_cache
from pathlib import Path

def scan(N):
    n=len(N);I=range(n)
    dual=[]
    for i in I:
        candidates=[j for j in I if N[i][j][0]]
        assert len(candidates)==1 and N[i][candidates[0]][0]==1
        dual.append(candidates[0])
    support=[[tuple(k for k in I if N[i][j][k]) for j in I] for i in I]
    dots={}
    for a,b,c,d in product(I,repeat=4):
        dots[a,b,c,d]=sum(N[a][b][k]*N[c][d][k] for k in I)
    def one_of(*keys):return any(dots[key]==1 for key in keys)
    @lru_cache(None)
    def spectrum(i4,i7,i5s,i8,i6,i9s):
        vals=[N[i4][i7][k]*N[i5s][i8][k]*N[i6][i9s][k] for k in I]
        total=sum(vals)
        return total,vals.index(1) if total==1 else None
    witness={0:None,1:None};passes={0:0,1:0}
    for i4,i5,i7,i9 in product(I,repeat=4):
      for i1 in support[i7][i9]:
       for i2 in support[i5][i4]:
        for i6 in support[i4][i1]:
         for i8 in support[i2][i7]:
          total,i0=spectrum(i4,i7,dual[i5],i8,i6,dual[i9])
          if total>1:continue
          if total==1:
              if not one_of((i5,i0,i2,i7),(i5,i4,i8,dual[i7]),(i2,dual[i4],i8,dual[i0])):continue
              if not one_of((i4,i1,i0,i9),(i4,i7,i6,dual[i9]),(i0,dual[i7],i6,dual[i1])):continue
          for i3 in support[i5][i6]:
           if not N[i8][i9][i3] or N[i2][i1][i3] != 1-total:continue
           if total==0:
               if not one_of((i5,i6,i2,i1),(i5,i4,i3,dual[i1]),(i2,dual[i4],i3,dual[i6])):continue
               if not one_of((i2,i1,i8,i9),(i2,i7,i3,dual[i9]),(i8,dual[i7],i3,dual[i1])):continue
           else:
               if not one_of((i5,i6,i8,i9),(i5,i0,i3,dual[i9]),(i8,dual[i0],i3,dual[i6])):continue
           passes[total]+=1
           if witness[total] is None:
               witness[total]=([i0] if total else [])+[i1,i2,i3,i4,i5,i6,i7,i8,i9]
    return {'rank':n,'dual':dual,'witness_count':passes,'first_witness':witness,
            'spectrum_cache_entries':spectrum.cache_info().currsize,
            'exhaustive':True}

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    output={}
    for name in ['rank6_one_only','rank6_zero_only']:
        start=time.monotonic()
        result=scan(json.loads((root/f'data/{name}.json').read_text()))
        result['seconds']=round(time.monotonic()-start,3)
        output[name]=result
        print(name,json.dumps(result),flush=True)

