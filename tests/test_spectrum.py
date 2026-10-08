"""Independent regression/oracle tests for the reusable spectrum interface.

Run with python -S tests/test_spectrum.py.
This intentionally does not import the optimized scanner's witness evaluator
when constructing the dense oracle.
"""
from pathlib import Path
import copy
import itertools
import json
import subprocess
import sys
import tempfile
import unittest
REPO=Path(__file__).resolve().parents[1]
ROOT=REPO.parent
sys.path.insert(0,str(REPO))
from fusion_ring import FusionRingError,load_fusion_ring,validate_fusion_ring
from spectrum import scan,evaluate_witness


def cyclic(n):
    return [[[int((i+j)%n==k) for k in range(n)]for j in range(n)]for i in range(n)]


def symmetric_group_three():
    elements=list(itertools.permutations(range(3)))
    return [[[int(tuple(a[b[x]] for x in range(3))==c) for c in elements]
             for b in elements]for a in elements]


REP_S3=[[[1,0,0],[0,1,0],[0,0,1]],[[0,1,0],[1,0,0],[0,0,1]],[[0,0,1],[0,0,1],[1,1,1]]]
FIB=[[[1,0],[0,1]],[[0,1],[1,1]]]
ISING=[[[1,0,0],[0,1,0],[0,0,1]],[[0,1,0],[1,0,0],[0,0,1]],[[0,0,1],[0,0,1],[1,1,0]]]


def dense_oracle(N):
    n=len(N);I=range(n)
    star=[next(j for j in I if N[i][j][0]==1)for i in I]
    result={0:set(),1:set()}
    def scalar(a,b,c,d):
        return sum(N[a][b][t]*N[c][d][t]for t in I)
    for labels in itertools.product(I,repeat=9):
        a,b,c,d,e,f,g,h,j=labels
        if not (N[d][a][f] and N[e][d][b] and N[e][f][c] and N[g][j][a] and N[b][g][h] and N[h][j][c]):
            continue
        spectrum=sum(N[d][g][t]*N[star[e]][h][t]*N[f][star[j]][t]for t in I)
        if spectrum==0 and N[b][a][c]==1:
            if (scalar(e,f,b,a)==1 or scalar(e,d,c,star[a])==1 or scalar(b,star[d],c,star[f])==1) and (scalar(b,a,h,j)==1 or scalar(b,g,c,star[j])==1 or scalar(h,star[g],c,star[a])==1):
                result[0].add(labels)
        if spectrum==1 and N[b][a][c]==0:
            for q in I:
                if not N[d][g][q]==N[star[e]][h][q]==N[f][star[j]][q]==1:
                    continue
                if not (scalar(e,q,b,g)==1 or scalar(e,d,h,star[g])==1 or scalar(b,star[d],h,star[q])==1):
                    continue
                if not (scalar(e,f,h,j)==1 or scalar(e,q,c,star[j])==1 or scalar(h,star[q],c,star[f])==1):
                    continue
                if not (scalar(d,a,q,j)==1 or scalar(d,g,f,star[j])==1 or scalar(q,star[g],f,star[a])==1):
                    continue
                result[1].add((q,)+labels)
    return result


class GenericSpectrumTests(unittest.TestCase):
    def test_known_categorifiable_positive_controls(self):
        for name,N,d in [('unit',cyclic(1),[1]),('Z2',cyclic(2),[1,1]),
                         ('Z3',cyclic(3),[1,1,1]),('Fibonacci',FIB,None),
                         ('Ising',ISING,None),('Rep(S3)',REP_S3,[1,1,2]),
                         ('Z[S3]',symmetric_group_three(),[1]*6)]:
            with self.subTest(name=name):
                report=validate_fusion_ring(N,d)
                self.assertEqual(scan(N)['witness_count'],{0:0,1:0})
                self.assertEqual(report['commutative'],name!='Z[S3]')
        self.assertEqual(validate_fusion_ring(cyclic(3))['dual'],[0,2,1])
        self.assertNotEqual(validate_fusion_ring(symmetric_group_three())['dual'],list(range(6)))

    def test_dense_rank_two_and_three_oracle(self):
        for N in (cyclic(2),FIB,cyclic(3),REP_S3,ISING):
            with self.subTest(rank=len(N),tensor=N):
                dense=dense_oracle(N)
                fast=scan(N,mode='all')
                found={k:{tuple(w['indices'])for w in fast['witness_details'][k]}for k in (0,1)}
                self.assertEqual(found,dense)
                self.assertTrue(fast['exhaustive'])

    def test_all_frozen_rank_six_witnesses_and_modes(self):
        for filename,expected in [('rank6_zero_only',{0:12,1:0}),('rank6_one_only',{0:0,1:96})]:
            N=json.loads((REPO/'data'/f'{filename}.json').read_text())
            dense=dense_oracle(N)
            fast=scan(N,mode='all')
            self.assertEqual({k:len(v)for k,v in dense.items()},expected)
            self.assertEqual(fast['witness_count'],expected)
            self.assertEqual({k:{tuple(w['indices'])for w in fast['witness_details'][k]}for k in (0,1)},dense)
            for criterion,kind in [('zero',0),('one',1)]:
                selected=scan(N,criterion,mode='count')
                self.assertEqual(selected['witness_count'],{kind:expected[kind]})
                first=scan(N,criterion,mode='first')
                self.assertEqual(first['witness_count'][kind],int(expected[kind]>0))
                self.assertEqual(first['exhaustive'],not bool(expected[kind]))
                self.assertEqual(first['witness_count_is_complete'],first['exhaustive'])
                if expected[kind]:
                    witness=first['first_witness'][kind]
                    self.assertIn(tuple(witness),dense[kind])
                    self.assertTrue(evaluate_witness(N,witness,criterion)['passes'])
                    bad=witness.copy();bad[-1]=0
                    self.assertFalse(evaluate_witness(N,bad,criterion)['passes'])
                    if kind:
                        wrong_i0=witness.copy();wrong_i0[0]=1
                        self.assertFalse(evaluate_witness(N,wrong_i0,'one')['passes'])
            # Both-mode cannot certify a complete count when it stopped early;
            # these inputs have only one kind, so it must finish the full scan.
            self.assertTrue(scan(N,mode='first')['exhaustive'])

    def test_malformed_rules_rejected(self):
        candidates=[None,[],[[[1,0]]]]
        for value in (-1,True,1.0):
            N=copy.deepcopy(FIB);N[1][1][1]=value;candidates.append(N)
        N=copy.deepcopy(FIB);N[0][0][0]=0;candidates.append(N)
        N=copy.deepcopy(FIB);N[1][1][0]=0;candidates.append(N)
        N=copy.deepcopy(FIB);N[1][1][0]=2;candidates.append(N)
        # This tensor has unit, self-duality, nonnegative coefficients, and
        # reciprocity, but (a*a)*b != a*(a*b). It must fail associativity.
        N=[[[1,0,0],[0,1,0],[0,0,1]],[[0,1,0],[1,0,0],[0,0,0]],[[0,0,1],[0,0,0],[1,0,0]]]
        with self.assertRaisesRegex(FusionRingError,'Associativity'):
            validate_fusion_ring(N)
        for N in candidates:
            with self.subTest(N=N),self.assertRaises(FusionRingError):
                validate_fusion_ring(N)
        for d in ([1,1],[1,1.618],[1,True],['1','1/0'],[1,-1]):
            with self.subTest(d=d),self.assertRaises(FusionRingError):
                validate_fusion_ring(FIB,d)
        self.assertEqual(validate_fusion_ring(REP_S3,['1','1/1','2'])['dimensions'],[1,1,2])

    def test_cli_raw_and_metadata_inputs(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'input.json';output=root/'output.json'
            for data in (REP_S3,{'fusion_matrices':REP_S3,'labels':['1','sign','standard'],
                               'dimensions':[1,1,2],'name':'Rep(S3)','unit_index':0}):
                source.write_text(json.dumps(data))
                N,metadata=load_fusion_ring(source)
                self.assertEqual(N,REP_S3)
                done=subprocess.run([sys.executable,'-S',str(REPO/'criteria.py'),'spectrum',str(source),
                                     '--mode','all','--output',str(output)],cwd='/tmp',text=True,capture_output=True)
                self.assertEqual(done.returncode,0,done.stderr)
                result=json.loads(output.read_text())
                self.assertEqual(result['witness_count'],{'0':0,'1':0})
                self.assertIn('does not establish',result['interpretation'])
            for data in ({'fusion_matrices':REP_S3,'unit_index':1},
                         {'fusion_matrices':REP_S3,'labels':['x','x','y']},
                         {'fusion_matrices':REP_S3,'dimensions':[1,1,3]},
                         None):
                source.write_text(json.dumps(data))
                done=subprocess.run([sys.executable,str(REPO/'criteria.py'),'spectrum',str(source)],text=True,capture_output=True)
                self.assertEqual(done.returncode,2,done.stdout+done.stderr)
                self.assertNotIn('Traceback',done.stderr)
            source.write_text('{"fusion_matrices": NaN}')
            with self.assertRaises(FusionRingError):load_fusion_ring(source)


if __name__=='__main__':unittest.main(verbosity=2)
