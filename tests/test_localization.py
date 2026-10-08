"""Regression against the paper and exact positive controls for the general API."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from localization import (generate_localization,generate_from_spec,combine_localizations,
                          to_singular,_decode,_normalize)
from polynomials import parse,substitute


def example(name):return json.loads((ROOT/'examples'/f'localization_{name}.json').read_text())
def canonical(poly):return tuple(sorted(_normalize(poly).items()))
def coordinate(n,i,c=1):
    e=[0]*n;e[i]=1;return {tuple(e):Q(c)}
def evaluate(poly,values):
    return sum(c*__import__('math').prod(value**power for value,power in zip(values,e))for e,c in poly.items())


class GenericLocalizationTests(unittest.TestCase):
    def test_f210_exactly_recovers_both_twelve_equation_systems(self):
        data=example('f210');N=data['fusion_matrices'];spec=data['localization']
        reference=json.loads((ROOT/'data/f210-equations.json').read_text())
        expected={canonical(parse(s,reference['original_variables']))for s in reference['original_equations']}
        maps=[
            ({(3,3):0,(5,3):1,(6,3):2},{(1,1):3,(1,3):4,(3,3):5,(1,5):6,(3,5):7,(1,6):8,(3,6):9}),
            ({(2,2):0,(5,2):1,(6,2):2},{(3,3):3,(2,3):4,(2,2):5,(3,5):6,(2,5):7,(3,6):8,(2,6):9})]
        for center,(xmap,ymap) in zip(spec['centers'],maps):
            system=generate_localization(N,categorical_dimensions=spec['categorical_dimensions'],**center)
            self.assertEqual(system['counts']['variables'],10)
            self.assertEqual(system['counts']['equations'],12)
            self.assertEqual(system['counts']['instances'],27)
            images=[]
            for name in system['variables']:
                meta=system['variable_metadata'][name];key=tuple(meta['indices'])
                images.append(coordinate(10,(xmap if meta['kind']=='x'else ymap)[key],Q(1,25 if meta['kind']=='x' else 5)))
            actual={canonical(substitute(_decode(row['polynomial']),images))for row in system['equations']}
            self.assertEqual(actual,expected)

    def test_f210_coupling_is_the_exact_mixed_polynomial(self):
        data=example('f210');system=generate_from_spec(data['fusion_matrices'],data['localization'])
        self.assertEqual(system['counts']['variables'],20)
        self.assertEqual(system['counts']['equations'],25)
        self.assertEqual(len(system['couplings']),1)
        xm={1:{(3,3):0,(5,3):1,(6,3):2},3:{(2,2):10,(5,2):11,(6,2):12}}
        ym={1:{(1,1):3,(1,3):4,(3,3):5,(1,5):6,(3,5):7,(1,6):8,(3,6):9},
            3:{(2,2):13,(2,3):14,(3,3):15,(2,5):16,(3,5):17,(2,6):18,(3,6):19}}
        images=[]
        for name in system['variables']:
            meta=system['variable_metadata'][name];key=tuple(meta['indices']);center=meta['center']
            images.append(coordinate(20,(xm if meta['kind']=='x'else ym)[center][key],Q(1,25 if meta['kind']=='x' else 5)))
        actual=substitute(_decode(system['couplings'][0]['polynomial']),images)
        target={(0,)*20:1}
        for i,j,c in [(0,15,5),(1,17,7),(2,19,7)]:
            e=[0]*20;e[i]+=1;e[j]+=1;target[tuple(e)]=c
        e=[0]*20;e[0]=1;target[tuple(e)]=-5
        self.assertEqual(canonical(actual),canonical(target))

    def test_known_s3_category_solves_full_and_corollary_systems(self):
        data=example('s3');N=data['fusion_matrices'];Y=[[Q(1,2)]*3,[Q(1,2),Q(1,2),Q(-1,2)],[Q(1,2),Q(-1,2),Q(0)]]
        for mode,selected in [('full',None),('corollary',None),('corollary',[1])]:
            system=generate_localization(N,k=2,subset=selected,categorical_dimensions=[1,1,2],mode=mode)
            values=[]
            for name in system['variables']:
                meta=system['variable_metadata'][name]
                if meta['kind']=='y':i,j=meta['indices'];values.append(Y[i][j])
                else:
                    args=meta['indices'] if mode=='full'else[meta['indices'][0],meta['indices'][1],meta['indices'][1]]
                    values.append(sum(d*__import__('math').prod(Y[i][a]for a in args)for i,d in enumerate([1,1,2])))
            self.assertTrue(all(evaluate(_decode(row['polynomial']),values)==0 for row in system['equations']))
            self.assertIn('ring r=0,',to_singular(system))

    def test_symbolic_fibonacci_solution_in_exact_quadratic_field(self):
        data=example('fibonacci');system=generate_from_spec(data['fusion_matrices'],data['localization'])
        # Pairs (a,b) represent a+b*t modulo t^2-t-1.
        def mul(x,y):return(x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
        values=[]
        for name in system['variables']:
            kind=system['variable_metadata'][name]['kind']
            values.append({'categorical_dimension':(Q(0),Q(1)),'inverse_dimension':(Q(-1),Q(1)),'y':(Q(-2),Q(1))}[kind])
        for row in system['equations']:
            total=(Q(0),Q(0))
            for e,c in _decode(row['polynomial']).items():
                term=(c,Q(0))
                for power,value in zip(e,values):
                    for _ in range(power):term=mul(term,value)
                total=(total[0]+term[0],total[1]+term[1])
            self.assertEqual(total,(0,0))
        self.assertEqual(system['dimensions']['mode'],'symbolic')
        families={origin['family']for row in system['equations']for origin in row['origins']}
        self.assertIn('nonzero_dimension',families);self.assertIn('dimension_character',families)

    def test_coupled_s3_tensor_square_has_an_exact_common_solution(self):
        # This checks a nonzero ordered compatibility equation on an actual
        # tensor-product category, including shared symbolic dimensions.
        S=example('s3')['fusion_matrices'];d=[1,1,2]
        N=[[[S[a//3][b//3][c//3]*S[a%3][b%3][c%3]for c in range(9)]for b in range(9)]for a in range(9)]
        D=[d[a//3]*d[a%3]for a in range(9)]
        Y=[[Q(1,2)]*3,[Q(1,2),Q(1,2),Q(-1,2)],[Q(1,2),Q(-1,2),Q(0)]]
        X=[[sum(d[t]*Y[t][a]*Y[t][b]**2 for t in range(3))for b in range(3)]for a in range(3)]
        for dimensions,count in [(D,24),(None,87)]:
            system=generate_from_spec(N,{'mode':'corollary','categorical_dimensions':dimensions,
                'centers':[{'k':8,'subset':[0,6,8]},{'k':6,'subset':[0,3,6]}],
                'couplings':[[8,6]],'compatible_bases':True})
            self.assertEqual(system['counts']['equations'],count)
            self.assertTrue(system['couplings'][0]['polynomial'])
            values=[]
            for name in system['variables']:
                meta=system['variable_metadata'][name];kind=meta['kind'];indices=meta['indices']
                if kind=='categorical_dimension':value=Q(D[indices[0]])
                elif kind=='inverse_dimension':value=Q(1,D[indices[0]])
                else:
                    a,b=indices;table=X if kind=='x'else Y
                    value=table[a//3][b//3]
                    if meta['center']==8:value*=table[a%3][b%3]
                    else:self.assertEqual((a%3,b%3),(0,0))
                values.append(value)
            self.assertTrue(all(evaluate(_decode(row['polynomial']),values)==0 for row in system['equations']))
            self.assertEqual(evaluate(_decode(system['couplings'][0]['polynomial']),values),0)

    def test_proper_subset_uses_all_summands_and_need_not_contain_unit(self):
        data=example('f210');system=generate_localization(data['fusion_matrices'],k=1,subset=[3],categorical_dimensions=data['localization']['categorical_dimensions'])
        self.assertEqual(system['localizations'][0]['support'],[0,1,3,5,6])
        self.assertEqual(system['localizations'][0]['selected_support'],[3])
        self.assertIn('k1_y_3_5',system['variables']);self.assertIn('k1_y_3_6',system['variables'])
        self.assertIn('0,3',system['localizations'][0]['coordinates']['y'])
        self.assertTrue(any('k1_y_3_5' in row['text'] for row in system['equations']))

    def test_invalid_hypotheses_and_dimension_specializations(self):
        s3=example('s3')['fusion_matrices'];fib=example('fibonacci')['fusion_matrices']
        with self.assertRaisesRegex(ValueError,'dimension character'):generate_localization(s3,k=2,categorical_dimensions=[1,1,3])
        with self.assertRaisesRegex(ValueError,'floats'):generate_localization(s3,k=2,categorical_dimensions=[1,1,2.0])
        with self.assertRaisesRegex(ValueError,'nonzero'):generate_localization(s3,k=2,categorical_dimensions=[1,1,0])
        self.assertEqual(generate_localization(s3,k=2,categorical_dimensions=[Q(1),Q(1),Q(2)])['dimensions']['values'],['1','1','2'])
        with self.assertRaisesRegex(ValueError,'characteristic-zero'):generate_from_spec(s3,{'centers':[{'k':2}],'ground_characteristic':5})
        with self.assertRaisesRegex(ValueError,'Unknown'):generate_from_spec(s3,{'centers':[{'k':2}],'categorical_dimension':[1,1,2]})
        # A nonzero negative character is a legitimate specialization, even when
        # the resulting necessary system may have no solution.
        self.assertEqual(generate_localization(s3,k=2,categorical_dimensions=[1,1,-1])['dimensions']['values'],['1','1','-1'])
        with self.assertRaisesRegex(ValueError,'supplied odd'):generate_localization(s3,k=2,odd_witness=0)
        with self.assertRaisesRegex(ValueError,'duplicate'):generate_localization(s3,k=2,subset=[0,0])
        with self.assertRaisesRegex(ValueError,'Full theorem'):generate_localization(s3,k=2,subset=[1],mode='full')
        with self.assertRaisesRegex(ValueError,'support'):generate_localization(fib,k=0,subset=[1])
        c2=[[[1,0],[0,1]],[[0,1],[1,0]]]
        with self.assertRaisesRegex(ValueError,'odd-multiplicity'):generate_localization(c2,k=1)
        c3=[[[int(k==(i+j)%3)for k in range(3)]for j in range(3)]for i in range(3)]
        with self.assertRaisesRegex(ValueError,'self-dual'):generate_localization(c3,k=1)
        # Tambara-Yamagami fusion rules for C3: self-dual center, non-self-dual support.
        ty=[[[0]*4 for _ in range(4)]for _ in range(4)]
        for i in range(3):
            for j in range(3):ty[i][j][(i+j)%3]=1
            ty[i][3][3]=ty[3][i][3]=1
        ty[3][3]=[1,1,1,0]
        with self.assertRaisesRegex(ValueError,'summand'):generate_localization(ty,k=3)

    def test_ordered_couplings_require_both_domains_and_compatible_bases(self):
        data=example('f210');N=data['fusion_matrices'];d=data['localization']['categorical_dimensions']
        A=generate_localization(N,k=1,subset=[0,1,3],categorical_dimensions=d)
        B=generate_localization(N,k=3,subset=[0,2,3],categorical_dimensions=d)
        with self.assertRaisesRegex(ValueError,'compatible_bases'):combine_localizations([A,B],[[1,3]])
        with self.assertRaisesRegex(ValueError,'both selected'):combine_localizations([A,B],[[3,1]],compatible_bases=True)
        Bsmall=generate_localization(N,k=3,subset=[0,2],categorical_dimensions=d)
        with self.assertRaisesRegex(ValueError,'both selected'):combine_localizations([A,Bsmall],[[1,3]],compatible_bases=True)
        with self.assertRaisesRegex(ValueError,'different'):combine_localizations([A,A])

    def test_noncommutative_input_and_unit_center_are_supported(self):
        elements=list(permutations(range(3)));lookup={g:i for i,g in enumerate(elements)}
        N=[[[int(k==lookup[tuple(g[h[t]]for t in range(3))])for k in range(6)]for h in elements]for g in elements]
        self.assertNotEqual(N[1][2],N[2][1])
        system=generate_localization(N,k=0,categorical_dimensions=[1]*6)
        self.assertEqual(system['localizations'][0]['support'],[0])
        self.assertEqual(system['counts']['equations'],0)
        self.assertEqual(system['variables'],[])
        self.assertIn('dummy',to_singular(system))

if __name__=='__main__':unittest.main()
