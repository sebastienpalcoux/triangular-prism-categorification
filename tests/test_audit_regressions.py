"""Regression tests for the defects reproduced in the independent audit."""
from fractions import Fraction
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fusion_ring import FusionRingError, validate_fusion_ring
from localization import _Equations
from polynomials import parse
from spectrum import evaluate_witness


def cyclic_three():
    return [[[int(c == (a + b) % 3) for c in range(3)]
             for b in range(3)] for a in range(3)]


class AuditRegressionTests(unittest.TestCase):
    def test_wrong_duality_cannot_certify_a_categorifiable_ring(self):
        N = cyclic_three()
        indices = [0, 1, 1, 0, 1, 0, 0, 1, 0]
        self.assertFalse(evaluate_witness(N, indices, 'zero')['passes'])
        self.assertFalse(evaluate_witness(N, indices, 'zero', dual=[0, 2, 1])['passes'])
        for dual in ([0, 1, 2], [0, 2], [0, True, 2]):
            with self.subTest(dual=dual), self.assertRaisesRegex(FusionRingError, 'dual'):
                evaluate_witness(N, indices, 'zero', dual=dual)

    def test_mixed_row_containers_preserve_commutativity(self):
        N = cyclic_three()
        mixed = copy.deepcopy(N)
        mixed[0][1] = tuple(mixed[0][1])
        self.assertEqual(validate_fusion_ring(N), validate_fusion_ring(mixed))
        self.assertTrue(validate_fusion_ring(mixed)['commutative'])

    def test_certificate_parser_rejects_incomplete_or_ambiguous_syntax(self):
        for expression in ('', '-', '1+', 'U1--U2', 'U1+-U2', 'U1*', '*U1', 'U1^'):
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                parse(expression, ['U1', 'U2'])
        self.assertEqual(parse('U1-U2', ['U1', 'U2']), {(1, 0): 1, (0, 1): -1})
        self.assertEqual(parse('0', ['U1']), {})
        self.assertEqual(parse('1/2*U1', ['U1'], rational=True), {(1,): Fraction(1, 2)})
        with self.assertRaises(ValueError):
            parse('1/0*U1', ['U1'], rational=True)

    def test_negative_singleton_equations_have_the_same_primitive_row(self):
        equations = _Equations(['x'])
        for coefficient in (1, -1, Fraction(-3, 2)):
            equations.add({(1,): coefficient}, {'coefficient': str(coefficient)})
        self.assertEqual(len(equations.rows), 1)
        self.assertEqual(equations.rows[0]['text'], 'x')
        self.assertEqual(len(equations.rows[0]['origins']), 3)

    def test_cli_rejects_all_input_and_output_aliases_before_writing(self):
        data = json.loads((ROOT / 'examples/localization_fibonacci.json').read_text())
        for alias_kind in ('same', 'symlink', 'hardlink'):
            for target in ('input', 'config', 'output'):
                for destination in ('json', 'singular'):
                    with self.subTest(alias=alias_kind, target=target, destination=destination):
                        with tempfile.TemporaryDirectory() as folder:
                            root = Path(folder)
                            source = root / 'ring.json'
                            source.write_text(json.dumps(data))
                            config = root / 'config.json'
                            config.write_text(json.dumps(data['localization']))
                            output = root / 'output.json'
                            output.write_text('existing output\n')
                            protected = {'input': source, 'config': config, 'output': output}[target]
                            alias = protected if alias_kind == 'same' else root / 'alias'
                            if alias_kind == 'symlink':
                                alias.symlink_to(protected)
                            elif alias_kind == 'hardlink':
                                os.link(protected, alias)
                            arguments = [sys.executable, '-S', str(ROOT / 'criteria.py'),
                                         'localize', str(source), '--config', str(config)]
                            if target == 'output':
                                arguments += ['--output', str(output), '--singular', str(alias)]
                            elif destination == 'json':
                                arguments += ['--output', str(alias)]
                            else:
                                arguments += ['--output', str(root / 'new.json'),
                                              '--singular', str(alias)]
                            paths = [source, config, output]
                            before = [path.read_bytes() for path in paths]
                            done = subprocess.run(arguments, text=True, capture_output=True)
                            self.assertEqual(done.returncode, 2, done.stdout + done.stderr)
                            self.assertIn('distinct paths', done.stderr)
                            self.assertEqual([path.read_bytes() for path in paths], before)
                            self.assertFalse((root / 'new.json').exists())


if __name__ == '__main__':
    unittest.main()
