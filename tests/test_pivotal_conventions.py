"""Finite pointed regressions for the manuscript's pivotal conventions.

These exact scalar checks use Vec(C2) with its ordinary left duality and
the spherical pivotal character chi(g) = -1. They guard the right-dual
normalization and the internal pivotal factor found by this example.
They are not a general categorical proof or a general PE--TPE converter.
"""

from fractions import Fraction as Q
from itertools import product
import unittest


def chi(degree):
    """The nontrivial pivotal character of the additive group C2."""
    return Q((-1) ** (degree % 2))


def admissible_prism_labels():
    """Enumerate edge degrees for which all six external vertices exist."""
    for labels in product((0, 1), repeat=9):
        x1, x2, x3, x4, x5, x6, x7, x8, x9 = labels
        # Duality fixes each C2 degree. These are the six ordered targets
        # of alpha_1, ..., alpha_6 in the simplified prism equation.
        vertex_degrees = (
            x4 + x1 + x6,
            x5 + x2 + x4,
            x6 + x3 + x5,
            x9 + x1 + x7,
            x7 + x2 + x8,
            x8 + x3 + x9,
        )
        if all(degree % 2 == 0 for degree in vertex_degrees):
            yield labels


class PointedPivotalConventionTests(unittest.TestCase):
    def test_all_c2_prisms_retain_the_internal_pivotal_factor(self):
        # Evaluations, coevaluations and rotations of homogeneous lines
        # have coefficient 1. Thus each tetrahedron is the product of
        # its four vertex coefficients. Unequal coefficients also check
        # that each internal basis vector pairs with its reciprocal.
        a1, a2, a3, a4, a5, a6 = map(Q, (2, 3, 5, 7, 11, 13))
        b0, b1, b2, b3 = map(Q, (17, 19, 23, 29))
        cases = omitted_coupon_failures = 0

        for labels in admissible_prism_labels():
            with self.subTest(edge_degrees=labels):
                cases += 1
                x1, x2, x3, x4, x5, x6, x7, x8, x9 = labels
                simple_candidates = {
                    (x4 + x7) % 2,
                    (x5 + x8) % 2,
                    (x6 + x9) % 2,
                }
                self.assertEqual(len(simple_candidates), 1)
                internal_degree = simple_candidates.pop()

                # Z0 and Z1,Z2,Z3 must have degree zero for their basis
                # vectors to exist; all four are homogeneous unit lines.
                self.assertEqual((x3 + x2 + x1) % 2, 0)
                for target_degree in (
                    x4 + internal_degree + x7,
                    x5 + internal_degree + x8,
                    x6 + internal_degree + x9,
                ):
                    self.assertEqual(target_degree % 2, 0)

                lhs = (a3 * a1 * a2 / b0) * (a4 * a6 * a5 * b0)
                rhs_tetrahedra = (
                    (a1 * a4 * b1 / b3)
                    * (a2 * a5 * b2 / b1)
                    * (a3 * a6 * b3 / b2)
                )
                dimension = chi(internal_degree)
                pivotal_coupon = chi(internal_degree)
                corrected_rhs = dimension * pivotal_coupon * rhs_tetrahedra
                omitted_coupon_rhs = dimension * rhs_tetrahedra

                self.assertEqual(lhs, corrected_rhs)
                self.assertEqual(
                    omitted_coupon_rhs == lhs,
                    internal_degree == 0,
                )
                omitted_coupon_failures += omitted_coupon_rhs != lhs

        self.assertEqual(cases, 16)
        self.assertEqual(omitted_coupon_failures, 8)

    def test_right_dual_uses_q_a_c_and_not_the_bare_trace_dual(self):
        for a, b in product((0, 1), repeat=2):
            c = (a + b) % 2
            for mu in (Q(2), Q(3, 5), Q(-7)):
                with self.subTest(a=a, b=b, mu=mu):
                    dimension_c = chi(c)
                    trace_dual_q = 1 / (dimension_c * mu)
                    trace_pairing = dimension_c * mu * trace_dual_q
                    self.assertEqual(trace_pairing, 1)

                    # coev_C has coefficient 1. The right invariant dual
                    # is (q a_C tensor id_C*) coev_C, so its evaluation
                    # pairing with tilde(mu) has coefficient mu*q*chi(c).
                    right_invariant_dual = trace_dual_q * chi(c)
                    self.assertEqual(mu * right_invariant_dual, 1)
                    self.assertEqual(mu * trace_dual_q, chi(c))
                    self.assertEqual(mu * trace_dual_q == 1, c == 0)


if __name__ == "__main__":
    unittest.main()
