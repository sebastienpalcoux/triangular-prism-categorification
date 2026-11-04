"""Exact scalar checks in complex C2-graded vector spaces.

We use the usual left duality and the pivotal character chi(g)=-1 in
skeletal coordinates. The pentagon test constructs its changes of basis
from fusion maps and then applies the coordinate conversion of §5.2.
The remaining tests evaluate the two sides of particular duality pairings.
"""

from fractions import Fraction as Q
from itertools import product
import unittest


def chi(degree):
    """Categorical dimension of the homogeneous line of this C2 degree."""
    return Q((-1) ** (degree % 2))


class PointedPivotalConventionTests(unittest.TestCase):
    def test_dimension_in_converted_pentagon(self):
        # With the ordinary associator, a change of fusion bases produces
        # the ratio of the two products of vertex coefficients.
        basis = {(0, 0): Q(2), (0, 1): Q(3), (1, 0): Q(5), (1, 1): Q(7)}

        def associator(a, b, c):
            return basis[a, b] * basis[a ^ b, c] / (basis[b, c] * basis[a, b ^ c])

        def tetrahedron(a, b, c):
            # The right intermediate label of F(a,b,c) is b+c in C2.
            return associator(a, b, c) / chi(b ^ c)

        unequal_without_dimension = 0
        for a, b, c, d in product((0, 1), repeat=4):
            with self.subTest(input_degrees=(a, b, c, d)):
                left = [(a ^ b, c, d), (a, b, c ^ d)]
                right = [(a, b, c), (a, b ^ c, d), (b, c, d)]
                self.assertEqual(
                    associator(*left[0]) * associator(*left[1]),
                    associator(*right[0]) * associator(*right[1]) * associator(*right[2]),
                )
                two = tetrahedron(*left[0]) * tetrahedron(*left[1])
                three = tetrahedron(*right[0]) * tetrahedron(*right[1]) * tetrahedron(*right[2])
                self.assertEqual(two, chi(b ^ c) * three)
                unequal_without_dimension += two != three

        self.assertEqual(unequal_without_dimension, 8)

    def test_oriented_cap_retains_pivotal_morphism(self):
        # ev_{X*,-} = ev_{X**,+} (a_{X*} tensor id_{X**}).
        # For an odd line the two orientations have different scalars.
        for degree in (0, 1):
            with self.subTest(degree=degree):
                ordinary_evaluation = Q(1)
                opposite_evaluation = chi(degree)
                pivotal_on_dual = chi(degree)
                h = Q(7, 5)
                self.assertEqual(opposite_evaluation * h, ordinary_evaluation * pivotal_on_dual * h)
                self.assertEqual(opposite_evaluation * h == ordinary_evaluation * h, degree == 0)

    def test_trace_dual_and_left_evaluation_in_skeletal_coordinates(self):
        for a, b in product((0, 1), repeat=2):
            c = a ^ b
            for mu in (Q(2), Q(3, 5), Q(-7)):
                with self.subTest(a=a, b=b, mu=mu):
                    trace_dual_q = 1 / (chi(c) * mu)
                    self.assertEqual(chi(c) * mu * trace_dual_q, 1)
                    # Closing with the ordinary evaluation alone has
                    # coefficient 1, whereas the pivotal trace includes chi(c).
                    evaluation_dual = trace_dual_q * chi(c)
                    self.assertEqual(mu * evaluation_dual, 1)
                    self.assertEqual(mu * trace_dual_q == 1, c == 0)


if __name__ == "__main__":
    unittest.main()
