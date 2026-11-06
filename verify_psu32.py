#!/usr/bin/env python3
"""Verify the order-72 positivity counterexample using ordinary complex characters.

Realize PSU(3,2) = F_3^2 ⋊ Q8 on nine affine points. The two characters
are those in Yilong Wang, “On the positivity conjecture”, Sections 2–3.
"""
from itertools import product


def verify():
    """Check irreducibility, both indicators and the multiplicity of Y in X* ⊗ X."""
    points = tuple(product(range(3), repeat=2))

    def permutation(action):
        images = (action(x, y) for x, y in points)
        return tuple(3 * (x % 3) + y % 3 for x, y in images)

    def compose(left, right):
        return tuple(left[i] for i in right)

    identity = tuple(range(9))
    minus = permutation(lambda x, y: (-x, -y))
    a = permutation(lambda x, y: (-x - y, -x + y))
    b = permutation(lambda x, y: (-y, x))
    ab = compose(a, b)
    if compose(a, a) != minus or compose(b, b) != minus or compose(b, a) != compose(minus, ab):
        raise ValueError("The generators do not satisfy the quaternion relations")

    positive = (identity, a, b, ab)
    quaternion = positive + tuple(compose(minus, q) for q in positive)
    translations = [permutation(lambda x, y: (x + u, y + v)) for u, v in points]
    linear_part = {compose(t, q): q for t in translations for q in quaternion}
    order = len(linear_part)
    if order != 72:
        raise ValueError(f"Expected 72 distinct affine permutations, found {order}")

    # X is the permutation representation with its constant line removed.
    # Y is inflated from the degree-two character of Q8, with values 2,-2,0.
    chi_x = {g: sum(i == j for i, j in enumerate(g)) - 1 for g in linear_part}
    quaternion_trace = {identity: 2, minus: -2}
    chi_y = {g: quaternion_trace.get(q, 0) for g, q in linear_part.items()}
    if (chi_x[identity], chi_y[identity]) != (8, 2):
        raise ValueError("Unexpected dimensions for X and Y")
    if any(sum(value**2 for value in chi.values()) != order for chi in (chi_x, chi_y)):
        raise ValueError("The characters of X and Y must each have norm one")

    # Character sums are integral here. Compare before dividing by |G|.
    squares = [compose(g, g) for g in linear_part]
    indicator_x = sum(chi_x[square] for square in squares)
    indicator_y = sum(chi_y[square] for square in squares)
    multiplicity = sum(chi_x[g]**2 * chi_y[g] for g in linear_part)
    sums = (indicator_x, indicator_y, multiplicity)
    if sums != (order, -order, 2 * order):
        raise ValueError(f"Unexpected indicator and multiplicity numerators: {sums}")
    result = tuple(total // order for total in sums)
    print(f"PSU(3,2): {order} elements, irreducible X and Y of degrees 8 and 2")
    print(f"nu_2(X) = {result[0]}, nu_2(Y) = {result[1]}, multiplicity of Y in X* x X = {result[2]}")
    return result


if __name__ == "__main__":
    verify()
