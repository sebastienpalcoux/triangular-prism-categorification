# Pivotal conventions in the triangular prism equations

## The dimension on the interior edge

In the general pivotal equation of §4.2, the factor on a simple interior edge labelled by `X` is its categorical dimension `d_X`. The pivotal maps appearing within the tetrahedral coefficients are actual morphisms. In the general spherical equation of §4.3, reversing the middle edge trades an evaluation of one orientation for the opposite evaluation composed with the pivotal map on `X*`. The external corner changes, including the displayed applications of `σ⁻²`, are part of the same equation. Extracting a scalar sign from one of these morphisms without changing the evaluation convention would change the coefficient.

The simplified triangular prism equation uses a strictly pivotal realization, with `a_X = id_X` and double duality the identity on objects **and** morphisms. Pivotal equivalence preserves traces, dimensions and indicators. The dimension on the interior edge remains `d_X`. A self-dual simple object still requires a chosen isomorphism to its dual when the calculation identifies their labels.

## The pairing used in §5.2

For a basis vector `μ: A ⊗ B → C`, its trace-dual is the map `q_μ: C → A ⊗ B` characterized by `tr(ν ∘ q_μ) = δ(ν,μ)` on the chosen basis. In the strictly pivotal realization, the bent tensor

    (q_μ ⊗ id_{C*}) ∘ coev_C : 1 → A ⊗ B ⊗ C*

is dual to the bent basis for the evaluation pairing. If `Z = C ⊗ B* ⊗ A*`, the bent basis vectors are maps from the unit to `Z`, and their partners are maps from the unit to `Z*`. The categorical dual of a bent basis vector instead has type `Z* → 1`.

Schur's lemma in the coefficient extraction gives `F = d_i T`, where `i` is the simple object at the right intermediate edge of that associator coefficient. In the two routes of the pentagon the common factors are `d_{i_1} d_{i_6}`. Their cancellation leaves `d_{i_0}` inside the sum over the interior simple object. The localization generator uses these categorical dimensions, specified exactly or left symbolic. Neither strictification nor trace duality imposes positivity.

## What the pointed tests check

Run the three scalar tests from this directory:

```sh
python3 -B -S -m unittest discover -s tests -p test_pivotal_conventions.py -v
```

They use complex `C2`-graded vector spaces, with ordinary left duality and the pivotal character `χ(g) = −1`. Thus the odd line has categorical dimension `−1`. This skeletal coordinate choice has nonidentity pivotal map on the odd line. The formulas in §5.2 for an already strictly pivotal realization cannot be applied to these bare ordinary cups without transporting the other data.

- For each of the 16 assignments of degrees to the four factors of a pentagon, the test obtains `F` as the ratio of two products of chosen fusion-basis coefficients in the trivial associator. Converting each coefficient to `T = F/d_i` gives equality with the single interior factor `d_{i_0}`. Omitting that factor fails for exactly the eight assignments with odd interior label. This checks the dimension cancellation in coordinates; it does not independently derive the graphical change of basis.
- The cap test uses the two orientations of evaluation. Their coefficients on a homogeneous line of degree `j` are `1` and `χ(j)` in these coordinates. It checks the local edge reversal of §4.3 after composing with the pivotal map, including the case `χ(j) = −1`.
- For each pair of input degrees and three nonzero vertex coefficients, the last test compares the ordinary evaluation pairing with the pivotal trace. In these skeletal scalar coordinates the trace-dual has coefficient `1/(χ(c) μ)`, and its partner for bare evaluation has coefficient `1/μ`, where `c` is the output degree. The tests concern one pointed category and particular coordinates. The general statements are proved in the manuscript.
