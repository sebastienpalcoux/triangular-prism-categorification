# Pivotal factors and contraction-dual normalization

This note explains the normalization used in the manuscript's graphical prism equations and the PE–TPE correspondence. The categorical proofs are in Sections 4–6; the finite scalar regression tests described below complement those proofs.

## A single pivotal coupon

In a realization with `X** = X` on objects and `a_X² = id_X`, a simple object has `a_X = epsilon_X id_X`, where `epsilon_X² = 1`. These assumptions do not imply that the pivotal structure is the identity, or that double duality is the identity on morphisms.

The general prism equation contains both the categorical dimension `d_X` and a single pivotal coupon on an internal `X*` strand. In this realization that coupon contributes `epsilon_X`. After it is extracted from the tetrahedral contraction, the weight is therefore

```text
D_X = epsilon_X d_X.
```

The assumption `a_X² = id_X` removes pairs of pivotal coupons, not this single coupon. The graphical spherical and simplified prism formulas retain this weight. The general equation, before these object identifications, retains the actual pivotal morphism.

For example, take complex `C2`-graded vector spaces with ordinary left duality and pivotal character `chi(g) = -1`. Then `d_g = -1` and `epsilon_g = -1`, while every ordinary evaluation and coevaluation of the chosen homogeneous lines has coefficient `1`. With prism edge labels

```text
X1 = X2 = X3 = X7 = X8 = X9 = 1,
X4 = X5 = X6 = g,
```

and unit vertex coefficients, the two-tetrahedron side is `1`. The unique internal object on the three-tetrahedron side is `g`; the weight `epsilon_g d_g = 1` gives the same value. The weight `d_g` alone would give `-1`.

## Trace-dual maps and invariant dual tensors

Let `mu: A tensor B -> C` have trace-dual `q: C -> A tensor B`, so `tr(nu q) = delta(nu,mu)`. Set `Z = C tensor B* tensor A*`. In the same object realization, the **right** invariant-dual tensor is

```text
'tilde(mu) = (q a_C tensor id_C*) coev_C.
```

The order of the pairing matters:

```text
b_(Z*) (tilde(nu), 'tilde(mu))
  = tr(a_C^(-1) nu q a_C)
  = tr(nu q).
```

The last equality uses cyclicity. The left dual for `b_Z`, with the arguments in the opposite order, instead uses `(a_(A tensor B) q tensor id_C*) coev_C`. No equality between these two expressions is presumed when double duality acts nontrivially on morphisms.

For a simple output, write `qbar = q a_C`, and let `tau_C` denote closure using left duality without a pivotal coupon. Then `tau_C = epsilon_C tr_C`. Coefficient extraction in the manuscript consequently gives `F = D_i6 T`, and substitution in the five pentagon factors leaves the internal weight `D_i0`. This agrees with the graphical prism formula.

In the sign example with `A = 1`, `B = C = g` and `mu = 1`, the trace-dual map has coefficient `q = -1`; multiplication by `a_C` makes the right invariant-dual coefficient `+1`. The worked example in Section 5 tracks this normalization and the subsequent corner coupons separately.

## Why the localization generator still uses categorical dimensions

Localization is performed in an equivalent **strictly pivotal realization**, where every `epsilon_i = 1` and hence `D_i = d_i`. This preserves the categorical dimensions; it does not assert that they are positive. Selfduality is handled separately: the localization hypotheses imply indicator one for the relevant simple objects, permitting symmetric isomorphisms `kappa_i: X_i -> X_i*`. Cups, caps, incident vertex morphisms and their contraction-dual bases are transported through the same choices. Coupled localizations use one common collection of these choices.

Thus the general localization equations and their implementation continue to use the specified or symbolic **categorical dimensions**. Strictification alone does not identify a simple object literally with its chosen dual, and selfduality up to isomorphism alone does not identify `epsilon_i` with its Frobenius–Schur indicator.

## Exact finite regression checks

`tests/test_pivotal_conventions.py` uses only Python's standard library. It checks:

- all 16 admissible `C2` prism labellings, with unequal exact rational vertex and internal-basis coefficients: the retained pivotal factor gives equality in every case, while omitting it fails in exactly 8 cases;
- right invariant-dual normalization in all four `C2` fusion channels with three different nonzero rational basis coefficients.

Run these checks with the complete verification command `python3 verify.py`, or alone with:

```sh
python3 -B -S -m unittest discover -s tests -p test_pivotal_conventions.py -v
```

These are concrete scalar controls for the stated pointed category, not a general implementation or proof of graphical categorical calculus.
