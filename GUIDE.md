# Reader and referee guide

Run `python3 verify.py` from this directory. It is also safe to invoke the script by its absolute path from another directory. The Python standard library is sufficient; `python3 -S verify.py` also works. Do not use `python3 -O`, which disables assertion-based checks.

## Manuscript-to-file map

| Mathematical claim | Inputs and certificate | Verification |
|---|---|---|
| F210 fusion rules and dimensions | `data/f210.json` | Fusion axioms and dimension character in `verify.py` |
| Localized F210 polynomial systems | `data/f210-equations.json` | `verify_localization.py`: all 27 corollary instances yield precisely the 12 listed nonzero equations, up to nonzero rational scaling; both localizations are checked |
| F210 in characteristic zero | `certificates/f210-basis.txt`, `certificates/f210-lift.txt.gz` | `verify_f210.py`: three linear eliminations; 28 integer identities; 14 spanning monomials; a 196-by-196 matrix with determinant 85 modulo 101 |
| Formal codegrees of F210 | `data/f210.json` | Exact characteristic polynomial of the Casimir multiplication matrix |
| F210 in exceptional positive characteristics | `data/f210.json` | Exact products `B=175R`, `R^2=210R`, `B^2=36750B` |
| F660 over any algebraically closed field | `data/f660.json` | Every numerical hypothesis of the displayed zero-spectrum witness |
| Independence of zero/one-spectrum criteria | `data/rank6_zero_only.json`, `data/rank6_one_only.json` | `spectrum.py`: exhaustive admissible-tuple enumeration, without identifying tuples by symmetry |
| Census-based unitary classification | `census/` | `census/verify_census.py`: exact primary-3 witnesses and group-ring basis maps |
| Frobenius–Schur counterexample and bounded group search | `optional/check_small_groups.g`, `certificates/small-groups-output.txt` | Optional independent GAP search; the recorded output was reproduced with GAP 4.12.1 |

## Expected numerical output

- F210 formal codegrees: `5,5,6,7,7,7,210`.
- F660 witness indices `(i1,...,i9)=(1,3,4,1,1,3,4,2,2)`, using zero-based labels. All six required fusion coefficients equal 1. The eight spectrum terms vanish. The two alternative-sum triples are `(5,1,8)` and `(1,8,5)`.
- `rank6_zero_only`: 12 zero-spectrum witnesses and no one-spectrum witness.
- `rank6_one_only`: no zero-spectrum witness and 96 one-spectrum witnesses.
- F210 certificate: 28 ideal-membership identities, 14 spanning monomials, matrix size 196, determinant residue 85 in characteristic 101.
- Census: 33 inputs, 27 negative primary-3 certificates, F210 and F660, and four group character rings.
- Optional bounded GAP search: the exceptional groups are `SmallGroup(72,41)` and `SmallGroup(128,i)` for `i=764,801,802`. Every displayed exceptional multiplicity is 2.

## What is and is not being certified

The F210 checker verifies polynomial identities over the integers before using a modular calculation solely to prove that a **rational matrix** is nonsingular. It does not infer characteristic-zero inconsistency from a modular Gröbner basis. The precise argument is in `CERTIFICATE.md`.

For positive characteristic, no algebraic-closure enumeration is needed. The identity `B^2=36750B` forces the product of all basis-character values to vanish in characteristics 2, 3, 5, and 7. The lifting argument for the other characteristics is mathematical and appears in the manuscript.

The computational files verify the specified algebraic and combinatorial statements. The deduction from the graphical identities to categorification obstructions is proved in the manuscript. Enumeration completeness and the primary-3 theorem are cited mathematical inputs, not conclusions of the checker.

## Optional regeneration

`optional/regenerate_f210.py` invokes Singular's existing `liftstd` routine and saves its output in a separate directory. `optional/check_generated.py` verifies those output files with the same integer arithmetic and modular nonsingularity checks. A different Singular version can produce a different lifting matrix; the identities, rather than a software-specific textual output, are the certificate.

GAP is retained for the two group-database computations. The bounded search inspects every group in the SmallGroups library through order 128 and every relevant pair of irreducible characters. The core Python checker does not reproduce a finite-group database.
