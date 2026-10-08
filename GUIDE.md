# Reader and referee guide

The supplement has two purposes: applying the paper's criteria to **new fusion rings**, and verifying its worked applications. All entry points use Python 3.10 or later and the standard library. Run commands from this directory; the scripts also resolve their own frozen inputs when called by absolute path.

## 1. Supply a fusion ring

The basic input is a JSON tensor of nonnegative integers:

```
N[i][j][k] = coefficient of b_k in b_i b_j.
```

Its shape is `rank × rank × rank`, and `b_0` is the unit. Relabel all three indices consistently if the original unit has another label. The programs derive the dual involution from the unit coefficients. They check the unit, involution, Frobenius reciprocity, compatibility of duality with products, and associativity. They **do not assume commutativity**.

Either supply the tensor directly, as in `data/f660.json`, or wrap it in an object. Here is a complete input for the Fibonacci fusion ring, with `X*X = 1 + X`:

```json
{
  "name": "Fibonacci",
  "labels": ["1", "X"],
  "unit_index": 0,
  "fusion_matrices": [
    [[1, 0], [0, 1]],
    [[0, 1], [1, 1]]
  ]
}
```

Optional `dimensions` supplies a **positive rational dimension character** for an additional consistency check. Use integers or exact rational strings such as `"3/2"`, never floating-point approximations. Omit this field when the Frobenius–Perron dimensions are irrational, as for Fibonacci. It is separate from the categorical dimensions used in localization.

Validation establishes the stated ring axioms, not the existence of a categorification.

## 2. Search the zero- and one-spectrum criteria

```sh
python3 criteria.py spectrum my-ring.json --criterion both --mode first
python3 criteria.py spectrum data/rank6_zero_only.json --criterion both --mode count --output spectrum-report.json
python3 criteria.py spectrum data/rank6_one_only.json --criterion one --mode all --output all-witnesses.json
```

- `--criterion zero`, `one`, or `both` chooses the theorem to apply; `both` is the default.
- `--mode first` stops after finding a witness for each requested criterion. If a requested criterion has no witness, the search must finish to establish that fact.
- `--mode count` exhaustively counts witnesses and retains the first detailed witness for each criterion; this is the default.
- `--mode all` additionally retains every witness and may produce a large file.

Search counts refer to **ordered tuples**, without quotienting by symmetry. Zero-spectrum indices are `(i1,...,i9)`; one-spectrum indices are `(i0,i1,...,i9)`. Every returned detailed witness lists the six required coefficients, all spectrum terms, the opposite coefficient, the alternative dimension sums, and a separate Boolean check of every theorem hypothesis. The output states whether the search was exhaustive; partial counts are not total counts.

A verified witness excludes categorification over every algebraically closed field, by the corresponding theorem in the manuscript. These criteria require neither a pivotal structure nor categorical dimensions. Finding no witness is inconclusive. The exhaustive search has many index choices and can be expensive for large or dense rings; `first` is often the appropriate initial search.

Python API:

```python
from fusion_ring import load_fusion_ring, validate_fusion_ring
from spectrum import scan, evaluate_witness

N, metadata = load_fusion_ring("my-ring.json")
validate_fusion_ring(N, metadata.get("dimensions"))
report = scan(N, criterion="both", mode="first")
# To verify a known zero-spectrum tuple directly:
# detail = evaluate_witness(N, [i1,i2,i3,i4,i5,i6,i7,i8,i9], "zero")
```

## 3. Generate localization equations

This implementation generates equations over the rational numbers for **characteristic-zero spherical categorification**, under the additional hypotheses of the localization theorem. Their solutions may have algebraic, rather than rational, coordinates. The generator produces **necessary polynomial equations**; it does not decide consistency automatically. The geometry determines small, overlapping subsystems, which can be studied separately and then coupled.

The input object can contain a `localization` specification. Alternatively, place that specification in a separate JSON file and use `--config`:

```sh
python3 criteria.py localize examples/localization_f210.json --output f210-system.json
python3 criteria.py localize examples/localization_s3.json --output s3-system.json
python3 criteria.py localize examples/localization_fibonacci.json --output fibonacci-system.json
python3 criteria.py localize my-ring.json --config my-localization.json --output system.json --singular system.sing
```

`--singular` exports the same exact equations as a Singular input file. Exporting the file uses Python; solving it with Singular is a separate optional step.

For each center, provide its basis index `k`. The center and every element of its square support must be self-dual, and the square must be multiplicity-free. The program finds a basis object `b` for which the multiplicity of `b_k` in `b*b*` is odd, or checks the optional supplied `odd_witness`. This search also detects whether any reducible object can be an odd witness: for self-dual `k`, the cross terms from distinct simple summands occur in pairs, so an odd total requires an odd diagonal term. The program checks all these fusion-rule conditions. They are hypotheses of this particular localization criterion; the general spectrum criteria above do not require them.

Set `mode` to `"full"` for the full localization theorem, retaining three-index `x` and two-index `y` coordinates on the entire square support. Set it to `"corollary"` (the default) for the smaller system with two-index `x(a,b)=x(a,b,b)` coordinates. In corollary mode, each center's optional `subset` selects any subset of its square support; it need not contain the unit or the center. If omitted, it defaults to the whole support. Summation indices always run over the whole support. The first and second positions of corollary `x` are not identified merely by interchanging them. The generator uses separate variable names for different centers.

A minimal rational specification for the `Rep(S3)` input is:

```json
{
  "mode": "full",
  "categorical_dimensions": [1, 1, 2],
  "centers": [{"k": 2}]
}
```

For the symbolic Fibonacci input, use `{"mode": "full", "centers": [{"k": 1}]}`. A selected corollary subsystem of `Rep(S3)` is specified by `{"mode": "corollary", "categorical_dimensions": [1, 1, 2], "centers": [{"k": 2, "subset": [1]}]}`. The optional `ground_characteristic` field must be `0`. Unknown specification fields are rejected to catch misspelled input keys.

Categorical dimensions are essential extra data:

- With explicit `categorical_dimensions`, give a nonzero rational value for every basis element. The program verifies the unit, dual symmetry, and all dimension-character identities. These values must be justified as the categorical dimensions of the spherical categorification under investigation. A positive dimension vector may be substituted when a pseudo-unitarity argument warrants it, as for F210.
- With `categorical_dimensions` omitted or `null`, the dimensions remain symbolic. The generator includes the dimension-character equations, dual symmetry, unit normalization, and equations imposing nonzero dimensions by introducing inverse variables. This permits irrational dimension assignments, including Fibonacci, without approximating them.

Inconsistency with one chosen numerical dimension character excludes only the corresponding class of spherical categorifications unless a mathematical argument shows that all relevant categorifications admit that choice. Symbolic dimensions avoid imposing a single rational character.

Couplings use an ordered list such as `"couplings": [[1, 3]]`, together with `"compatible_bases": true`. For a pair `[k,l]`, the centers must be distinct and `l` must belong to both selected supports. The program checks these requirements and generates the stated coupling equation, with its full intersection sum. The declaration records the simultaneous compatible formal basis choice guaranteed under a hypothetical categorification; it does not require a pre-existing categorification or numerical basis data. Numerical solutions computed in independently fixed bases cannot simply be combined without the compatible change of bases. No additional equality between coordinates of different centers is imposed. The example F210 specification exhibits the two centers and coupling used in the paper.

The complete, executable specifications in `examples/` are the recommended starting points. Their generated JSON includes variable names, exact rational coefficients, equations, and the checked hypothesis data. Preserve this information with any certificate obtained from the equations. A claimed inconsistency should be accompanied by an exact certificate or another independently checkable mathematical argument. A solution of a localized system alone is not a construction of a fusion category.

Python API:

```python
from fusion_ring import load_fusion_ring
from localization import generate_from_spec, to_singular

N, metadata = load_fusion_ring("examples/localization_fibonacci.json")
system = generate_from_spec(N, metadata["localization"])
singular_input = to_singular(system)
```

## 4. Verify the paper's applications

```sh
python3 verify.py
```

Run without `-O`, which disables assertion-based certificate checks. A clean run with site packages disabled (`python3 -S verify.py`) also checks that verification uses the standard library.

| Mathematical claim | Inputs and certificate | Verification |
|---|---|---|
| F210 fusion rules and dimensions | `data/f210.json` | Fusion axioms and dimension character |
| Localized F210 systems | `data/f210-equations.json` | `verify_localization.py`: all 27 corollary instances give exactly the 12 nonzero equations up to nonzero rational scaling, for each localization |
| F210 in characteristic zero | `certificates/f210-basis.txt`, `certificates/f210-lift.txt.gz` | `verify_f210.py`: three linear eliminations; 28 integer ideal-membership identities; 14 spanning monomials; determinant 85 modulo 101 for the 196-by-196 coupling matrix |
| Exactly 14 distinct solutions of each separate F210 localization | `certificates/f210-solutions.json` | `verify_f210_solutions.py`: squarefree degree-14 parameter polynomial, exact coordinate substitution, and the certified upper bound |
| F210 character table and dual-Burnside property | `data/f210.json` | `verify_f210_characters.py`: exact cyclotomic arithmetic checks all character products, codegrees, and zero entries |
| F210 in exceptional positive characteristics | `data/f210.json` | Exact identities `B=175R`, `R^2=210R`, `B^2=36750B` |
| F660 over any algebraically closed field | `data/f660.json` | Every numerical hypothesis of the displayed zero-spectrum witness |
| Independence of the two spectrum criteria | `data/rank6_zero_only.json`, `data/rank6_one_only.json` | Exhaustive admissible-tuple searches, without symmetry identifications |
| Census-based unitary classification | `census/` | Exact primary-3 witnesses and group-ring basis maps |
| Frobenius–Schur counterexamples and bounded group search | `optional/check_small_groups.g`, `certificates/small-groups-output.txt` | Optional GAP search through order 128 |
| Smallest finite simple-group counterexample | `optional/check_simple_groups.g`, `certificates/simple-groups-output.txt` | Optional GAP character-table checks for all 34 nonabelian simple groups of order at most 126000; the only counterexample is PSU(3,5) |

Expected numerical results:

- F210 character codegrees, in displayed column order: `210,6,7,7,7,5,5`.
- F660 displayed witness: `(i1,...,i9)=(1,3,4,1,1,3,4,2,2)`. The six required coefficients are all 1; the eight spectrum terms vanish; the two alternative-sum triples are `(5,1,8)` and `(1,8,5)`.
- `rank6_zero_only`: 12 zero-spectrum witnesses and no one-spectrum witness.
- `rank6_one_only`: no zero-spectrum witness and 96 one-spectrum witnesses.
- F210 coupling certificate: 28 identities, 14 spanning monomials per localization, matrix size 196, determinant residue 85 modulo 101.
- Verified census subset: 33 inputs, comprising 27 primary-3 exclusions, F210, F660, and four group character rings.
- Bounded group search, with both `X` and `Y` simple: `SmallGroup(72,41)` and `SmallGroup(128,i)` for `i=764,801,802`; the exceptional multiplicities displayed by the search are 2. The order-72 minimality statement has this simple-`X` restriction; it is not a minimality statement allowing arbitrary reducible `X`.

## 5. Scope of the certificates and optional regeneration

The F210 checker first verifies polynomial identities over the integers. The modular determinant calculation then proves that a specified **rational matrix** is nonsingular, with all denominators controlled at the prime 101. This is not an inference of characteristic-zero inconsistency from a modular Gröbner basis. See [CERTIFICATE.md](CERTIFICATE.md) for the complete argument.

For exceptional positive characteristics, the identity `B^2=36750B` forces the product of basis-character values to vanish in characteristics 2, 3, 5, and 7. The lifting argument in the other characteristics is a mathematical input proved in the manuscript. The checker verifies the specified arithmetic and combinatorial statements; the manuscript proves their implications for categorification.

The census verification concerns the supplied 33 rings under the bounds stated in `census/README.md`. Completeness of the classification and the primary-3 theorem are cited inputs. The broader historical 505-ring census described in the introduction is not a claim that this program independently repeats every enumeration in that census.

`optional/regenerate_f210.py` invokes Singular's `liftstd` and writes to `generated-certificates/`. `optional/check_generated.py` checks the resulting identities using exact Python arithmetic. Regeneration may produce a different lifting matrix; the identities are what matters. GAP is used by the optional group-database regeneration scripts. Their recorded versions are listed in the repository README.

The optional historical check `gap -q optional/check_simple_groups.g` uses GAP's CTblLib character-table package (recorded version 1.3.7). In addition to the bounded simple-group search, it checks the dual-Burnside property for the ten smallest sporadic groups, with exceptions M22 and M24. The general finite-simple-group statement in the manuscript is supported by the cited mathematical source; this finite check does not replace that theorem.
