# Triangular prism equations and categorification

Computational material for **Triangular prism equations and categorification**, by Zhengwei Liu, Sebastien Palcoux, Yunxiang Ren, and Gert Vercleyen.

The programs apply the paper's **zero-spectrum, one-spectrum, and localization criteria to user-supplied fusion rings**. F210 and F660 are worked applications. Start with [GUIDE.md](GUIDE.md) for the input format, assumptions, examples, and connection to the manuscript. [CERTIFICATE.md](CERTIFICATE.md) explains the exact characteristic-zero certificate for F210.

For a new fusion ring, the workflow is: search the spectrum criteria using the fusion rules alone; when the localization hypotheses hold, generate small necessary systems with justified or symbolic categorical dimensions; then add the coupling equations for compatible localizations. The last step can produce an obstruction even when the individual local systems have solutions. The guide separates these steps and explains what each output establishes.

## Apply the criteria to another fusion ring

Use Python 3.10 or later with its standard library. Supply exact fusion coefficients in JSON, with `N[i][j][k]` the coefficient of basis element `k` in the product `i*j`, and unit index `0`.

```sh
python3 criteria.py spectrum my-ring.json --criterion both --mode first
python3 criteria.py localize examples/localization_fibonacci.json --output system.json
```

The spectrum search accepts noncommutative fusion rings and reports every numerical hypothesis of each returned witness. A witness excludes fusion categorification over any field, without assuming pivotality or pseudo-unitarity. The argument uses semisimplicity. Both criteria first detect an obstruction in rank four, on the fusion ring of Proposition 6.10. [The minimum-rank proof](SPECTRUM_RANK.md) gives the theoretical argument and the two witnesses. Run `python3 verify_spectrum_rank.py` for the small-rank verification in §6.3. Localization additionally needs a self-dual center, suitable multiplicity-free support, odd-multiplicity data, and categorical dimensions (specified exactly or left symbolic). The guide explains these choices and how to couple localizations. A spectrum witness is an obstruction; generating localization equations is a first step toward checking their consistency. Absence of an obstruction does not prove categorifiability.

## Verify the paper's applications

Run without Python's optimization flag, since the frozen certificate checks use assertions:

```sh
python3 verify.py
```

This command uses exact arithmetic in the Python standard library. Its largest input is a 24 MB compressed integer certificate, which allows direct checking of identities that were more expensive to discover. It checks the displayed fusion rules, both F210 localizations and their coupling, the F210 character table and positive-characteristic identities, the F660 witness, the rank-six spectrum examples, the order-72 character calculation in Remark 3.6, and the stated census exclusions. General-purpose implementations also have tests on categorifiable fusion rings, including a noncommutative group ring and Fibonacci fusion rules.

To reproduce Remark 3.6 alone, run `python3 -S verify_psu32.py`. It constructs the 72 affine permutations of PSU(3,2) and computes character norms, Frobenius–Schur indicators and the stated multiplicity. The larger searches over groups through order 128 and over finite simple groups have separate GAP scripts and archived outputs.

The general triangular prism equation retains its pivotal morphisms. In the strictly pivotal coordinates of §5.2, the PE–TPE change of variables uses the categorical dimension on the interior edge. The pointed tests check dimension cancellation, reversal of a cap and the evaluation pairing in explicit scalar coordinates. See [PIVOTAL_CONVENTIONS.md](PIVOTAL_CONVENTIONS.md).

The census-based conclusion concerns nonpointed, unitary, 1-Frobenius, simple integral fusion categories of rank at most eight and Frobenius–Perron dimension at most 20000, up to Grothendieck equivalence. Completeness of the underlying census is an input from the cited classification. See [census/README.md](census/README.md) and its provenance file for the exact checks, the treatment of the dimension bound and the unrestricted rank-four result.

## Near-group motivation

A concrete application of the geometric formulation is Huixuan He, Zhengwei Liu, Fan Lu, Sebastien Palcoux, and Yunxiang Ren, *Triangular prism equations for near-group categories: algebraic solutions and a unitary conjecture*, [arXiv:2610.11652](https://arxiv.org/abs/2610.11652) (2026).

Applying TPE to `G+|G|` for finite abelian `G`, [Theorem 2.3](https://arxiv.org/html/2610.11652v1) derives a variant of Izumi's equations **without assuming unitarity**: it characterizes the existence of spherical categorifications with trivial Frobenius–Schur indicators over fields in which `|G|` is invertible. Remark 2.2 explains the useful form: the essential identity equates a weighted quadratic sum to a cubic term and a boundary correction, whereas Izumi's equation sums cubic terms. In characteristic two it gives `β(2h) = β(h)^4`. Radchenko and Wheeler obtain solutions by proving the quadratic-sum identity for their finite quantum dilogarithms and then passing to Izumi's subsystem. The near-group TPE give the quadratic-sum equation directly, so that the last conversion is unnecessary for the categorical construction. See *Real quadratic fields and finite quantum dilogarithms I*, [arXiv:2609.21892](https://arxiv.org/abs/2609.21892).

These constructions motivate the formulation; this repository implements the general localization and spectrum criteria described in [GUIDE.md](GUIDE.md).

## Optional independent regeneration

The checked certificates can also be regenerated using existing mathematical software:

```sh
python3 optional/regenerate_f210.py       # Singular: generate the lifting certificate
python3 optional/check_generated.py      # Python: verify the regenerated certificate

gap -q optional/check_small_groups.g     # group search through order 128
gap -q optional/check_simple_groups.g    # historical simple-group checks; requires CTblLib
gap -q census/regenerate_group_rings.g   # regenerate four group tensors
```

The F210 generator writes to `generated-certificates/`. The recorded regeneration used Singular 4.3.2 and GAP 4.12.1.

## Conventions and provenance

Python fusion tensors and witness tuples use **zero-based** indices. Subtract 1 from each printed tuple entry in Propositions 6.9 and 6.10. Theorem 7.9 already uses zero-based indices. The fusion matrix `N[i]` has row index `j` and column index `k`. `data/manifest-sha256.json` records hashes of the frozen input data and certificates. Third-party census PDFs retain their original rights; their inclusion and provenance do not assert a new license for them.
