# Imported census and exact checks

Run `python census/verify_census.py` from the repository root. Verification uses Python 3.10 or later and its standard library.

This directory documents the external enumeration used in the rank-at-most-eight classification consequence. It contains the 33 rings of rank at most eight in the Bruns–Palcoux census, their exact integer fusion matrices and dimension vectors, and exact primary-3 exclusion witnesses for 27 of them.

## What is verified here

`verify_census.py` checks all 33 inputs for nonnegative integer coefficients, both unit identities, duality and its involution, Frobenius reciprocity, associativity, commutativity, the positive dimension character, the dimension bound and 1-Frobenius divisibility, and simplicity. Positivity of the common dimension eigenvector identifies its entries as the Frobenius–Perron dimensions. The script also checks that the two matrices used in this article agree entry by entry with the corresponding census inputs.

For each of 27 exclusions, `primary3_negative_witnesses.json` contains an integer vector `q` and an exact negative integer. With the row convention `M[i][j][k] = N[i,j]^k`, dimension vector `d`, and `L = lcm(d)`, the verifier evaluates

    q^T [sum_i (L/d_i) (M_i tensor M_i tensor M_i)] q < 0.

This is an exact certificate that the primary 3-matrix is not positive semidefinite. The primary criterion therefore excludes a unitary categorification; see Huang–Liu–Palcoux–Wu, *Complete Positivity of Comultiplication and Primary Criteria for Unitary Categorification*, IMRN 2024, 817–860, Theorem 3.18, DOI `10.1093/imrn/rnad214`. Integer tensor contractions perform the verification without floating point arithmetic. The discovery of the vectors used numerical eigenspaces followed by rounding; the independent exact verification does not depend on that discovery procedure.

The six inputs left by these 27 exclusions are:

| Input ID | Rank | Global dimension | Identification in the imported census |
|---|---:|---:|---|
| `S60_0` | 5 | 60 | Character ring of A5 = PSL(2,4) |
| `S168_0` | 6 | 168 | Character ring of PSL(2,7) |
| `S210_1` | 7 | 210 | F210, excluded separately in the article |
| `S360_1` | 7 | 360 | Character ring of A6 = PSL(2,9) |
| `S660_5` | 8 | 660 | F660, excluded separately in the article |
| `S660_11` | 8 | 660 | Character ring of PSL(2,11) |

The identifiers retain the source list name and its zero-based position. No claim of positive semidefiniteness for these six matrices is needed for the exclusion argument.

## Enumeration provenance and limits

The underlying exhaustive enumeration is an imported result. This verifier does **not** independently re-enumerate all fusion rings or prove that the imported list is exhaustive. Nor does checking one representative establish that the census contains no pair of isomorphic entries. Those assertions retain their stated dependence on the Bruns–Palcoux census.

The four character-ring identifications were checked independently using GAP 4.12.1: construct `AlternatingGroup(5)`, `PSL(2,7)`, `AlternatingGroup(6)`, and `PSL(2,11)`, compute their irreducible characters, and compute every tensor coefficient by the character inner product. `group_character_rings.json` freezes the resulting integer tensors and explicit basis bijections. The Python verifier checks all coefficients under these bijections. The optional command `gap -q census/regenerate_group_rings.g > regenerated_group_rings.json` regenerates the group tensors from the named groups; a GAP installation is needed only for this independent regeneration. The four frozen tensors matched the census tensors in their displayed order. Regenerated character orderings may depend on the GAP version, so dimension-preserving basis permutations may be needed to compare a later version.

The original documents are frozen under `source/` and their source commit, Git blob hashes, and SHA-256 digests are recorded in `PROVENANCE.json`. The source is the `SimpleIntegralPaper` directory of `sebastienpalcoux/Fusion-Categories`, commit `06f3a37b086c0591e086a3aed9a9458f4270165f`. Its summary is dated August 2025. The summary gives 505 rings across all rank-dependent bounds; only the 33 relevant rings are extracted here. The matrices were extracted from the data PDF using `pdftotext -layout` and parsed with `ast.literal_eval`, not executed as code. The original PDF remains available for direct comparison.

The August 2025 census summary states its rank-dependent Frobenius–Perron dimension bounds with ≤. The manuscript follows that convention. The headers in the archived data file instead use <. For rank eight the stated bound is 20000. Independently, a fusion category over the complex numbers of dimension 20000 = 2⁵ · 5⁴ is solvable and therefore contains a nontrivial invertible object (Etingof–Nikshych–Ostrik, *Weakly group-theoretical and solvable fusion categories*, Theorem 1.6 and Proposition 4.5(iv)). If the category is simple, the fusion subcategory generated by this object must be the whole category, which is then pointed. Thus the endpoint contributes no nonpointed simple fusion category.

The manuscript’s census table begins at rank five. In rank four, Dong–Palcoux–Plessis, *Classification of vanishing fusion rings of rank four*, Theorem 1.7, exclude simple integral 1-Frobenius fusion rings in every dimension. [arXiv:2610.22000](https://arxiv.org/abs/2610.22000)

## Errata in the archived source PDFs

In `source/1FrobSimpleIntegralFusionData.pdf`, page 53, the header for the 48 dimension-3366 tensors gives an incorrect type. The tensors have the dimension vector `(1,11,11,17,17,18,18,18,22,33)`, as correctly stated in the separate summary.

In `source/1FrobSimpleIntegral.pdf`, page 1, the dimension-1260 factorization should read `2^2 * 3^2 * 5 * 7`. Its “following pages” means the separate full-data PDF. The strict/inclusive bound difference is discussed above.

The preliminary list in Liu–Palcoux–Wu contains two rank-eight, dimension-1320 tensors representing the same based ring. Exchanging the two degree-15 basis elements (zero-based indices 5 and 6) carries the second tensor to the first. The present census retains one representative, `S1320_0`.

## Files

- `rank_at_most_8.json`: all 33 relevant imported inputs, with exact dimension vectors.
- `primary3_negative_witnesses.json`: 27 exact exclusion certificates.
- `verify_census.py`: a standard-library verifier.
- `group_character_rings.json`: independent GAP character-ring tensors and basis bijections.
- `regenerate_group_rings.g`: optional group-character regeneration using GAP.
- `PROVENANCE.json`: frozen source identifiers and integrity digests.
- `source/1FrobSimpleIntegral.pdf`: original summary.
- `source/1FrobSimpleIntegralFusionData.pdf`: original full 505-ring data document.
