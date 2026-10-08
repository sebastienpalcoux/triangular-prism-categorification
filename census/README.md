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

There is a boundary discrepancy in the original documents: the summary states a weak dimension bound, while the data-file header states strict bounds. The manuscript handles global dimension 20000 separately: since 20000 = 2^5 * 5^4, Etingof–Nikshych–Ostrik, *Weakly group-theoretical and solvable fusion categories*, Theorem 1.6, makes such a category solvable. By Proposition 4.5(iv) it has a nontrivial invertible object. A nonpointed simple fusion category cannot have such an object, because it generates a nontrivial proper pointed fusion subcategory. Thus that boundary contributes no categorification to the conclusion.

## Files

- `rank_at_most_8.json`: all 33 relevant imported inputs, with exact dimension vectors.
- `primary3_negative_witnesses.json`: 27 exact exclusion certificates.
- `verify_census.py`: a standard-library verifier.
- `group_character_rings.json`: independent GAP character-ring tensors and basis bijections.
- `regenerate_group_rings.g`: optional group-character regeneration using GAP.
- `PROVENANCE.json`: frozen source identifiers and integrity digests.
- `source/1FrobSimpleIntegral.pdf`: original summary.
- `source/1FrobSimpleIntegralFusionData.pdf`: original full 505-ring data document.
