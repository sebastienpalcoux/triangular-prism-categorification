# Triangular prism equations and categorification

Computational material for **Triangular prism equations and categorification**, by Zhengwei Liu, Sebastien Palcoux, Yunxiang Ren, and Gert Vercleyen.

Start with [GUIDE.md](GUIDE.md) for the connection between the manuscript and the files. [CERTIFICATE.md](CERTIFICATE.md) explains the characteristic-zero certificate for F210.

## Verify the results

Use Python 3.10 or later, without optimization (`-O`):

```sh
python3 verify.py
```

No packages, SageMath, Singular, GAP, network connection, or floating-point arithmetic are required for this command. A clean run with Python 3.12.14 took approximately 7–13 seconds on the verification machine. The largest file is a 24 MB compressed exact integer certificate; it replaces a potentially expensive algebraic search with direct arithmetic checks.

The command checks:

- The four displayed fusion tensors, their axioms, and the stated integral dimension vectors.
- The two localized F210 systems, the rescaling and linear eliminations, 28 exact ideal-membership identities, and the final nonsingularity certificate.
- The F210 formal codegrees and the ring identities that settle the exceptional positive characteristics.
- Every hypothesis of the displayed F660 zero-spectrum witness.
- Exhaustive zero/one-spectrum checks for the two rank-six examples.
- The 33 supplied census inputs, 27 exact primary-3 obstruction witnesses, and the four remaining group character-ring identifications.

The completeness of the underlying census is imported from the cited classification; the code does not independently repeat that enumeration. See [census/README.md](census/README.md) and its provenance file.

## Optional independent regeneration

Existing mathematical software is used only to generate certificates or access finite-group databases:

```sh
python3 optional/regenerate_f210.py       # requires Singular
python3 optional/check_generated.py      # checks its output with standard Python

gap -q optional/check_small_groups.g     # optional exhaustive group check through order 128
gap -q census/regenerate_group_rings.g   # optional regeneration of four group tensors
```

The F210 generator writes to `generated-certificates/`, keeping the frozen certificate separate. Singular 4.3.2 and GAP 4.12.1 were used in the audit. No SageMath installation is needed. Reimplementing the finite-group databases in Python would add complexity without improving the verification.

## Conventions and provenance

All indices in the data and code start at **0**, and `0` denotes the unit. The tensor entry `N[i][j][k]` is the coefficient of `b_k` in `b_i b_j`. Thus fusion matrix `N[i]` has row index `j` and column index `k`.

`data/manifest-sha256.json` records hashes of the frozen input data and certificates. Third-party census PDFs retain their original rights; their inclusion and provenance do not assert a new license for them. This repository contains computational material and public source data, not editorial correspondence.
