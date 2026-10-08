# Triangular prism equations and categorification

Computational material for **Triangular prism equations and categorification**, by Zhengwei Liu, Sebastien Palcoux, Yunxiang Ren, and Gert Vercleyen.

The programs apply the paper's **zero-spectrum, one-spectrum, and localization criteria to user-supplied fusion rings**. F210 and F660 are worked applications. Start with [GUIDE.md](GUIDE.md) for the input format, assumptions, examples, and connection to the manuscript. [CERTIFICATE.md](CERTIFICATE.md) explains the exact characteristic-zero certificate for F210.

For a new ring, the workflow is: search the spectrum criteria using the fusion rules alone; when the localization hypotheses hold, generate small necessary systems with justified or symbolic categorical dimensions; then add the coupling equations for compatible localizations. The last step can produce an obstruction even when the individual local systems have solutions. The guide separates these steps and explains what each output establishes.

## Apply the criteria to another ring

Use Python 3.10 or later with its standard library. Supply exact fusion coefficients in JSON, with `N[i][j][k]` the coefficient of basis element `k` in the product `i*j`, and unit index `0`.

```sh
python3 criteria.py spectrum my-ring.json --criterion both --mode first
python3 criteria.py localize examples/localization_fibonacci.json --output system.json
```

The spectrum search accepts noncommutative rings and reports every numerical hypothesis of each returned witness. Localization additionally needs a self-dual center, suitable multiplicity-free support, odd-multiplicity data, and categorical dimensions (specified exactly or left symbolic). The guide explains these choices and how to couple localizations. A spectrum witness is an obstruction; generating localization equations is a first step toward checking their consistency. Absence of an obstruction does not prove categorifiability.

## Verify the paper's applications

Run without Python's optimization flag, since the frozen certificate checks use assertions:

```sh
python3 verify.py
```

This command uses exact arithmetic in the Python standard library. Its largest input is a 24 MB compressed integer certificate, which allows direct checking of identities that were more expensive to discover. It checks the displayed fusion rules, both F210 localizations and their coupling, the F210 character table and positive-characteristic identities, the F660 witness, the rank-six spectrum examples, and the stated census exclusions. General-purpose implementations also have tests on known categorifiable rings, including a noncommutative group ring and Fibonacci fusion rules.

[PIVOTAL_CONVENTIONS.md](PIVOTAL_CONVENTIONS.md) explains the graphical normalization and its relation to the localization dimensions. Exact pointed-category regression tests check the single internal pivotal factor and the right contraction-dual normalization.

The census-based conclusion concerns nonpointed, unitary, 1-Frobenius, simple integral fusion categories of rank at most eight and Frobenius–Perron dimension at most 20000, up to Grothendieck equivalence. Completeness of the underlying census is an input from the cited classification. See [census/README.md](census/README.md) and its provenance file for the exact checks and the separate argument at dimension 20000.

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

All indices start at **0**. The fusion matrix `N[i]` has row index `j` and column index `k`. `data/manifest-sha256.json` records hashes of the frozen input data and certificates. Third-party census PDFs retain their original rights; their inclusion and provenance do not assert a new license for them.
