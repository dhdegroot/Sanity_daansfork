# Tests

This directory contains a simple test for Sanity.

## Overview

The test runs Sanity once for each `-v_m` method (MAP, EAP, MLE, MARG), each time with `-e 1`, on a small count matrix (`count_table.tsv`). Each run's output is compared against the reference results in the matching `default_output_<METHOD>` folder.

The comparison for each method can come out as:
- Identical results
- Results that differ but are numerically close (checked with `numpy.isclose`, using its default tolerances)
- Results that differ significantly

## Requirements

Python script `compare.py` requires the `numpy` library (`pip3 install numpy` if you don't have it).

## Running the Test

1. **Compile Sanity** in the `src` project directory
2. **Navigate to the tests directory**:
    ```bash
    cd tests
    ```
3. **Run the test script**
    ```bash
    python3 compare.py
    ```
    Details on differences are saved in `compare.log`, which is truncated at the start of each run and covers all four methods. It is not tracked in git.

Note: `compare.py` exits with code 0 if every method passed, and 1 if any method FAILED, so it can be used in a script or CI. A method whose results differ but are within the `numpy.isclose` tolerances is reported as "PASSED with acceptable differences" and counts as a pass.

## Matrix Market reader

After compiling, run `python3 compare_mtx.py`. It compares complete Bonsai NPY
outputs and identifier order against TSV input using the default thread setting
and `-n 1`, `-n 2`, and `-n 12`. It covers the parallel integral-count path,
fractional counts, gzip input, and removal of an empty gene. Temporary inputs
and outputs are removed automatically.

Run `python3 compare_mtx_index.py` to check every index entry and total,
partition/buffer boundaries, the thread limit, reader fallbacks and malformed
input rejection. It compiles a small test driver using `CXX` (default `g++`),
requires NumPy and removes its temporary files automatically.

## Per-gene input rows

Run `python3 compare_rows.py` to check exact counts for integer/fractional MTX,
scientific notation, large values, duplicates, gzip, TSV, reverse reads and
buffer boundaries. It also checks concurrent reader instances, limited thread
teams, stale row indexes, malformed entries and truncated rows. The script
uses NumPy, compiles a temporary driver with `CXX` (default `g++`), and removes
its temporary files.

## Section 2 reference

Run `python3 compare_section2.py` after compiling. It is an independent NumPy
implementation of Section 2 of E. van Nimwegen's note "Getting proper
error-bars in Sanity" (the formulas Sanity's fit implements): q_g is solved by
bisection with its own Lambert W iteration, L(v) comes from Eqs. 36–38, the
estimates at the selected grid point from Eqs. 32, 33, 39, 41, 43, 44 and 46,
and the MARG averages from Eqs. 47–50. The script compares the `-e 1 -npy`
output of all four `-v_m` methods with these, on a synthetic 24 x 60 matrix or
on a TSV given with `--counts`, and checks the identities the note relies on
(`sum_c f*_c = 1`, `<x*> = log(n_g) - q_g`, Eqs. 40, 41 and 43 against the
numerically inverted precision matrix of Eq. 35, and the constant offset
between Eqs. 36–38 and the previous likelihood, see
`docs/likelihood_equivalence.md`). It exits with 1 if any check fails.

## Inference changes

Compile the previous version separately, then run:

```bash
python3 compare_inference.py --reference /path/to/previous/bin/Sanity
```

An optional `--candidate` selects a different candidate executable. This checks
all four variance methods on TSV and fractional MTX inputs, a broad count range,
empty genes and a custom variance grid.
Extended text, float64 arrays, Bonsai arrays and identifier order must match
exactly. Command metadata is excluded because the executable paths differ.
A reference built before `-npy` wrote every selected result as NPY (before the
output clean-up) writes a different set of `-npy` files, so those comparisons fail.
Inputs and outputs are generated in a temporary directory and removed.

The inference comparison also covers cached fitting constants, including
custom variance grids.
