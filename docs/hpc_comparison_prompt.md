# Prompt for the HPC agent: compare three versions of Sanity on simulated data

You are working on the CBE cluster (Slurm). Your task is to compare three versions of Sanity on simulated
single-cell count matrices with known ground truth, and to write a short report. The code is in my fork
https://github.com/dhdegroot/Sanity_daansfork (or the clone I point you to), branch `clean-up-meas`. Do not change
the code. Ask me before running anything on the data.

The three versions:

| version | source | command-line difference | what its LTQ files hold |
|:--|:--|:--|:--|
| **A** previous | commit `01df9080` (the last commit before the new fit; it is on the history of `clean-up-meas` and is the tip of branch `clean-up` of the upstream repository jmbreda/Sanity) | none | the previous implementation: `mu_g + delta_gc` with the SI 2021 error bars |
| **B** new, posterior-mode | branch `clean-up-meas` | default (`--measured-expansion posterior-mode`) | posterior LTQs `x*_gc` with posterior error bars, and measurement LTQs `x^m_gc` from the note's Eq. 44/46 (expansion around the posterior mode) |
| **C** new, likelihood-maximum | same binary as B | `--measured-expansion likelihood-maximum` | the same posterior LTQs; measurement LTQs `log(n_gc/N_c) ± 1/sqrt(n_gc)` for nonzero counts, Eq. 44/46 for zero counts |

B and C are one binary run with two settings and differ only in the two `measured_ltq*` files. Background (all on
`clean-up-meas`): `docs/output_logic.md` lists the formula behind every output file,
`docs/measured_ltq_linearisation.md` explains why Eq. 44 fails for surprising counts and why C exists,
`docs/measured_expansion_plan.md` describes the option.

## 1 Build the two binaries and record the commits

Clone the fork once and make two checkouts, one per binary:

```bash
git clone https://github.com/dhdegroot/Sanity_daansfork.git sanity-new
cd sanity-new && git checkout clean-up-meas
git worktree add ../sanity-previous 01df9080
```

Build each with GCC and OpenMP (`module load` a recent GCC first):

```bash
cd src && make Sanity
```

The binary is `bin/Sanity` in each checkout. Record `git rev-parse HEAD` of both checkouts and `bin/Sanity -v` in
the report. In `sanity-new`, check that `bin/Sanity -h` lists `--measured-expansion`; if it does not, you have the
wrong commit: stop and tell me. Then run `python3 tests/compare_section2.py` there once (needs NumPy); it must
print PASSED.

## 2 How to find the data

1. The paths of the simulated count matrices (TSV with a header line of cell names and one gene per row, or
   Matrix Market `.mtx` with `-mtx_genes` and `-mtx_cells` files; either may be gzipped) can be found in `/groups/degroot/lab/data_public/scrnaseq/simulated_datasets` with `simulated_binary_10_gens_samplingnoise_seed_1231` as an example, but you can pick 3 more datasets to test from that folder. Note that they contain a Gene_table.txt matrix, which you need and the ground truth is in delta_true.txt and variance_true.txt
2. Ask me how much compute I want to spend (number of threads per job, partition, time limit) and where to put the output.
3. Make sure to start off slurm on all jobs. I need to leave soon and will then close down this agent. Write a follow-up prompt in the docs/ folder of this repo.

## 3 How Sanity is run

All versions take the same command line; C adds one option. Run every dataset with every variance method and
full-precision output:

```bash
# A (sanity-previous/bin/Sanity) and B (sanity-new/bin/Sanity)
bin/Sanity -f <counts.tsv> -d <out_dir>/ -n <threads> -v_m <METHOD> -e 1 -npy
# C (sanity-new/bin/Sanity)
bin/Sanity -f <counts.tsv> -d <out_dir>/ -n <threads> -v_m <METHOD> -e 1 -npy --measured-expansion likelihood-maximum
```

with `METHOD` in `MAP` (the default), `MLE`, `EAP`, `MARG`. For now, only try `EAP` and `MAP`.
Each run goes in its own output directory. That is
3 versions × 2 methods per dataset; submit them as Slurm jobs (one job per dataset × version × method, or a job
array), with `--cpus-per-task` equal to `-n`. Record for every run the wall time and peak memory from `sacct`
(`Elapsed`, `MaxRSS`) and from `/usr/bin/time -v` if available. MARG is the slowest; in A it is several times slower
than in B and C. B and C cost the same (the fit is identical; only the second pass differs).

Genes with zero total count are dropped by Sanity; `geneID.txt` lists the genes that were kept, in output order.
Every numerical output is float64 NPY (`numpy.load(path, allow_pickle=False)`), one row per gene in the order of
`geneID.txt`, matrix columns in the order of `cellID.txt`. Error bars are standard deviations, not variances.

### Output of A (`-e 1 -npy`)

| file | shape | content |
|:--|:--|:--|
| `log_transcription_quotients.npy` | genes × cells | LTQ estimate `mu_g + delta_gc` |
| `ltq_error_bars.npy` | genes × cells | `sqrt(d_mu_g^2 + d_delta_gc^2)` |
| `delta.npy` | genes × cells | log fold change `delta_gc` from the gene mean |
| `d_delta.npy` | genes × cells | error bar on `delta_gc` (asymmetric zero-count variant for cells with count 0) |
| `mu.npy` | genes × 1 | gene mean `mu_g = psi(n_g) - q_g` (digamma) |
| `d_mu.npy` | genes × 1 | `sqrt(trigamma(n_g))` (MARG: plus the spread over the grid) |
| `variance.npy` | genes × 1 | the selected `v_g` (MARG: posterior mean of `v_g`) |
| `likelihood.npy` | genes × bins | posterior over the variance grid |
| `likelihood_grid.npy` | bins × 1 | the variance grid (160 points from 0.001 to 50, log-spaced) |

### Output of B and C (`-e 1 -npy`)

| file | shape | content |
|:--|:--|:--|
| `log_transcription_quotients.npy` | genes × cells | posterior LTQ `x*_gc` (MARG: averaged over the grid) |
| `ltq_error_bars.npy` | genes × cells | posterior error bar `eps_gc` (includes the uncertainty on `mu_g`) |
| `measured_ltq.npy` | genes × cells | B: `x^m_gc` of Eq. 44 for every cell. C: `log(n_gc/N_c)` where `n_gc > 0`, Eq. 44 where `n_gc = 0` |
| `measured_ltq_error_bars.npy` | genes × cells | B: `1/sqrt(n_g f*_gc)` (Eq. 46), can be very large. C: `1/sqrt(n_gc)` where `n_gc > 0`, Eq. 46 where `n_gc = 0` |
| `delta.npy` | genes × cells | `x*_gc - mu_g` |
| `d_delta.npy` | genes × cells | identical to `ltq_error_bars.npy` |
| `mu.npy` | genes × 1 | `mu_g = mean over cells of x*_gc` |
| `d_mu.npy` | genes × 1 | `sigma_mu` (note Eq. 43; MARG: Eq. 50) |
| `variance.npy` | genes × 1 | the selected `v_g` (MARG: posterior mean) |
| `likelihood.npy`, `likelihood_grid.npy` | as above | as above |

`sanity_command.txt` in each output directory records the version, method, the command line and, for B and C,
`# LTQ estimate:` and `# Measured expansion:`; copy it into the report. For MARG, the zero-count entries of C's
measured files are Eq. 44 averaged over the grid, as in B.

## 4 Compare each version with the ground truth

For each dataset, version and method:

- **`v_g`**: the estimate versus the true `v_g`. Report the bias of `log(v_est / v_true)`, its spread, and the fraction
  of genes within a factor 2, split by expression level (bin the genes by their total count `n_g`, e.g. decades).
  Also note how many genes sit at the grid edges (0.001 or 50). B and C are identical here; report once.
- **`mu_g`**: estimate minus the true mean (`mu.npy`), by expression level. Note that A uses `psi(n_g) - q_g` and
  B/C use the mean of `x*_gc`, so a shift of `log(n_g) - psi(n_g)` (about `1/(2 n_g)`) between A and B/C is
  expected; compare each with the truth.
- **LTQs**, four estimates: A's `mu_g + delta_gc` (its `log_transcription_quotients.npy`); B's `x*_gc`
  (`log_transcription_quotients.npy`, identical in C); B's `x^m_gc` (`measured_ltq.npy`, Eq. 44); C's `x^m_gc`
  (`measured_ltq.npy`, likelihood maximum). Against the true `x_gc`: RMSE and bias, by expression level and
  separately for cells with zero and with nonzero counts; for the nonzero counts also by the count itself
  (`n_gc` = 1, 2, 3–4, 5–9, ≥ 10).
- **Error-bar calibration**: for each of the four estimates with its error bar (A's `ltq_error_bars.npy`, B's
  `ltq_error_bars.npy`, B's `measured_ltq_error_bars.npy`, C's `measured_ltq_error_bars.npy`), compute
  `z = (estimate - truth) / error bar`. Report the fraction of entries with `|z| < 1`, `< 2`, `< 3` (expected
  0.683, 0.954, 0.997 if calibrated) and the distribution of `z` (histogram or quantiles), split by expression level,
  by zero versus nonzero counts, and for the nonzero counts by `n_gc` as above. Do the same for `mu_g` with `d_mu.npy`
  and, if the truth has it, for the gene variance.
- **What to expect for the measurement LTQs.** B's `x^m_gc` is expected to be far from the truth for cells whose count
  greatly exceeds the expected count (`x^m - x* = n_gc/(n_g f*_gc) - 1`): quantify how many entries have
  `|x^m - x*| > 1`, `> 5`, how many have `x^m > 0` (a fraction above 1, impossible), and whether their error bars
  cover the truth. C's `x^m_gc` is `≤ 0` for every nonzero count by construction; its Gaussian is the Laplace
  approximation of a skewed likelihood, so some under-coverage at `n_gc = 1` and `2` (of order 10 %) is expected and
  should shrink with `n_gc`. Treat all of this as findings, not errors.

## 5 Compare the versions with each other

For the same dataset and method:

- **A versus B.** `variance.npy` should be identical for MAP, MLE and EAP (the selected grid point is the same by
  construction) and equal to about 1e-9 for MARG; `likelihood.npy` should agree to about 1e-8. Report the number of
  genes where the selected `v_g` differs, if any. `delta.npy` should agree to about 1e-8. `mu.npy` differs by
  `log(n_g) - psi(n_g)`; `ltq_error_bars.npy` and `d_delta.npy` differ on purpose (new posterior error bars versus
  the old ones, and no asymmetric zero-count bar any more). Show the ratio of B's to A's error bars by expression
  level and by zero versus nonzero counts.
- **B versus C.** Every file except the two `measured_ltq*` files must be bit-identical (`numpy.array_equal`); so
  must the zero-count entries of the two measured files. Report any difference as a bug. For the nonzero entries show
  the distribution of `x^m(C) - x^m(B)` and of the ratio of the error bars, by the expected count at the gene mean
  `N_c exp(mu_g)` (bins below 0.01, 0.01–0.1, 0.1–1, above 1) and by `n_gc`.
- Runtime and peak memory of the three versions per method.

## 6 Report

Write a short report (Markdown, with tables and figures as PNG files) with: the commit hashes and build details, the
datasets and which truth they contain, the comparisons of Sections 4 and 5, a recommendation on which of B and C
should become the default for the measurement LTQs with the evidence for it, and a clear list of what you could not
check (missing truth quantities, failed jobs, datasets you skipped) and why. Put the scripts you used next to it so
the analysis can be rerun. Use only `numpy`, `matplotlib` and the standard library unless I say otherwise, and load
the large matrices with `mmap_mode="r"` when memory is tight.
