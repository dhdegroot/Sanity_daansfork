# Prompt for the HPC agent: compare Sanity `clean-up` with `new-direct-ltqs` on simulated data

You are working on the CBE cluster (Slurm). Your task is to compare two branches of the Sanity repository
(https://github.com/jmbreda/Sanity, or the clone I point you to) on simulated single-cell count matrices with known
ground truth, and to write a short report. Do not change either branch. Ask me before running anything on the data.

## 1 Build both branches and record the commits

Clone the repository twice (or use two worktrees), one per branch:

- `clean-up`: the previous implementation (log fold changes, SI 2021 formulas).
- `new-direct-ltqs`: the new implementation (direct LTQs with posterior error bars, Section 2 of E. van Nimwegen's
  note "Getting proper error-bars in Sanity"). Its `docs/output_logic.md` lists the formula behind every output file,
  and `docs/measured_ltq_linearisation.md` explains a known property of the measurement LTQs.

Build each with GCC and OpenMP (`module load` a recent GCC first; Apple clang is not an issue here):

```bash
cd src && make Sanity
```

The binary is `bin/Sanity`. Record `git rev-parse HEAD` of each checkout and `bin/Sanity -v` in the report.
Run `python3 tests/compare_section2.py` in the `new-direct-ltqs` checkout once (needs NumPy); it must print PASSED.

## 2 Ask me before running

Before any run, ask me for:

1. The paths of the simulated count matrices (TSV with a header line of cell names and one gene per row, or
   Matrix Market `.mtx` with `-mtx_genes` and `-mtx_cells` files; both may be gzipped).
2. The paths of their ground truth, and which quantities it contains and in which layout. Ask explicitly whether it
   has the true LTQs `x_gc` (log fraction of the cell's mRNA from gene `g`), the true gene variances `v_g`, the true
   gene means `mu_g`, and how genes and cells are identified (names or order). Check that every dataset's gene and
   cell order matches the truth before using it.
3. How much compute I want to spend (number of threads per job, partition, time limit) and where to put the output.

## 3 How Sanity is run

Both branches have the same command line. Run every dataset with every variance method and full-precision output:

```bash
bin/Sanity -f <counts.tsv> -d <out_dir>/ -n <threads> -v_m <METHOD> -e 1 -npy
```

with `METHOD` in `MAP` (the default), `MLE`, `EAP`, `MARG`. Each run goes in its own output directory. Submit the runs as
Slurm jobs (one job per dataset × branch × method is fine, or a job array), with `--cpus-per-task` equal to `-n`.
Record for every run the wall time and peak memory from `sacct` (`Elapsed`, `MaxRSS`) and from `/usr/bin/time -v`
if available. MARG is the slowest; on `clean-up` it is several times slower than on `new-direct-ltqs`.

Genes with zero total count are dropped by Sanity; `geneID.txt` lists the genes that were kept, in output order.
Every numerical output is float64 NPY (`numpy.load(path, allow_pickle=False)`), one row per gene in the order of
`geneID.txt`, matrix columns in the order of `cellID.txt`. Error bars are standard deviations, not variances.

### Output of `clean-up` (`-e 1 -npy`)

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

### Output of `new-direct-ltqs` (`-e 1 -npy`)

| file | shape | content |
|:--|:--|:--|
| `log_transcription_quotients.npy` | genes × cells | posterior LTQ `x*_gc` (MARG: averaged over the grid) |
| `ltq_error_bars.npy` | genes × cells | posterior error bar `eps_gc` (includes the uncertainty on `mu_g`) |
| `measured_ltq.npy` | genes × cells | measurement LTQ `x^m_gc` (prior divided out, note Eq. 44) |
| `measured_ltq_error_bars.npy` | genes × cells | `1/sqrt(n_g f*_gc)` (note Eq. 46); can be very large |
| `delta.npy` | genes × cells | `x*_gc - mu_g` |
| `d_delta.npy` | genes × cells | identical to `ltq_error_bars.npy` |
| `mu.npy` | genes × 1 | `mu_g = mean over cells of x*_gc` |
| `d_mu.npy` | genes × 1 | `sigma_mu` (note Eq. 43; MARG: Eq. 50) |
| `variance.npy` | genes × 1 | the selected `v_g` (MARG: posterior mean) |
| `likelihood.npy`, `likelihood_grid.npy` | as above | as above |

`sanity_command.txt` in each output directory records the version, method, `# LTQ estimate:` (new branch only) and
the command line; copy it into the report.

## 4 Compare each branch with the ground truth

For each dataset, branch and method:

- **`v_g`**: the estimate versus the true `v_g`. Report the bias of `log(v_est / v_true)`, its spread, and the fraction
  of genes within a factor 2, split by expression level (bin the genes by their total count `n_g`, e.g. decades).
  Also note how many genes sit at the grid edges (0.001 or 50).
- **`mu_g`**: estimate minus the true mean (`mu.npy`), by expression level. Note that `clean-up` uses `psi(n_g) - q_g`
  and `new-direct-ltqs` uses the mean of `x*_gc`, so a shift of `log(n_g) - psi(n_g)` (about `1/(2 n_g)`) between the
  branches is expected; compare both with the truth.
- **LTQs**: on `new-direct-ltqs`, `x*_gc` (`log_transcription_quotients.npy`) and `x^m_gc` (`measured_ltq.npy`); on
  `clean-up`, `mu_g + delta_gc` (which is `log_transcription_quotients.npy` there). Against the true `x_gc`: RMSE and
  bias, by expression level and separately for cells with zero and with nonzero counts.
- **Error-bar calibration**: for each LTQ estimate with its error bar (`ltq_error_bars.npy`,
  `measured_ltq_error_bars.npy`, and `clean-up`'s `ltq_error_bars.npy`), compute `z = (estimate - truth) / error bar`.
  Report the fraction of entries with `|z| < 1`, `< 2`, `< 3` (expected 0.683, 0.954, 0.997 if calibrated) and the
  distribution of `z` (histogram or quantiles), split by expression level and by zero versus nonzero counts. Do the same
  for `mu_g` with `d_mu.npy` and, if the truth has it, for the gene variance.
- The measurement LTQs `x^m_gc` of the new branch are expected to be far from the truth for cells whose count greatly
  exceeds the expected count (`x^m - x* = n_gc/(n_g f*_gc) - 1`); quantify how many entries have `|x^m - x*| > 5` and
  whether their error bars cover the truth. Treat this as a finding, not an error.

## 5 Compare the two branches with each other

For the same dataset and method:

- `variance.npy` should be identical for MAP, MLE and EAP (the selected grid point is the same by construction) and
  equal to about 1e-9 for MARG; `likelihood.npy` should agree to about 1e-8. Report the number of genes where the
  selected `v_g` differs, if any.
- `delta.npy` should agree to about 1e-8.
- `mu.npy` differs by `log(n_g) - psi(n_g)`; `ltq_error_bars.npy` and `d_delta.npy` differ on purpose (new posterior
  error bars versus the old ones, and no asymmetric zero-count bar any more). Show the ratio of new to old error bars by
  expression level and by zero versus nonzero counts.
- Runtime and peak memory of both branches per method.

## 6 Report

Write a short report (Markdown, with tables and figures as PNG files) with: the commit hashes and build details, the
datasets and which truth they contain, the comparisons of Sections 4 and 5, and a clear list of what you could not
check (missing truth quantities, failed jobs, datasets you skipped) and why. Put the scripts you used next to it so
the analysis can be rerun. Use only `numpy`, `matplotlib` and the standard library unless I say otherwise, and load the
large matrices with `mmap_mode="r"` when memory is tight.
