## [3.0.0] - 2026-11-16
### Changed
- Every output file is named after what it holds, whatever the options of the run: `ltqs_posterior`,
  `ltqs_posterior_error_bars`, `ltqs_measurement_approx`, `ltqs_measurement_approx_error_bars`, `gene_means`,
  `gene_means_error_bars`, `gene_variances`, `gene_variances_likelihood` (and `gene_variances_grid.npy` with
  `-npy -e 1`), replacing `log_transcription_quotients`, `ltq_error_bars`, `measured_ltq`, `measured_ltq_error_bars`,
  `mu`, `d_mu`, `variance`, `likelihood` and `likelihood_grid`. `delta` and `d_delta` are no longer written: they
  were `ltqs_posterior` minus `gene_means` and a copy of `ltqs_posterior_error_bars`. The per-gene vectors,
  `geneID.txt` and `cellID.txt` are always written.
- `-npy` writes the per-cell matrices and the likelihood table as float32, which carries the six decimals of the
  text output, and the per-gene vectors as float64, in every mode.
- `--bonsai-output-mode` is the preset `-npy -v_m EAP --ltq-estimate measured` and nothing else: no fixed file
  set, no override of `-e` or `-gz`, no folder check, and EAP instead of MAP (an explicit other `-v_m` is an
  error). `sanity_command.txt` keeps `# Output mode: bonsai` and `# LTQ estimate:`; the prose `# Output:` line is
  gone, since the file names say what they hold.
- `Sanity_distance` and `Sanity_gene_correlation` read the new file names and take the prior-free values from
  `ltqs_measurement_approx` instead of dividing the prior out themselves (SI Eqs. 62, 63), which was inexact once
  the posterior error bars carried the uncertainty on the gene mean.
- The fit estimates the log transcription quotients x_gc directly instead of the log fold changes.
  The selected gene variance v_g is unchanged, and so is
  delta_gc = x*_gc - mu_g; mu_g is now the mean of x*_gc over the cells (a shift of log(n_g) - psi(n_g)), and the
  error bars on x_gc and mu_g come from the full posterior covariance over the cells.
- The asymmetric error bar for zero-count cells is gone; every cell gets the posterior error bar. MARG is several
  times faster as a result.
- The posterior error bar on x_gc no longer combines the error bars on mu_g and delta_gc of 2.0: it is the
  posterior standard deviation itself, which already includes the uncertainty on mu_g.

### Added
- The measurement LTQs x^m_gc (the posterior with the prior divided out) with error bars
  1/sqrt(n_g f*_gc): `ltqs_measurement_approx` and `ltqs_measurement_approx_error_bars`, written with
  `-e 1`, in Bonsai output mode, or selected with the new `--ltq-estimate` (`-ltq`) option (`posterior`, `measured`
  or `both`); `sanity_command.txt` records `# LTQ estimate:`.
- `tests/compare_numpy_fit.py`: an independent NumPy implementation of the fit's equations that the output is
  checked against.

## [2.0.0] - 2026-07-29
### Added
- New methods for calculating v_g (`-v_m` option):
    - MAP (**default**) maximum a posteriori estimate of v_g
    - EAP expected a posteriori estimate of v_g
    - MLE maximum likelihood estimate of v_g
    - MARG marginal likelihood estimate of v_g
- Lazy reading of the input file, and lazy output writing which significantly reduces memory usage for large datasets.
- Support for gzipped input files (`.gz` extension).
- Timestamped run command is now saved by default in the output directory (`sanity_command.txt`).

### Changed
- By default, *Sanity* will use the MAP method for calculating v_g, if not specified otherwise.
- No `*_vmax` output files. For all methods the output files are named the same.
- Text input is now strictly tab-separated.
Commas and spaces used to be accepted as separators in some parts of the reader but not in others, so a comma- or space-separated file could be read inconsistently.
They are now part of a field everywhere, which also means gene and cell names may contain them.
A file with no tabs at all is rejected with an error instead of being processed as if it had zero cells.
- Fixed a sign error in the trigamma asymptotic series, which slightly affected the error bars on mu (and hence `d_mu.txt` and `ltq_error_bars.txt`) for genes with very low total counts, by under 1e-4 relative.
