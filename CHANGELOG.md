## [Unreleased]
### Changed
- The fit estimates the log transcription quotients x_gc directly instead of the log fold changes, following
  Section 2 of E. van Nimwegen's note "Getting proper error-bars in Sanity" (see `docs/output_logic.md`).
  The selected gene variance v_g is unchanged (`docs/likelihood_equivalence.md`), and so is
  delta_gc = x*_gc - mu_g; mu_g is now the mean of x*_gc over the cells (a shift of log(n_g) - psi(n_g)), and the
  error bars on x_gc and mu_g come from the full posterior covariance over the cells (note Eqs. 41, 43).
- The asymmetric error bar for zero-count cells is gone; every cell gets the posterior error bar. MARG is several
  times faster as a result.
- `ltq_error_bars` no longer combines `d_mu` and `d_delta`: it is the posterior error bar itself, which already
  includes the uncertainty on mu_g. `d_delta` holds the same values.

### Added
- The measurement LTQs x^m_gc (the posterior with the prior divided out, note Eq. 44) with error bars
  1/sqrt(n_g f*_gc) (Eq. 46): `measured_ltq` and `measured_ltq_error_bars`, written with `-e 1` or selected with the
  new `--ltq-estimate` (`-ltq`) option (`posterior`, `measured` or `both`).
- Bonsai output mode keeps its file names but holds the measurement estimates by default
  (`--ltq-estimate posterior` gives the previous content); `sanity_command.txt` records `# LTQ estimate:`.
- `tests/compare_section2.py`: an independent NumPy implementation of the fit's equations that the output is
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
