# Plan: measurement LTQs expanded around the likelihood maximum

Branch `new-direct-ltqs`. Equation numbers refer to Section 2 of E. van Nimwegen's note; (B), (C), (E) to
`docs/measured_ltq_linearisation.md`. Implemented as planned (commits "Add --measured-expansion …", "Check both
measured expansions …", "Document --measured-expansion"); Section 4 is what the HPC comparison should confirm.

## 1 What changes

Today the measurement LTQ of every cell comes from the quadratic expansion of the Poisson log-likelihood around the
posterior mode $x^*_c$ (Eq. 44, 46):

$$
x^m_c = x^*_c + \frac{n_c - n f^*_c}{n f^*_c}, \qquad (\epsilon^m_c)^2 = \frac{1}{n f^*_c}.
$$

The alternative expands the Poisson log-likelihood of a **nonzero** count around its own maximum
$\hat x_c = \log(n_c/N_c)$, where the slope vanishes and the curvature is $-n_c$:

$$
n_c > 0: \quad x^m_c = \log(n_c / N_c), \qquad (\epsilon^m_c)^2 = \frac{1}{n_c}. \tag{E}
$$

**Zero counts keep Eq. 44/46**: $x^m_c = x^*_c - 1$, $(\epsilon^m_c)^2 = 1/(n f^*_c)$. Their likelihood
$e^{-N_c e^x}$ has no maximum, and the expansion at $x^*_c$ never blows up for them (the bar is wide when the zero is
uninformative, of order 1 when it is not; Section 6.2 of the linearisation note).

Everything else is untouched: $f^*_c$, $q_g$, $L(v)$, the selected $v_g$, $x^*_c$, $\epsilon_c$, $\mu_g$,
$\sigma_\mu$ and the posterior over the grid are identical whichever expansion is used. For a nonzero count (E) does
not depend on $v$, so under MARG the average over the grid returns (E) unchanged with zero spread; zero-count cells
are averaged over the grid as now.

## 2 Interface

A new option that chooses the expansion point of the measurement LTQs, orthogonal to `--ltq-estimate` (which
chooses whether the measured files are written at all):

| option | values | meaning |
|:--|:--|:--|
| `-mx`, `--measured-expansion` | `posterior-mode` (default) | Eq. 44/46 for every cell, as now |
| | `likelihood-maximum` | (E) for $n_c > 0$, Eq. 44/46 for $n_c = 0$ |

- Applies to `measured_ltq`, `measured_ltq_error_bars`, and to `delta.npy` / `d_delta.npy` in Bonsai output mode when
  `--ltq-estimate measured` (there `delta.npy` $= x^m_c - \mu_g$, i.e. $\log(n_c/N_c) - \mu_g$ for nonzero counts).
- `sanity_command.txt` records `# Measured expansion: posterior-mode|likelihood-maximum`; the Bonsai `# Output:` line
  names it too.
- The file names do not change and one run writes one variant. Getting both variants of one dataset costs two full
  runs (the fit is the whole cost; the second pass is negligible). Writing both variants from one run would mean two
  more output files, against the goal of a minimal output, so it is not proposed; the HPC comparison runs the branch
  twice instead.
- The default stays `posterior-mode` until the HPC comparison has decided. Switching the default afterwards is a
  one-line change plus regenerated `tests/default_output_*`; dropping the option altogether is then also possible.

## 3 Code changes, in commit order

Each commit builds and passes the applicable tests; `-n 1/4/8` stay bit-identical (the change is per cell, no new
reductions).

1. **`src/sanity_options.hpp/.cpp`** – `enum class MeasuredExpansion { PosteriorMode = 0, LikelihoodMaximum = 1 }`,
   a `MeasuredExpansion measured_expansion = MeasuredExpansion::PosteriorMode;` member of `SanityOptions`, the option
   table entry `{Option::MeasuredExpansionOption, {"-mx", "--measured-expansion", "--measured_expansion"}, true}`,
   `parse_measured_expansion` next to `parse_ltq_estimate`, `measured_expansion_name()`, and one usage line.
   No interaction with the other options (Bonsai mode accepts both values).
2. **`src/gene_fit.hpp/.cpp`** – `get_gene_expression_level` takes the enum (one more argument; `fit_fractions`
   unchanged). `cell_error_bars` gets the cell's count, its $\log N_c$ and the enum; for
   `count > 0 && expansion == LikelihoodMaximum` it returns `measured_ltq = log(count) - log_cell_total` and
   `measured_ltq_error_sq = 1.0 / count`, otherwise Eq. 44/46 as now. Both call sites
   (`estimate_at_selected_bin`, `marginalise_over_variance`) already loop over the cells with `umi_counts_c` and
   `log_cell_totals_c` in scope, so they pass `umi_counts_c[c]` and `log_cell_totals_c[c]`. No change to the West
   updates: for a value constant over the grid the weighted mean is that value and the spread is exactly zero.
   Header comment of `measured_ltq_c` / `measured_ltq_error_sq_c` updated with (E).
3. **`src/sanity_main.cpp`** – passes `options.measured_expansion` into the per-gene call (the only place the
   fit is invoked).
4. **`src/sanity_output.cpp`** – the `# Measured expansion:` line in `sanity_command.txt`; the Bonsai `# Output:`
   line; the text-file header of the two measured matrices names the expansion.
5. **`tests/compare_section2.py`** – the reference gains the (E) variant (`np.log(n_c / N_c)`, `1 / n_c` where
   `n_c > 0`, else the Eq. 44/46 values; for MARG the zero-count cells are averaged over the grid as the posterior-mode
   reference already does). The script runs Sanity once more with `-mx likelihood-maximum` per method and checks:
   the measured files against the reference (tolerance 1e-7 as now), the eight non-measured files bit-identical to
   the `posterior-mode` run, and the zero-count entries of the measured files identical between the two runs.
   `tests/README.md` describes the extra run. `tests/default_output_*` do not change (default unchanged);
   `compare.py` must still report identical output.
6. **Docs** – `README.md` (option list, usage text, the Output paragraph on the two kinds of estimates, the Bonsai
   section), `docs/output_logic.md` (the files table gets the two formulas for `measured_ltq*`, the Bonsai table a
   note), `docs/measured_ltq_linearisation.md` (last paragraph: option 1 is now selectable), `CHANGELOG.md`.
7. **Verification before the HPC run** – `tests/compare.py`, `compare_section2.py` (synthetic and
   `count_table.tsv`), the three reader tests, `-n 1/4/8` bit-identity for MAP and MARG with both flag values,
   `Sanity_distance` on `-e 1` output. On the 3000 × 3000 benchmark: no measurable runtime or memory difference
   is expected (one branch per cell in pass 2).

Size: about 60 lines of C++, 40 lines of Python, and the documentation; one day including the checks.

## 4 Expected behaviour, to be confirmed on the HPC data

- Nonzero entries: $x^m_c \le 0$ always (a fraction at most 1), against 0.9–1.2 % of nonzero entries with
  $x^m_c > 0$ under Eq. 44 in cells of realistic size. The displacement from $x^*_c$ becomes $\log(n_c/(n f^*_c))$
  instead of $n_c/(n f^*_c) - 1$: 4.6 instead of 96 for the cell of Figure 1 of the linearisation note.
- Calibration against the truth: for $n_c \ge 5$ the Gaussian (E) is a good approximation of the log-gamma
  likelihood and $|z| < 1, 2, 3$ should come out near 0.68, 0.95, 0.997; for $n_c = 1$ and $2$ the likelihood is
  skewed (mean $\psi(n_c) - \log N_c$, variance $\psi_1(n_c)$: $-0.58$ and $1.64$ at $n_c = 1$ against 0 and 1 for
  (E)), so under-coverage of order 10 % at $n_c = 1$ would not be a defect of the implementation.
- Zero-count entries are identical to the `posterior-mode` run; their calibration is unchanged.
- Combined with Sanity's own prior $N(\mu_g, v_g)$, (E) reproduces $x^*_c$ to 0.05 (median) and 0.3–0.5 (90 %)
  for single counts, 2.6 at worst (Section 6.1 of the linearisation note); Eq. 44 reproduces it exactly. This is the
  price of a prior-free measurement and is the quantity to watch if Bonsai's results move.
