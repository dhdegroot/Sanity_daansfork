# What Sanity writes, and the formula behind every number

Branch `clean-up-meas`. Equation numbers refer to Section 2 of E. van Nimwegen, "Getting proper error-bars in
Sanity" (the note), unless marked SI (the 2021 Supplementary Information). Per gene $g$ (index dropped below):
counts $n_c$, $n = \sum_c n_c$, cell totals $N_c$, $C$ cells, variance grid $v_k$ ($k = 0,\ldots,K-1$).

## 1 What is computed per gene

**Stage 1, every grid point $v_k$** (kept per grid point: four numbers, no per-cell arrays):

| kept | formula |
|:--|:--|
| $q_k$ | fixed by $\sum_c f^*_c = 1$ with $f^*_c = \frac{1}{nv} W\!\left(nv\, e^{v n_c + \log N_c - q}\right)$ (Eq. 31) |
| $L_k = L(v_k)$ | Eqs. 36–38 (log posterior at the optimum minus $\tfrac12 \log\det M$), see `docs/likelihood_equivalence.md` |
| $\mu(v_k) = \langle x^* \rangle$ | mean over cells of $x^*_c = \log f^*_c + \log n - \log N_c$ (Eqs. 32, 33) |
| $\sigma_\mu^2(v_k)$ | $\dfrac{\langle v \rangle}{C}\,\dfrac{1}{1 - \langle v\rangle / v}$ with $v_c = \dfrac{v}{n f^*_c v + 1}$ (Eqs. 39, 43) |

Then $P_k = e^{L_k} / \sum_{k'} e^{L_{k'}}$ (SI Eq. 34), the posterior over the grid under the $1/v$ prior implied by the
log-spaced grid.

**Stage 2a, MLE / MAP / EAP:** one grid point $k^*$ is chosen (MLE: $\arg\max L_k$; MAP: $\arg\max P_k / v_k$;
EAP: $v_k$ closest to $\sum_k P_k v_k$). At $v = v_{k^*}$, $f^*_c$ is rebuilt from $q_{k^*}$ and two passes over the
cells give

| quantity | formula | passes |
|:--|:--|:--|
| $x^*_c$ | $\log f^*_c + \log n - \log N_c$ (Eq. 32) | 1 |
| $\mu$ | $\langle x^* \rangle$ (Eq. 33) | 1 |
| $\langle v \rangle$, $v - \langle v \rangle$ | $\tfrac{v}{C}\sum_c \tfrac{1}{n f^*_c v + 1}$ and $\tfrac{v}{C}\sum_c \tfrac{n f^*_c v}{n f^*_c v + 1}$ (no cancellation) | 1 |
| $\epsilon_c^2$ | $v_c + \dfrac{1}{C}\dfrac{v_c^2}{v - \langle v \rangle}$ (Eq. 41) | 2 |
| $\sigma_\mu^2$ | Eq. 43 | after pass 1 |
| $x^m_c$ | $x^*_c + \dfrac{x^*_c - \langle x^* \rangle}{n f^*_c v}$ (Eq. 44) | 2 |
| $(\epsilon^m_c)^2$ | $1 / (n f^*_c)$ (Eq. 46) | 2 |
| $v$ | $v_{k^*}$ | |

The renormalised variances $v_c$ are computed in both passes and never stored (Section 2.10).

**Stage 2b, MARG:** the same per-cell quantities are computed at every grid point with $P_k > 0$ and averaged
bin by bin with West's weighted update (only the $x^*_c(v_k)$ of the current bin are held):

| quantity | formula |
|:--|:--|
| $x_c$ | $\sum_k P_k\, x^*_c(v_k)$ (Eq. 47) |
| $\epsilon_c^2$ | $\sum_k P_k \left[\epsilon_c^2(v_k) + (x^*_c(v_k) - x_c)^2\right]$ (Eq. 48) |
| $x^m_c$, $(\epsilon^m_c)^2$ | the same two forms applied to $x^m_c(v_k)$ and $1/(n f^*_c(v_k))$; **not in the note** |
| $\mu$ | $\sum_k P_k\, \mu(v_k)$ (Eq. 49) |
| $\sigma_\mu^2$ | $\sum_k P_k \left[\sigma_\mu^2(v_k) + (\mu(v_k) - \mu)^2\right]$ (Eq. 50) |
| $v$ | $\sum_k P_k v_k$ |

**Not computed any more:** the asymmetric zero-count error bar (SI Eq. 38), and $\mu = \psi(n) - q$,
$\sigma_\mu^2 = \psi_1(n)$ (SI Eqs. 44–47). $\mu$ therefore shifts by $\log n - \psi(n)$ relative to `clean-up`
($0.58$ at $n = 1$, $0.05$ at $n = 10$, $0.005$ at $n = 100$). $\delta_c = x^*_c - \mu$ is unchanged.

## 2 Files

Every numerical file is one row per gene, and every file name says what the file holds, whatever the run's
options. Text files have six decimals; matrices carry the gene name in the first column and the cell names in a
header line; the per-gene files hold one number per line. `-gz` gzips the text files, `-npy` writes `.npy` instead:
matrices `(genes, cells)` and `gene_variances_likelihood` `(genes, bins)` as float32, the per-gene files `(genes, 1)`
and `gene_variances_grid.npy` `(bins, 1)` as float64. `geneID.txt`, `cellID.txt` and `sanity_command.txt` are
always written.

| file | value | error bar is |
|:--|:--|:--|
| `ltqs_posterior` | $x^*_c$ (MARG: $x_c$) | |
| `ltqs_posterior_error_bars` | $\epsilon_c$ (Eq. 41 / 48) | standard deviation; includes the $\mu$ uncertainty (the $1/C$ term) |
| `ltqs_measurement_approx` | $x^m_c$ (Eq. 44; MARG: averaged) | |
| `ltqs_measurement_approx_error_bars` | $1/\sqrt{n f^*_c}$ (Eq. 46; MARG: averaged as Eq. 48) | standard deviation |
| `gene_means` | $\mu$ (Eq. 33 / 49) | |
| `gene_means_error_bars` | $\sigma_\mu$ (Eq. 43 / 50) | standard deviation |
| `gene_variances` | $v$ | |
| `gene_variances_likelihood` | $P_k$, with the grid $v_k$ as header line (`gene_variances_grid.npy` with `-npy`) | |

### Which run writes which files

`--ltq-estimate` (`-ltq`) chooses the per-cell matrices; `-e 1` writes everything whatever the flag says; the
per-gene vectors are always written.

| mode | files |
|:--|:--|
| default (`-ltq posterior`) | `ltqs_posterior`, `ltqs_posterior_error_bars` |
| `-ltq measured` | `ltqs_measurement_approx`, `ltqs_measurement_approx_error_bars` |
| `-ltq both` | the four files above |
| `-e 1` (any `-ltq`) | the four files above + `gene_variances_likelihood` |
| always | `gene_means`, `gene_means_error_bars`, `gene_variances`, `geneID.txt`, `cellID.txt`, `sanity_command.txt` |

`sanity_command.txt` records `# Sanity version:`, `# Method:`, `# LTQ estimate:` and, in Bonsai output mode,
`# Output mode: bonsai`.

### Bonsai output mode

`--bonsai-output-mode` is the preset `-npy -v_m EAP --ltq-estimate measured`, nothing more: Bonsai reads
`ltqs_measurement_approx.npy`, `ltqs_measurement_approx_error_bars.npy`, `gene_means.npy` and `gene_variances.npy`
by name, and checks `# Output mode: bonsai` and `# LTQ estimate: measured`. An explicit `-v_m` other than `EAP` is
an error. Bonsai does not divide the prior out of these files, since Sanity has already done so exactly (Eqs. 44,
46); dividing it out of the posterior files cell by cell would be wrong, because $\epsilon_c$ carries the shared
uncertainty of $\mu$ (`docs/measured_ltq_linearisation.md`, Section 7).

## 3 Redundancy, for the discussion on a minimal output

- The log fold changes $\delta_c = x^*_c - \mu$ and their bars of versions before 3.0 are `ltqs_posterior` − `gene_means`
  and `ltqs_posterior_error_bars`; `Sanity_distance` and `Sanity_gene_correlation` compute them on reading.
- `ltqs_posterior_error_bars` and `gene_means`/`gene_means_error_bars` are consistent by construction: $\epsilon_c^2 = v_c + v_c^2 / (C(v - \langle v\rangle))$
  and $\sigma_\mu^2$ is the mean of the full covariance (Eq. 40), so $\sigma_\mu$ can be recomputed from $v$ and the $v_c$.
  The posterior variance with $\mu$ held fixed, $v_c = v / (n f^*_c v + 1)$ (Eq. 39), is not written; it is
  $1 / v_c = 1 / (\epsilon^m_c)^2 + 1 / v$ from the measurement error bar and the gene variance.
- `ltqs_measurement_approx` and its error bar are functions of $x^*_c$, $\mu$, $v$ and $n f^*_c = N_c e^{x^*_c}$ (Eqs. 44, 46),
  so at a fixed $v$ they are derivable from the posterior files plus the count totals (see `docs/measured_ltq_linearisation.md`).
  For MARG they are not.
- `gene_variances_likelihood` is only needed to reconstruct the posterior over $v$.
