# What Sanity writes, and the formula behind every number

Branch `new-direct-ltqs`. Equation numbers refer to Section 2 of E. van Nimwegen, "Getting proper error-bars in
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

Every numerical file is one row per gene. Text files have six decimals; matrices carry the gene name in the first
column and the cell names in a header line; the per-gene files hold one number per line. `-gz` gzips the text files,
`-npy` writes float64 `.npy` instead (matrices `(genes, cells)`, per-gene files `(genes, 1)`, `likelihood`
`(genes, bins)`, plus `likelihood_grid.npy` `(bins, 1)`); `geneID.txt` and `cellID.txt` are written whenever the
files do not carry names themselves (`-npy` or `-e 1`).

| file | value | error bar is |
|:--|:--|:--|
| `log_transcription_quotients` | $x^*_c$ (MARG: $x_c$) | |
| `ltq_error_bars` | $\epsilon_c$ (Eq. 41 / 48) | standard deviation; includes the $\mu$ uncertainty (the $1/C$ term) |
| `measured_ltq` | $x^m_c$ (Eq. 44; MARG: averaged) | |
| `measured_ltq_error_bars` | $1/\sqrt{n f^*_c}$ (Eq. 46; MARG: averaged as Eq. 48) | standard deviation |
| `mu` | $\mu$ (Eq. 33 / 49) | |
| `d_mu` | $\sigma_\mu$ (Eq. 43 / 50) | standard deviation |
| `variance` | $v$ | |
| `delta` | $x^*_c - \mu$ | |
| `d_delta` | $\epsilon_c$, identical to `ltq_error_bars` | kept for `Sanity_distance` / `Sanity_gene_correlation` |
| `likelihood` | $P_k$, with the grid $v_k$ as header line | |

### Which run writes which files

`--ltq-estimate` (`-ltq`) chooses the *default* set; `-e 1` writes everything whatever the flag says.

| mode | files |
|:--|:--|
| default (`-ltq posterior`) | `log_transcription_quotients`, `ltq_error_bars` |
| `-ltq measured` | `measured_ltq`, `measured_ltq_error_bars` |
| `-ltq both` | the four files above |
| `-e 1` (any `-ltq`) | the four files above + `mu`, `d_mu`, `variance`, `delta`, `d_delta`, `likelihood`, `geneID.txt`, `cellID.txt` |

`sanity_command.txt` records `# Method:` and `# LTQ estimate:`.

### Bonsai output mode

`--bonsai-output-mode` (requires MAP; ignores `-e`, `-gz`, `-npy`) writes the fixed set `delta.npy`, `d_delta.npy`
(float32, `(genes, cells)`), `mu.npy`, `variance.npy` (float64, `(genes, 1)`), `geneID.txt`, `cellID.txt`,
`sanity_command.txt`. The file names are the same as before, the content depends on `--ltq-estimate`:

| `--ltq-estimate` | `delta.npy` | `d_delta.npy` |
|:--|:--|:--|
| `measured` (default in this mode) | $x^m_c - \mu$ | $1/\sqrt{n f^*_c}$ |
| `posterior` | $x^*_c - \mu$ | $\epsilon_c$ |

`mu.npy` = $\mu$ and `variance.npy` = $v$ in both cases. `sanity_command.txt` says which it is
(`# LTQ estimate:` and the `# Output:` line).

## 3 Redundancy, for the discussion on a minimal output

- `d_delta` = `ltq_error_bars`, and `delta` = `log_transcription_quotients` − `mu`: both exist only for the two
  downstream tools, which read `delta.txt`, `d_delta.txt` and `variance.txt`.
- With `-e 1`, `ltq_error_bars` and `mu`/`d_mu` are consistent by construction: $\epsilon_c^2 = v_c + v_c^2 / (C(v - \langle v\rangle))$
  and $\sigma_\mu^2$ is the mean of the full covariance (Eq. 40), so $\sigma_\mu$ can be recomputed from $v$ and the $v_c$.
- `measured_ltq` and its error bar are functions of $x^*_c$, $\mu$, $v$ and $n f^*_c = N_c e^{x^*_c}$ (Eqs. 44, 46),
  so at a fixed $v$ they are derivable from the posterior files plus the count totals (see `docs/measured_ltq_linearisation.md`).
  For MARG they are not.
- `likelihood` is only needed to reconstruct the posterior over $v$.
