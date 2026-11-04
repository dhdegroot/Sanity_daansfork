# The measurement LTQs of Eq. 44 are a quadratic expansion around the posterior mode, and it breaks when the observed count is far from the expected count

Notation of Section 2 of "Getting proper error-bars in Sanity", one gene, index $g$ dropped: counts $n_c$,
$n = \sum_c n_c$, totals $N_c$, prior variance $v$, optimal fractions $f^*_c$ with $n f^*_c = N_c e^{x^*_c}$ the
expected count of cell $c$ at the posterior mode $x^*_c$ (Eq. 32).

## 1 What Eq. 44 is

Eq. 44 divides the prior out of the Laplace posterior. The log posterior at fixed $v$ is
$L(x|v) = \sum_c \ell_c(x_c) - \sum_c (x_c - \langle x \rangle)^2 / (2v)$ with the Poisson log-likelihood of one cell

$$
\ell_c(x) = n_c x - N_c e^{x}.
$$

At the mode $x^*$ the gradient vanishes, which is Eq. 28/30:

$$
\ell'_c(x^*_c) = n_c - N_c e^{x^*_c} = n_c - n f^*_c = \frac{x^*_c - \langle x^*\rangle}{v}. \tag{A}
$$

The Laplace approximation replaces $\ell_c$ by its second-order Taylor expansion around $x^*_c$, with
$\ell''_c(x^*_c) = -n f^*_c$:

$$
\ell_c(x) \approx \ell_c(x^*_c) + (n_c - n f^*_c)(x - x^*_c) - \tfrac12\, n f^*_c\,(x - x^*_c)^2. \tag{B}
$$

A quadratic in $x$ is a Gaussian "measurement": (B) equals $-\tfrac12\, n f^*_c\,(x - x^m_c)^2 + \text{const}$ with

$$
x^m_c = x^*_c + \frac{n_c - n f^*_c}{n f^*_c} = x^*_c + \frac{x^*_c - \langle x^*\rangle}{n f^*_c\, v}, \qquad
\langle (x_c - x^m_c)^2 \rangle = \frac{1}{n f^*_c},
$$

where the second form uses (A). These are exactly Eqs. 44 and 46. So, at the optimum,

$$
\boxed{\,x^m_c - x^*_c = \frac{n_c}{n f^*_c} - 1\,} \tag{C}
$$

the measurement LTQ lies above the posterior mode by the ratio of the observed count to the expected count, minus one.
(C) is exact, not a numerical artefact: the code computes Eq. 44 as written and the identity holds to $10^{-9}$.

## 2 Where it breaks

The expansion (B) is accurate only while $|x - x^*_c| \lesssim 1$, but its own maximum sits at distance
$n_c / (n f^*_c) - 1$ from $x^*_c$. Whenever the observed count exceeds the expected count by a large factor, Eq. 44
places the measurement far outside the region where the expansion is valid. Compare with the true Poisson
log-likelihood:

| | quadratic expansion (B), i.e. Eqs. 44/46 | true $\ell_c(x) = n_c x - N_c e^x$ |
|:--|:--|:--|
| maximum at | $x^*_c + n_c/(n f^*_c) - 1$ | $\hat x_c = \log(n_c / N_c) = x^*_c + \log\!\big(n_c / (n f^*_c)\big)$ |
| curvature there | $n f^*_c$ | $n_c$ |
| $n_c = 0$ | maximum at $x^*_c - 1$, width $1/\sqrt{n f^*_c}$ | no maximum: $\ell_c$ decreases monotonically, flat for $x \ll -\log N_c$ |

For a ratio $n_c/(n f^*_c) = r$ the two maxima differ by $r - 1$ versus $\log r$: for $r = 10^4$ that is $10^4$ versus $9.2$.

When is $r$ large? From (A), $x^*_c - \langle x^*\rangle = v\,(n_c - n f^*_c)$. If the prior is tight compared with the
sampling noise, $n f^*_c v \ll 1$, the mode hardly moves away from the gene mean, $x^*_c \approx \langle x^* \rangle + v n_c$,
so $n f^*_c \approx N_c e^{\langle x^*\rangle}$ is the count expected *at the gene mean*, and $r$ is the observed count
over that. A single count in a cell where the gene mean predicts $10^{-4}$ counts gives $r = 10^4$. This is the regime
of weakly expressed genes, small cells and small $v$, i.e. the bulk of a UMI matrix.

## 3 Numbers

`tests/count_table.tsv` (283 genes, 498 cells), MAP:

| | |
|:--|:--|
| entries with $\lvert x^m_c - x^*_c\rvert > 5$ | 1870 of 140 934 (1.3 %), in 247 genes, at all $v$ from $10^{-3}$ to $5$ |
| largest | gene `ENSG00000006047.12`: $n = 4$, $v = 0.001$, cell with $n_c = 1$ and $N_c = 3$ |
| | $x^*_c = -10.478$, $\langle x^* \rangle = -10.479$, $n f^*_c = 8.4 \cdot 10^{-5}$ |
| | $x^m_c = 11\,825$ with error bar $1/\sqrt{n f^*_c} = 109$; the true $\hat x_c = \log(1/3) = -1.1$ |

Synthetic matrix of 3000 genes × 3000 cells (Poisson counts from log-normal LTQs): 1.2 % of entries with
$\lvert x^m_c - x^*_c\rvert > 5$, largest $1116$.

For every zero-count cell, (C) gives $x^m_c = x^*_c - 1$ exactly, with the error bar $1/\sqrt{n f^*_c}$ of Section 2.8.

## 4 Why it matters downstream

A Gaussian measurement $x^m_c \pm 1/\sqrt{n f^*_c}$ handed to Bonsai or MARA is, by (B), the quadratic expansion of the
Poisson log-likelihood around *Sanity's* posterior mode. Combined with Sanity's own prior $N(\langle x^*\rangle, v)$ it
reproduces $x^*_c$ exactly, by construction. Combined with anything else it is only as good as (B) at the place where
the downstream posterior ends up:

- With a wider prior variance $v'$ the downstream mode moves by $\approx v'(n_c - n f^*_c)$ from $x^*_c$ (the linear
  term of (B) dominates while $n f^*_c (x - x^*_c) \ll n_c$), with nothing to stop it before $x^m_c$. The true likelihood
  stops at $\hat x_c = \log(n_c/N_c)$.
- Treated as data with weight $n f^*_c$, the example cell says "$x = 11\,825 \pm 109$": a 100-sigma statement that
  the true likelihood never makes. Any downstream step that averages measurements with their inverse variances
  (over cells of one type, over genes of a motif) sees these entries as extreme, confidently measured values.
- The problem is largest exactly where Section 2.8 expected trouble only for zero counts: it is the nonzero counts
  of weakly expressed genes that explode, while zero counts sit quietly at $x^*_c - 1$ with a huge error bar.

## 5 Possible ways out (not implemented, for discussion)

1. **Expand around the measurement's own optimum.** For $n_c > 0$ the Poisson log-likelihood has its maximum at
   $\hat x_c = \log(n_c / N_c)$ with curvature $n_c$, giving $x^m_c = \log(n_c/N_c)$, $(\epsilon^m_c)^2 = 1/n_c$. This
   does not depend on Sanity at all; the Sanity-specific information is then entirely in $\mu$, $v$, $x^*_c$ and
   $\epsilon_c$. Written relative to the mode, $x^m_c = x^*_c + \log\big(n_c / (n f^*_c)\big)$, i.e. (C) with $r - 1$
   replaced by $\log r$. Zero counts have no Gaussian representation; they are a one-sided constraint
   $x \lesssim -\log N_c$ (the likelihood is $e^{-N_c e^x}$).
2. **Keep Eq. 44 but recognise it as a Taylor expansion**, and pass on the two coefficients of (B), the gradient
   $n_c - n f^*_c$ and the curvature $n f^*_c$, instead of a mean and a variance; a downstream Gaussian update is
   the same, but the consumer can bound the step.
3. **Divide the prior out in the downstream model**, i.e. hand Bonsai $x^*_c$, $\epsilon_c$, $\mu$ and $v$ and let it
   subtract the Sanity prior from the posterior precision there, where the full Poisson likelihood (or a better
   approximation) can be used.

The code on `new-direct-ltqs` implements Eqs. 44 and 46 as written (`measured_ltq`, `measured_ltq_error_bars`, and
`delta.npy`/`d_delta.npy` in Bonsai output mode), so that the behaviour above can be measured on the simulated data.
