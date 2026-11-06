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

`tests/count_table.tsv` (283 genes, 498 cells, cell sizes $N_c$ from 1 to 3590 UMIs with median 47) and a synthetic
matrix of 3000 genes × 3000 cells (Poisson counts from log-normal LTQs, $N_c$ from 2405 to $3.4 \cdot 10^6$ with median
20 546), MAP: 1.3 % and 1.2 % of all entries have $\lvert x^m_c - x^*_c \rvert > 5$. The most extreme entries sit in
cells with a handful of UMIs (the largest, $x^m_c = 11\,825$, is a single count in a cell with 3 UMIs), but the effect
does not need small cells. What matters is how rare the gene is in the cell, $\lambda_c = N_c e^\mu$, and a gene seen a
handful of times in thousands of cells has $\lambda_c \sim 10^{-3}$ to $10^{-2}$ in a cell of 10 000 UMIs too.
Nonzero entries by cell size:

| $N_c$ | data | nonzero entries | $x^m - x^* > 1$ | $> 5$ | $> 100$ | largest |
|:--|:--|--:|--:|--:|--:|--:|
| 100–1000 | real | 8 153 | 18.5 % | 5.6 % | 0.23 % | 1124 |
| 1000–3000 | real | 5 892 | 9.3 % | 2.3 % | 0.02 % | 133 |
| 3000–10 000 | synthetic | 377 049 | 13.1 % | 3.5 % | 0.26 % | 1106 |
| 10 000–30 000 | synthetic | 2 482 586 | 10.3 % | 2.8 % | 0.13 % | 1116 |
| $> 30\,000$ | synthetic | 1 386 275 | 8.3 % | 2.0 % | 0.02 % | 307 |

Typical bad entries in cells of realistic size (MAP; $\lambda_c$ is the count expected at the gene mean):

| | real: `ENSG00000004799.7` | synthetic: `gene_214` | synthetic: `gene_2352` | synthetic: `gene_1117` |
|:--|--:|--:|--:|--:|
| gene: total count $n$, zeros | 8, 98.8 % | 7, 99.8 % | 30, 99 % | 175, 95 % |
| selected $v$ | 0.001 | 0.001 | 0.001 | 0.69 |
| cell size $N_c$, count $n_c$ | 1159, 3 | 10 460, 1 | 9216, 2 | 2405, 1 |
| $\lambda_c$ | 0.065 | 0.0009 | 0.0034 | 0.0052 |
| posterior mode $x^*_c$ | $-9.8$ | $-16.3$ | $-14.8$ | $-12.4$ |
| Eq. 44/46: $x^m_c \pm \epsilon^m_c$ | $+35.1 \pm 3.9$ | $+1100 \pm 33$ | $+575 \pm 17$ | $+83 \pm 10$ |
| likelihood peak $\log(n_c / N_c)$ | $-6.0$ | $-9.3$ | $-8.4$ | $-7.8$ |
| $(x^m_c - x^*_c)/\epsilon^m_c$ | 11.5 | 33 | 34 | 9.7 |

All four measurement LTQs are positive, i.e. they claim that the gene makes up more than all of the cell's mRNA, with
a small error bar. This is not rare: in the real cells with at least 1000 UMIs 57 of 6301 nonzero entries (0.9 %, in
41 genes) have $x^m_c > 0$, in the synthetic matrix 49 555 of 4.2 million (1.2 %, in 1284 of 3000 genes). Nor is it
limited to $v$ at the grid floor: `gene_1117` has 175 counts and $v = 0.69$.

![Prior, Poisson likelihood, posterior and the Gaussian approximation of Eqs. 44/46 for the `gene_1117` case above](figures/measured_ltq_example.png)

**Figure 1.** The cell of `gene_1117` (synthetic) from the table above. Top panel: Sanity's prior
$N(\mu, v)$, the Poisson likelihood $e^{n_c x - N_c e^x}$ of the cell's count, the posterior (their product, with
mode $x^*_c$), and the Gaussian of Eqs. 44/46, which touches the likelihood at $x^*_c$ with the same value, slope and
curvature (expansion (B)). Bottom panel: the same on a log scale; the Gaussian peaks at $x^m_c$, in the region $x > 0$,
while the true likelihood has collapsed at $\log(n_c/N_c)$. The gene's count distribution is printed above the panels.

For every zero-count cell, (C) gives $x^m_c = x^*_c - 1$ exactly, with the error bar $1/\sqrt{n f^*_c}$ of Section 2.8.

## 4 Why it matters downstream

A Gaussian measurement $x^m_c \pm 1/\sqrt{n f^*_c}$ handed to Bonsai or MARA is, by (B), the quadratic expansion of the
Poisson log-likelihood around *Sanity's* posterior mode. Combined with Sanity's own prior $N(\langle x^*\rangle, v)$ it
reproduces $x^*_c$ exactly, by construction. Combined with anything else it is only as good as (B) at the place where
the downstream posterior ends up:

- With a wider prior variance $v'$ the downstream mode moves by $\approx v'(n_c - n f^*_c)$ from $x^*_c$ (the linear
  term of (B) dominates while $n f^*_c (x - x^*_c) \ll n_c$), with nothing to stop it before $x^m_c$. The true likelihood
  stops at $\hat x_c = \log(n_c/N_c)$.
- Treated as data with weight $n f^*_c$, the real example cell of Section 3 (3 counts in 1159 UMIs) says
  "$x = 35 \pm 3.9$": an 11-sigma statement that the gene is $e^{35}$ times the cell's whole mRNA, which the true
  likelihood never makes; the cell of Figure 1 says "$x = 83 \pm 9.8$". Any downstream step that averages measurements with their inverse variances
  (over cells of one type, over genes of a motif) sees these entries as extreme, confidently measured values.
- The problem is largest exactly where Section 2.8 expected trouble only for zero counts: it is the nonzero counts
  of weakly expressed genes that explode, while zero counts sit quietly at $x^*_c - 1$ with a huge error bar.

## 5 Which entries and genes are affected

By (C) the displacement is a property of the entry, not of the gene. Writing $\lambda_c = N_c e^{\mu}$ for the count
the gene mean predicts in cell $c$, the mode equation (A) gives two regimes for a nonzero count:

$$
x^m_c - x^*_c \approx
\begin{cases}
n_c / \lambda_c - 1 & \text{prior wins: } v n_c \ll \log(n_c/\lambda_c) \\[2pt]
\log(n_c/\lambda_c) / (v n_c) & \text{likelihood wins: } v n_c \gg \log(n_c/\lambda_c)
\end{cases}
\tag{D}
$$

So an entry is bad (displacement above 1, outside the validity of (B)) when the observed count exceeds the count
expected at the gene mean by a factor of about 2 or more, and the prior is not overwhelmed. The minimum of the two
lines of (D), which needs only the counts and Sanity's $\mu$ and $v$, classifies nonzero entries at the threshold
"displacement $> 5$" correctly in 96 % (real) and 99 % (synthetic) of the cases; replacing $\lambda_c$ by the
counts-only $n N_c / N_{tot}$ still gets 85 % and 97 %. Zero counts always sit at exactly $x^*_c - 1$.

Per gene, the driver is the mean count per cell $n/C$ (synthetic 3000 × 3000 matrix, MAP; "bad" = displacement $> 1$):

| $n/C$ | genes | median $n$ | median zeros | median $v$ | median $\max_c (x^m_c - x^*_c)$ | nonzero entries bad | genes with an entry $> 5$ |
|:--|--:|--:|--:|--:|--:|--:|--:|
| $< 0.03$ | 467 | 23–52 | 98–99 % | 0.001 | 210–390 | 99.9 % | 100 % |
| 0.03–0.1 | 425 | 168 | 95 % | 0.001 | 60 | 99 % | 100 % |
| 0.1–0.3 | 383 | 503 | 86 % | 0.49 | 19 | 86 % | 99 % |
| 0.3–1 | 359 | 1669 | 65 % | 0.40 | 7 | 37 % | 78 % |
| 1–3 | 334 | 5174 | 33 % | 0.35 | 3 | 6 % | 13 % |
| 3–10 | 375 | 15955 | 10 % | 0.37 | 1.3 | 0.3 % | 0 % |
| $> 10$ | 657 | $> 5 \cdot 10^4$ | $< 1$ % | 0.35–0.40 | $< 0.5$ | 0 % | 0 % |

- **Sparse genes** ($n/C \lesssim 0.3$) are bad in every nonzero entry: tens to a few hundred counts in total, more than
  85 % zeros, $v$ at the grid floor. There $x^*_c \approx \mu$, so $n f^*_c \approx \lambda_c \ll 1$ and every single count
  is a surprise by the factor $1/\lambda_c$. The z-score of the measurement, $(x^m_c - x^*_c)/\epsilon^m_c =
  (n_c - n f^*_c)/\sqrt{n f^*_c}$, reaches 34 in the synthetic data and 11.5 in the
  real cells with at least 1000 UMIs (109 in a real cell of 3 UMIs).
- **Well-expressed genes** ($n/C \gtrsim 3$) are fine in the synthetic data: no entry above 5.
- **In between** it is a mixture, and in real data the mixture is wider because cell sizes are uneven. In the 283 × 498
  test matrix ($N_c$ from 1 to 3590, median 47) genes with thousands of counts and $n/C$ up to 12 still have
  displacements of 5 to 11, each time in a cell with 3 to 7 UMIs. There the bad entries are the small cells, not the
  gene; 10x data has a smaller spread but still a factor 10 to 50, so the mixed band persists for $n/C$ between about
  0.3 and 3.
- Total count, number of zeros and $v$ all correlate with badness, but only as faces of the same variable, the expected
  count per cell. In the real data "$v$ at the grid floor" is nearly equivalent to "every nonzero entry bad" (99 % of
  such genes versus 1 % of the others).

**Flagging.** The per-gene fraction of nonzero entries with $n_c > 2 n f^*_c$ is bimodal (synthetic: 1318 genes below
10 %, 1092 above 90 %), so a cut at 50 % is robust. It flags 47 % of the synthetic genes and 34 % of the real ones,
holding 0.5 % and 1.8 % of all UMIs but about half of the detected transcriptome, and leaves the mixed band (187 of
283 real genes with 10 to 90 % bad entries) unflagged. The natural flag is per entry: displacement above 1, available
today from the standard output as `measured_ltq` minus `log_transcription_quotients`, or from the counts as
$n_c / (N_c e^{x^*_c}) - 1$. The entries a gene flag would discard, the nonzero counts of lowly expressed genes, are
the only informative entries those genes have. The posterior output $x^*_c \pm \epsilon_c$ is unaffected throughout.

## 6 Possible ways out (not implemented, for discussion)

### 6.1 Expand around the measurement's own optimum (nonzero counts)

For $n_c > 0$ the Poisson log-likelihood $\ell_c(x) = n_c x - N_c e^x$ has its maximum at
$\hat x_c = \log(n_c/N_c)$ with curvature $\ell''_c(\hat x_c) = -n_c$, giving

$$
x^m_c = \log(n_c / N_c), \qquad (\epsilon^m_c)^2 = 1/n_c . \tag{E}
$$

This does not depend on Sanity at all; the Sanity-specific information is then entirely in $\mu$, $v$, $x^*_c$ and
$\epsilon_c$. Relative to the mode, $x^m_c = x^*_c + \log\big(n_c/(n f^*_c)\big)$, i.e. (C) with $r - 1$ replaced by
$\log r$: the displacement of the cell of Figure 1 drops from 96 to 4.6 (for the real cell of Section 3 from 45 to 3.8), and $x^m_c$
becomes $\log(n_c/N_c)$, the peak of the orange curve.

(E) is the Laplace approximation of a likelihood that is not Gaussian in $x$: under a flat prior in $x$ the posterior
of $\lambda = N_c e^x$ is $\mathrm{Gamma}(n_c, 1)$, so $x$ is log-gamma with mean $\psi(n_c) - \log N_c$ and variance
$\psi_1(n_c)$ (the SI 2021 formulas), skewed with an exponential left tail: for $x \ll \hat x_c$ the likelihood decays as
$e^{n_c (x - \hat x_c)}$, not as a Gaussian. The error this causes is bounded. Combining the Gaussian (E) with Sanity's
own prior $N(\mu, v)$ gives the mode $\tilde x_c = (n_c \hat x_c + \mu/v)/(n_c + 1/v)$; its distance from the exact
mode $x^*_c$ is (MAP, both matrices; the moment-matched Gaussian $\psi(n_c) - \log N_c \pm \sqrt{\psi_1(n_c)}$ is not
better and is listed for comparison):

| $n_c$ | cells (real / synthetic) | $\lvert\tilde x_c - x^*_c\rvert$ median | 90 % | max | moment-matched median | 90 % | max |
|:--|--:|--:|--:|--:|--:|--:|--:|
| 1 | 8 931 / 932 060 | 0.05 / 0.04 | 0.47 / 0.26 | 2.2 / 2.6 | 0.18 / 0.10 | 0.42 / 0.24 | 1.2 / 3.0 |
| 2 | 3 034 / 461 641 | 0.02 / 0.02 | 0.21 / 0.16 | 1.1 / 2.6 | 0.16 / 0.08 | 0.24 / 0.17 | 0.6 / 2.8 |
| 5–9 | 2 203 / 695 361 | 0.005 / 0.005 | 0.03 / 0.04 | 0.4 / 1.4 | 0.06 / 0.04 | 0.09 / 0.07 | 0.4 / 1.4 |
| $\geq 30$ | 666 / 644 044 | $< 10^{-3}$ | $< 10^{-2}$ | 0.01 / 0.5 | 0.01 / 0.01 | 0.02 / 0.01 | 0.02 / 0.5 |

The large values sit where (D) says: single counts with $\lambda_c \in [0.01, 0.1]$ and $v \approx 1$ (median error
0.34 to 0.39, 90 % at 0.6 to 0.8). For $v$ at the grid floor the error is below $0.005$ everywhere, because the prior
decides and the shape of the likelihood does not matter. Compare with Eq. 44, which is exact with Sanity's prior and
wrong by $r - 1$ with any other.

### 6.2 Zero counts

For $n_c = 0$ the likelihood is $e^{-N_c e^x}$: it has no maximum and, under a flat prior in $x$, no finite moments.
It is a soft one-sided bound,

$$
\ell_c(x) = -N_c e^x =
\begin{cases} \approx 0 & x \ll -\log N_c \\ \to -\infty & x \gg -\log N_c \end{cases}
$$

with the wall at $x = -\log N_c$, the LTQ at which one count would have been expected. No Gaussian represents it, so
option 1 needs a separate rule for zeros. The candidates:

1. **Drop the zero** ($\epsilon^m_c = \infty$). Prior-free, but it loses the bound. Whether the bound carries
   information depends on $\lambda_c = N_c e^\mu$, the count the gene mean predicts in that cell. The data say how
   much (MAP; $x^*_c - \mu$ is how far Sanity's posterior mode moves below the gene mean because of the zero):

   | $\lambda_c$ | zero cells (real / synthetic) | median $x^*_c - \mu$ | min | Eq. 44 bar $1/\sqrt{n f^*_c}$, median | min |
   |:--|--:|--:|--:|--:|--:|
   | $< 0.01$ | 44 942 / 669 001 | 0.00 / 0.00 | $-0.16$ / $-0.05$ | 19 / 13 | 10 / 10 |
   | 0.01–0.1 | 48 498 / 2 189 611 | $-0.03$ / 0.00 | $-0.68$ / $-0.49$ | 6.0 / 5.5 | 3.2 / 3.2 |
   | 0.1–0.5 | 18 846 / 1 253 810 | $-0.17$ / $-0.09$ | $-1.35$ / $-1.20$ | 2.5 / 2.4 | 1.4 / 1.4 |
   | 0.5–2 | 7 321 / 519 472 | $-0.61$ / $-0.28$ | $-2.2$ / $-2.3$ | 1.5 / 1.3 | 0.7 / 0.7 |
   | 2–10 | 1 013 / 113 020 | $-1.17$ / $-0.77$ | $-2.2$ / $-3.5$ | 1.1 / 0.8 | 0.5 / 0.3 |
   | $> 10$ | 2 / 6 831 | $-2.0$ / $-2.1$ | $-2.2$ / $-5.8$ | 0.7 / 0.7 | 0.7 / 0.2 |

   For $\lambda_c < 0.1$ (77 % of the real zeros, 60 % of the synthetic ones) the zero is uninformative and dropping it
   costs nothing. For $\lambda_c > 0.5$ (7 % and 13 % of the zeros) it moves the cell by 0.3 to 6 in log space, and
   dropping it hands those cells to the downstream prior, which puts them at the level of their neighbours.
2. **Keep Eq. 44 for zeros.** Expanding $-N_c e^x$ around $x^*_c$ gives exactly
   $x^m_c = x^*_c - 1$, $(\epsilon^m_c)^2 = 1/(N_c e^{x^*_c})$. Unlike the nonzero case this never blows up: at the
   mode $v\, N_c e^{x^*_c} = \mu - x^*_c$, so the precision $N_c e^{x^*_c}$ is at most a few (bar $\geq 0.2$ in the
   table) and the centre is within 1 of the mode. It is prior-dependent, but in the right way: when the zero is
   uninformative the bar is wide (10 and more), when it is informative the Gaussian sits about 1 below the wall with a
   bar of order 1, a fair local stand-in for a one-sided constraint. Its defect is symmetry: it penalises
   $x \ll x^*_c - 1$, which the true likelihood does not. Alternatively the wall itself can be used,
   $x^m_c = -\log N_c - 1$ with bar 1 (Eq. 44 at $n f^*_c = 1$), which is prior-free and has the same shape; it is
   right when the zero is informative and too narrow when it is not (a cell with $N_c = 100$ then claims
   $x \approx -5.6 \pm 1$ for a gene whose mean is $-12$).
3. **Pass the sufficient statistics** $(n_c, N_c)$ and let the downstream use the Poisson likelihood itself (option 3
   below). This is the only treatment that is exact for zeros and nonzeros alike.

### 6.3 Are zeros and nonzeros consistent then?

Not in every sense, and it is worth being precise about which.

- **Per cell, with Sanity's own prior.** Eq. 44 reproduces $x^*_c$ exactly for every cell by construction. Under
  option 1 with Eq. 44 kept for zeros, zero cells still reproduce $x^*_c$ exactly, nonzero cells to within the errors
  of the table in 6.1: 0.05 median and 0.5 at the 90th percentile for single counts, 2.6 at worst. With any other
  prior Eq. 44 is wrong by $r - 1$ for nonzero counts, option 1 by at most the log-gamma skew. If zeros are dropped,
  zero cells land at the downstream prior's mean instead of $x^*_c$, which is harmless for $\lambda_c < 0.1$ and
  wrong by 0.3 to 6 above.
- **Per gene.** The nonzero measurements (E) alone do not know where the gene sits. For sparse genes their mean,
  $\langle \log(n_c/N_c) \rangle_{n_c > 0}$, lies far above $\mu$: by 4.1 to 4.3 (median) for $n/C < 0.03$, by 1.9 to
  2.4 for $n/C$ between 0.03 and 0.3, by 0.6 to 1.2 for $n/C$ between 0.3 and 3, and by nothing above $n/C = 3$. It is
  the zeros, not the nonzeros, that hold $\mu$ down. So a downstream model that takes (E) for the nonzero cells must
  take the gene level from somewhere: from Sanity's $\mu$ and $v$ as its prior, or from zero-count measurements that
  carry the bound. Dropping the zeros and also not using $\mu$ would overestimate every sparse gene by those amounts.
- **Scale.** Both kinds of measurement are absolute LTQs, so they live on the same axis; the mix of a prior-free rule
  for nonzeros and a prior-dependent one for zeros is inconsistent in origin, not in units. The prior enters the zero
  cells only where the zero is informative, i.e. where $\lambda_c \gtrsim 0.5$; there it moves the Gaussian by a
  bounded amount (at most the $x^*_c - \mu$ of the table, typically 0.3 to 1).

Summary: the Gaussian-per-cell format cannot be consistent for zeros and nonzeros in all prior regimes, because the
zero-count likelihood has no Gaussian form. The pragmatic combination is (E) for nonzero counts and Eq. 44 for zeros,
with $\mu$ and $v$ passed alongside; its worst-case error is of order 1 in log space, against $r - 1$ for Eq. 44 as
it stands. Exact consistency needs the counts themselves.

### 6.4 Other options

- **Keep Eq. 44 but recognise it as a Taylor expansion**, and pass on the two coefficients of (B), the gradient
  $n_c - n f^*_c$ and the curvature $n f^*_c$, instead of a mean and a variance; a downstream Gaussian update is the
  same, but the consumer can bound the step. This treats zeros and nonzeros identically but does not remove the
  problem, it only makes it visible to the consumer.
- **Divide the prior out in the downstream model**, i.e. hand Bonsai $x^*_c$, $\epsilon_c$, $\mu$ and $v$ (or the
  counts) and let it subtract the Sanity prior from the posterior precision there, where the full Poisson likelihood
  can be used. This is the exact version of option 1.

The code on `clean-up-meas` implements both: `--measured-expansion posterior-mode` (the default) is Eqs. 44 and 46
as written, `--measured-expansion likelihood-maximum` is option 1 with Eq. 44/46 kept for zero counts (Sections 6.1
and 6.2). Both write `measured_ltq`, `measured_ltq_error_bars`, and `delta.npy`/`d_delta.npy` in Bonsai output mode,
so that the behaviour above can be measured on the simulated data (`docs/hpc_comparison_prompt.md`).
All numbers in Sections 3, 5 and 6 come from MAP runs on `tests/count_table.tsv` (283 × 498) and on a synthetic
3000 × 3000 matrix (Poisson counts from log-normal LTQs, $n/C$ spread evenly in log); the scripts are not part of the
repository.
