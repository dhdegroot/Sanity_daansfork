# The old and new L(v) differ by a constant

Notation of Section 2 of E. van Nimwegen, "Getting proper error-bars in Sanity" (the note), for one gene:
$C$ cells, counts $n_c$, $n = \sum_c n_c$, cell totals $N_c$, fractions $f_c$ with $\sum_c f_c = 1$, and
$q$ defined by $-q = \langle \log f \rangle - \langle \log N \rangle$. The old code (2021 SI Eqs. 20, 32, 33) uses
$\delta_c = \log f_c - \log N_c + q$ and computes

$$
L_{\rm old}(v) = -\frac{C}{2}\log v + \sum_c \left[n_c \delta_c - \frac{\delta_c^2}{2v}\right] - nq
- \frac12 \left[\sum_c \log\left(f_c + \tfrac{1}{nv}\right) + \log\left(1 - \sum_c \frac{f_c^2}{f_c + \frac{1}{nv}}\right)\right].
$$

The note's Eqs. 36–38 give, with $x^*_c = \log f_c + \log n - \log N_c$,

$$
L_{\rm new}(v) = -\frac{C-1}{2}\log v + \sum_c \left[n_c x^*_c - n f_c - \frac{(x^*_c - \langle x^*\rangle)^2}{2v}\right]
- \frac12 \left[-C\log v + \sum_c \log(n f_c v + 1) + \log \left\langle \frac{n f_c v}{n f_c v + 1} \right\rangle\right].
$$

**Step 1, the optimum.** $x^*_c = \delta_c + \log n - q$, and $\langle x^* \rangle = \log n - q$ because
$\langle \delta \rangle = \langle \log f\rangle - \langle \log N \rangle + q = 0$. Hence $x^*_c - \langle x^*\rangle = \delta_c$,
$\sum_c n_c x^*_c = \sum_c n_c \delta_c + n(\log n - q)$ and $\sum_c n f_c = n$, so

$$
\sum_c \left[n_c x^*_c - n f_c - \frac{(x^*_c - \langle x^*\rangle)^2}{2v}\right]
= \sum_c \left[n_c \delta_c - \frac{\delta_c^2}{2v}\right] - nq + n\log n - n .
$$

**Step 2, the sum of logarithms.** $n f_c v + 1 = nv\,(f_c + \tfrac{1}{nv})$, so
$\sum_c \log(n f_c v + 1) = C \log(nv) + \sum_c \log(f_c + \tfrac{1}{nv})$.

**Step 3, the last term.** Using $\sum_c f_c = 1$,

$$
1 - \sum_c \frac{f_c^2}{f_c + \frac{1}{nv}} = \sum_c f_c \left(1 - \frac{f_c}{f_c + \frac{1}{nv}}\right)
= \frac{1}{nv} \sum_c \frac{f_c}{f_c + \frac{1}{nv}} = \frac{1}{nv} \sum_c \frac{n f_c v}{n f_c v + 1}
= \frac{C}{nv} \left\langle \frac{n f_c v}{n f_c v + 1} \right\rangle,
$$

so $\log\langle \tfrac{n f_c v}{n f_c v + 1}\rangle = \log\left(1 - \sum_c \tfrac{f_c^2}{f_c + 1/(nv)}\right) + \log(nv) - \log C$.

**Step 4, the determinant.** Steps 2 and 3 give

$$
\log\det M = -C\log v + C\log(nv) + \sum_c \log\left(f_c + \tfrac{1}{nv}\right) + \log\left(1 - \sum_c \tfrac{f_c^2}{f_c + 1/(nv)}\right) + \log(nv) - \log C,
$$

and $-C\log v + C\log(nv) = C\log n$, so $\log\det M = \log\det(M/n)_{\rm old} + C\log n + \log v + \log n - \log C$.
The $\log(nv)$ that Step 3 produces is the $\log v$ that the determinant carries beyond the old form.

**Step 5, collect.** Substituting Steps 1 and 4 into $L_{\rm new}$:

$$
L_{\rm new}(v) = -\frac{C-1}{2}\log v - \frac12 \log v + \sum_c \left[n_c \delta_c - \frac{\delta_c^2}{2v}\right] - nq
- \frac12 \log\det(M/n)_{\rm old} + \left[n\log n - n - \frac{C+1}{2}\log n + \frac12 \log C\right].
$$

The first two terms are $-\frac{C}{2}\log v$, so

$$
L_{\rm new}(v) = L_{\rm old}(v) + n\log n - n - \frac{C+1}{2}\log n + \frac12\log C ,
$$

a constant in $v$. The extra $\frac12 \log v$ from $C \to C-1$ in Eq. 36 is cancelled by the $\log v$ hidden in the
last term of Eq. 38. Therefore the posterior over the variance grid, and the v selected by MLE, MAP and EAP, are the
same as before; the inferred variances are not "a little higher" as anticipated at the end of Section 2.5.

**Numerical check.** On `tests/count_table.tsv` (300 genes, 160 grid points) the binary with Eqs. 36–38 gives
`variance` identical to the `clean-up` build for MAP, MLE and EAP, and `likelihood` (the normalised posterior)
equal to a relative $1.2\cdot 10^{-10}$; the differences are rounding, as the two forms sum different terms.
