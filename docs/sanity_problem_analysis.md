Hi all,

I was implementing these new equations in Sanity for getting the likelihood-estimates directly, and on the first small test dataset some things seemed to go wrong. The inferred mean and error-bar for the measured LTQs can be far off, in particular when you have a very low expressed gene with only a few non-zeroes. The likelihood-means seem reasonable for the zeroes, but explode for the few non-zeroes.

I investigated a bit (aided by Claude and Erik's write-up), and now think that I understand where the problem is. _The important realisation for me was that the Gaussian approximation of the likelihood can be obtained by doing a 2nd order Taylor expansion of the Poisson log-likelihood **around the optimum of the posterior** ($x_c^*$)_. I wrote it out to convince myself in the following section, and I discuss the consequence in the section after

### The Gaussian likelihood as a Taylor-expansion
For a given $v$, we have the following log-posterior (Eq. 27 of Erik's notes):
$$
L(x|v) = -\frac{C-1}{2} \log(v) + \sum_c \left[n_c x_c -N_c e^{x_c} - \frac{(x_c - \langle x \rangle)^2}{ 2v}\right]
$$
Now I followed the derivation to get the optimal $x_c^*$ and the Laplace approximation of this posterior, after which we divide out the Gaussian prior to get a Gaussian approximation of the likelihood.

With $n f^*_c = N_c e^{x^*_c}$, we expand $L(x|v)$ around $x^*$:
- 0th order: $L(x^*|v)$.
<!-- - 1st order: $\sum_c \left[(n_c - n f^*_c) - \frac{x^*_c - \langle x^* \rangle}{v}\right] (x_c - x^*_c) = 0$ (this linear term is zero of course because we expand around the optimum, but it is useful to keep it here). -->
- 1st order: 0 because we expand around the optimum
- 2nd order: $-\tfrac12 \sum_c n f^*_c\, (x_c - x^*_c)^2 - \frac{1}{2v} \sum_c \big[(x_c - \langle x \rangle) - (x^*_c - \langle x^* \rangle)\big]^2.$

So the Laplace approximation of the posterior is
$$
L(x|v) \approx L(x^*|v) - \tfrac12 \sum_c n f^*_c\, (x_c - x^*_c)^2 - \frac{1}{2v} \sum_c \big[(x_c - \langle x \rangle) - (x^*_c - \langle x^* \rangle)\big]^2.
$$
What helped me is to see that, since the prior is quadratic, we can also get a quadratic approximation of the posterior by only expanding the likelihood around $x^*_c$, so the following is equivalent:
$$
L(x|v) \approx \sum_c \left[(n_c - n f^*_c)(x_c - x^*_c) - \tfrac12\, n f^*_c\,(x_c - x^*_c)^2 - \frac{(x_c - \langle x \rangle)^2}{2v}\right] + \text{const}.
$$
It is confusing at first, because now it looks like there are still linear terms, but these drop out. Now I think Erik was saying this all along, but this way of writing the Laplace approx. finally made me realise that it is really clean to "divide out the prior" because we just leave out the last prior-related bit. We thus get:
$$
L_m(x|v) \approx \sum_c \left[(n_c - n f^*_c)(x_c - x^*_c) - \tfrac12\, n f^*_c\,(x_c - x^*_c)^2\right] + \text{const}.
$$
where $L_m$ stands for measurement-likelihood. Now we can again complete the squares to show that this gives us the mean and variance of the likelihood-Gaussian:
$$
L_m(x|v) \propto -\tfrac12 \sum_c n f^*_c\, (x_c - x^m_c)^2,
\qquad
\text{with}
\qquad
x^m_c = x^*_c + \frac{n_c - n f^*_c}{n f^*_c} = x^*_c + \frac{x^*_c - \langle x^* \rangle}{n f^*_c\, v},
$$
which is what Erik had as Eq. 44.

### What can go wrong
So I think the problems occur when our Taylor expansion for this loglikelihood is not accurate near the optimum of this loglikelihood. The important thing is that we expand around $x^*_c$ not around the likelihood-maximum $\tilde x_c = \log(n_c/N_c)$. We know that the Taylor expansion is only valide when $|x-x_c^*|$ small, but $x_c^m - x_c^* = \frac{n_c}{nf_c^*} -1$ may not be that small. Especially when the "expected count", $n f_c^*$ is much lower than the actual count $n_c$ (while the problem is less severe when $n_c=0$). This can happen when we have a gene with a very small fitted variance $v$.

To illustrate what happened, I had Claude made a figure of the prior, the Poisson likelihood, the posterior, and then the Gaussian-approximated likelihood.

<img src="../../../IMBA_non_git/collaborations/nimwegen/example_of_misfit.png" alt="Example of the misfit" width="50%">

### A possible solution
Maybe this is very naive, but I think we mostly need to do the Taylor-expansion of the Poisson-likelihood around the max-posterior $x_c^*$ for zero counts, because for these zero-counts, the Poisson-likelihood has the maximum at minus infinity. And indeed, for those zero counts, it works quite well.

However, for the non-zero counts, maybe it is better to do the Taylor-expansion of the Poisson-likelihood directly around its own optimum. We would then just get
$$
n_c x_c - N_c e^{x_c} \approx -\tfrac12\, n_c\, (x_c - \tilde x_c)^2 + \text{const},
$$
since the slope vanishes at $\tilde x_c = \log(n_c/N_c)$ and the curvature there is $-N_c e^{\tilde x_c} = -n_c$.
So the measurement would be $x^m_c = \log(n_c/N_c)$ with variance $1/n_c$, independent of the prior.