#ifndef SANITY_GENE_FIT_HPP
#define SANITY_GENE_FIT_HPP

#include <vector>

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// "note" is E. van Nimwegen, "Getting proper error-bars in Sanity" (update_final_proper_error_bars.pdf), Section 2.
// Index letters at the end of a name give what an array runs over: g = genes, c = cells, b = grid points of the gene variance v_g.

// How the gene variance v_g is estimated. The values index the names returned by variance_method_name (sanity_options.hpp).
enum class VarianceMethod
{
    MARG = 0, // Average every estimate over the posterior of v_g.
    MLE = 1,  // Use the v_g that maximises the likelihood.
    MAP = 2,  // Use the v_g that maximises the posterior under a 1/v_g prior (the default).
    EAP = 3   // Use the grid point closest to the posterior mean of v_g.
};

// Where the Poisson log-likelihood of a cell is expanded to second order to get its measurement LTQ x^m_gc and
// error bar (--measured-expansion). The values index the names returned by measured_expansion_name (sanity_options.hpp).
enum class MeasuredExpansion
{
    PosteriorMode = 0,     // Around the posterior mode x*_gc for every cell: x^m_gc of note Eq. 44, error bar of Eq. 46 (the default).
    LikelihoodMaximum = 1  // For n_gc > 0 around the likelihood maximum log(n_gc/N_c): x^m_gc = log(n_gc/N_c), error bar 1/sqrt(n_gc);
                           // zero counts as PosteriorMode (their likelihood has no maximum). See docs/measured_ltq_linearisation.md.
};

// The results of the fit of one gene: the log transcription quotients (LTQs) x_gc of the gene in every cell.
struct RowComputation
{
    std::vector<double> ltq_c;                   // x*_gc, the posterior LTQ of each cell (note Eq. 32; MARG: Eq. 47).
    std::vector<double> ltq_error_sq_c;          // epsilon_gc^2, the posterior variance of x_gc (note Eq. 41; MARG: Eq. 48).
    std::vector<double> measured_ltq_c;          // x^m_gc, the measurement LTQ: x*_gc with the prior divided out (note Eq. 44), or
                                                 // log(n_gc/N_c) for n_gc > 0 with MeasuredExpansion::LikelihoodMaximum.
    std::vector<double> measured_ltq_error_sq_c; // The squared measurement error bar of x^m_gc: 1/(n_g f*_gc) (note Eq. 46), or
                                                 // 1/n_gc for n_gc > 0 with MeasuredExpansion::LikelihoodMaximum.
    double mean_ltq;                             // mu_g = <x*_g>, the mean LTQ of the gene over the cells (note Eq. 33; MARG: Eq. 49).
    double mean_ltq_error_sq;                    // sigma_mu^2, the variance of mu_g (note Eq. 43; MARG: Eq. 50).
    double gene_variance;                        // v_g, the prior variance of the LTQs across cells (note Eq. 25).
    std::vector<double> variance_likelihood_b;   // P(v_b | n_g), the posterior of v_g on the grid of variances (SI Eq. 34).
};

// Scratch space of one thread, allocated once and reused for every gene the thread fits.
struct GeneWorkspace
{
    GeneWorkspace(int num_cells, int num_bins, VarianceMethod variance_method);

    std::vector<double> lambert_log_arg_base_c;  // Q_c, the per-cell constants of the fit at one v_g (SI Eq. 24), used by fit_fractions.
    std::vector<double> fractions_c;             // f*_gc, the fraction of the gene's counts expected in each cell (note Eq. 31).
    std::vector<double> log_weighted_total_b;    // q_g(v_b) at every grid point (SI Eq. 22; note Eq. 31).
    std::vector<double> log_likelihood_b;        // L(v_b), the log marginal likelihood at every grid point (note Eq. 37).
    std::vector<double> mean_ltq_b;              // mu_g(v_b) = <x*_g(v_b)> at every grid point (note Eq. 33).
    std::vector<double> mean_ltq_error_sq_b;     // sigma_mu^2(v_b) at every grid point (note Eq. 43).
    std::vector<double> bin_ltq_c;               // x*_gc at the current grid point; MARG only.
    std::vector<double> ltq_spread_c;            // The spread of x*_gc(v_b) over the grid seen so far, in West's update; MARG only.
    std::vector<double> measured_ltq_spread_c;   // The same spread of x^m_gc(v_b); MARG only.
};

// Fits one gene: the LTQ of every cell with its error bar, the measurement LTQs, the mean LTQ and the gene variance.
// umi_counts_c = n_gc, the counts of this gene in each cell (SI Eq. 13); gene_total = n_g, their sum (SI Eq. 15).
// log_cell_totals_c = log(N_c), the log of the total counts of each cell (SI Eq. 9); sum_cell_totals = sum_c N_c.
// variance_grid_b = v_b, the candidate values of the gene variance (SI Section S1.1.4).
// measured_expansion chooses the formula of the measurement LTQs; nothing else depends on it.
RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, double gene_total, const std::vector<double> &variance_grid_b,
                                         VarianceMethod variance_method, MeasuredExpansion measured_expansion,
                                         const std::vector<double>& log_cell_totals_c, double sum_cell_totals, GeneWorkspace& workspace);

#endif
