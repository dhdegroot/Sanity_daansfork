#include "gene_fit.hpp"
#include "fit_frac.hpp"

#include <algorithm>
#include <optional>
#include <cmath>
#include <stdexcept>
#include <vector>

// "SI" is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.
// "note" is E. van Nimwegen, "Getting proper error-bars in Sanity" (update_final_proper_error_bars.pdf), Section 2.
// Index letters at the end of a name give what an array runs over: c = cells, b = grid points of the gene variance v_g.
// Where the SI has n_g + 1 (Eqs. 16, 19, 26, 31, 44, 46) this code uses n_g.
// This is intentional: the prior over alpha_g differs from the one in the SI.
//
// The fit estimates the LTQs x_gc of one gene directly (note Section 2) rather than the log fold changes delta_gc.
// get_gene_expression_level (at the end of this file) has two stages:
//   1. One pass over the grid of candidate variances v_b. At each grid point fit_fractions solves for q_g(v_b) and
//      the fractions f*_gc (note Eq. 31), from which the log likelihood L(v_b) (note Eqs. 36-38), the mean LTQ
//      mu_g(v_b) (Eq. 33) and its variance sigma_mu^2(v_b) (Eq. 43) follow. Only q_g, L, mu_g and sigma_mu^2 are
//      kept per grid point.
//   2. Either one grid point is selected and the cell estimates are made there (MLE, MAP, EAP), or a second pass over
//      the grid makes the cell estimates at every grid point and averages them over the posterior of v_g (MARG,
//      note Eqs. 47-50).
// Both branches make the cell estimates with the same functions: reconstruct_fractions (fit_frac.cpp), which gives
// f*_gc back from q_g, compute_ltqs, which turns f*_gc into x*_gc and the sums over cells that the error bars need,
// and cell_error_bars, which gives the error bars and the measurement LTQ of one cell. The renormalised variances
// v_gc of Eq. 39 are never stored: they are computed in the first pass for their average and again in the second.
// MARG averages the cell estimates over the grid as it goes, so no estimates are stored per grid point.

GeneWorkspace::GeneWorkspace(int num_cells, int num_bins, VarianceMethod variance_method)
    : lambert_log_arg_base_c(num_cells), fractions_c(num_cells),
      log_weighted_total_b(num_bins), log_likelihood_b(num_bins), mean_ltq_b(num_bins), mean_ltq_error_sq_b(num_bins)
{
    if (variance_method == VarianceMethod::MARG)
    {
        bin_ltq_c.resize(num_cells);
        ltq_spread_c.resize(num_cells);
        measured_ltq_spread_c.resize(num_cells);
    }
}

// The sums over cells at one v_g that the likelihood (note Eq. 38) and the error bars (Eqs. 41, 43) need.
// With w_gc = n_g f*_gc v_g + 1, the renormalised variance of Eq. 39 is v_gc = v_g / w_gc, and
// <v> = v_g sum_inverse_curvature / C and v_g - <v> = v_g sum_shrinkage / C, both without cancellation.
struct CellSums
{
    double sum_ltq;                // sum_c x*_gc
    double sum_shrinkage;          // sum_c (w_gc - 1) / w_gc = sum_c n_g f*_gc v_g / (n_g f*_gc v_g + 1)
    double sum_inverse_curvature;  // sum_c 1 / w_gc
};

// x*_gc of every cell at one v_g (note Eq. 32), written to ltq_c, and the sums over cells of CellSums.
// log_gene_total = log(n_g), total_times_variance = n_g v_g.
static CellSums compute_ltqs(const std::vector<double>& fractions_c, double log_gene_total, double total_times_variance,
                             const std::vector<double>& log_cell_totals_c, double* ltq_c)
{
    const int num_cells = static_cast<int>(fractions_c.size());
    CellSums sums{0.0, 0.0, 0.0};
    for (int c = 0; c < num_cells; ++c)
    {
        ltq_c[c] = std::log(fractions_c[c]) + log_gene_total - log_cell_totals_c[c]; // x*_gc (Eq. 32)
        sums.sum_ltq += ltq_c[c];
        const double scaled_expected_count = total_times_variance * fractions_c[c]; // n_g f*_gc v_g
        const double curvature = scaled_expected_count + 1.0;                        // w_gc = n_g f*_gc v_g + 1
        sums.sum_shrinkage += scaled_expected_count / curvature;
        sums.sum_inverse_curvature += 1.0 / curvature;
    }
    return sums;
}

// sigma_mu^2 at one v_g (note Eq. 43): <v> / (C (1 - <v>/v_g)), written with the sums of CellSums.
static inline double mean_ltq_error_sq_from_sums(const CellSums& sums, double gene_variance, int num_cells)
{
    return gene_variance * sums.sum_inverse_curvature / (num_cells * sums.sum_shrinkage);
}

// The estimates of one cell at one v_g besides x*_gc.
struct CellErrorBars
{
    double ltq_error_sq;          // epsilon_gc^2 (note Eq. 41)
    double measured_ltq;          // x^m_gc (note Eq. 44)
    double measured_ltq_error_sq; // 1/(n_g f*_gc) (note Eq. 46)
};

// The error bars and the measurement LTQ of one cell at one v_g, from its fraction f*_gc, its x*_gc, the mean
// <x*> and the sums of the first pass. expected_count = n_g f*_gc, gene_variance = v_g.
static inline CellErrorBars cell_error_bars(double expected_count, double ltq, double mean_ltq, double gene_variance,
                                            const CellSums& sums)
{
    const double scaled_expected_count = expected_count * gene_variance;                   // n_g f*_gc v_g
    const double renormalised_variance = gene_variance / (scaled_expected_count + 1.0);    // v_gc (Eq. 39)
    CellErrorBars bars;
    // Eq. 41: v_gc + (1/C) v_gc^2 / (v_g - <v>), with C (v_g - <v>) = v_g sum_shrinkage.
    bars.ltq_error_sq = renormalised_variance + renormalised_variance * renormalised_variance / (gene_variance * sums.sum_shrinkage);
    bars.measured_ltq = ltq + (ltq - mean_ltq) / scaled_expected_count; // Eq. 44
    bars.measured_ltq_error_sq = 1.0 / expected_count;                   // Eq. 46
    return bars;
}

// What one grid point contributes to stage 1 besides q_g.
struct BinSummary
{
    double log_likelihood;     // L(v_b) (note Eq. 37)
    double mean_ltq;           // mu_g(v_b) = <x*> (note Eq. 33)
    double mean_ltq_error_sq;  // sigma_mu^2(v_b) (note Eq. 43)
};

// L(v_b), the log marginal likelihood at one grid point, from f*_gc and q_g at that point (note Eqs. 36-38), with
// the mean LTQ and its variance at that point, which come from the same sums over cells.
// Eq. 36 is the log posterior at the optimum x*_gc = log(f*_gc) + log(n_g) - log(N_c) (Eq. 32); Eq. 37 adds minus
// half the log determinant of the precision matrix of Eq. 35, which is Eq. 38.
// The quadratic term of Eq. 36 is summed around log(n_g) - q_g, which is <x*> at the exact optimum (Eq. 33), and
// then moved to the actual mean <x*>: sum_c (x_c - <x>)^2 = sum_c (x_c - m)^2 - C (<x> - m)^2 for any m.
// log1p is used for log(n_g f*_gc v_b + 1) because n_g f*_gc v_b is often much smaller than one.
// This L(v_b) equals the one of the 2021 SI (Eqs. 20, 32, 33) up to a constant that does not depend on v_b, so the
// posterior over the grid and the selected v_g are unchanged; see docs/likelihood_equivalence.md.
static BinSummary bin_log_likelihood(const std::vector<double>& umi_counts_c, const std::vector<double>& fractions_c,
                                     const std::vector<double>& log_cell_totals_c, double gene_total, double log_gene_total,
                                     double bin_variance, double log_weighted_total)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const double total_times_variance = gene_total * bin_variance; // n_g v_b
    const double centre = log_gene_total - log_weighted_total;     // log(n_g) - q_g, the mean <x*> up to the Newton tolerance (Eq. 33).
    CellSums sums{0.0, 0.0, 0.0};
    double sum_counts_times_ltq = 0.0;   // sum_c n_gc x*_gc, the data term of Eq. 36
    double sum_fractions = 0.0;          // sum_c f*_gc, which is one up to the Newton tolerance (Eq. 29)
    double sum_centred_ltq_sq = 0.0;     // sum_c (x*_gc - centre)^2
    double sum_log_curvature = 0.0;      // sum_c log(n_g f*_gc v_b + 1), the second term of Eq. 38
    for (int c = 0; c < num_cells; ++c)
    {
        const double ltq = std::log(fractions_c[c]) + log_gene_total - log_cell_totals_c[c]; // x*_gc (Eq. 32)
        sums.sum_ltq += ltq;
        sum_counts_times_ltq += umi_counts_c[c] * ltq;
        sum_fractions += fractions_c[c];
        const double deviation = ltq - centre;
        sum_centred_ltq_sq += deviation * deviation;
        const double scaled_expected_count = total_times_variance * fractions_c[c]; // n_g f*_gc v_b
        const double curvature = scaled_expected_count + 1.0;                        // w_gc
        sum_log_curvature += std::log1p(scaled_expected_count);
        sums.sum_shrinkage += scaled_expected_count / curvature;                     // C times the average in the last term of Eq. 38
        sums.sum_inverse_curvature += 1.0 / curvature;
    }
    BinSummary summary;
    summary.mean_ltq = sums.sum_ltq / num_cells; // <x*> (Eq. 33)
    summary.mean_ltq_error_sq = mean_ltq_error_sq_from_sums(sums, bin_variance, num_cells);
    const double mean_offset = summary.mean_ltq - centre; // <x*> - centre
    const double sum_ltq_deviation_sq = sum_centred_ltq_sq - num_cells * mean_offset * mean_offset; // sum_c (x*_gc - <x*>)^2
    // Eq. 36, the log posterior at the optimum.
    const double log_posterior = -0.5 * (num_cells - 1) * std::log(bin_variance)
                                 + sum_counts_times_ltq - gene_total * sum_fractions
                                 - sum_ltq_deviation_sq / (2.0 * bin_variance);
    // Eq. 38, log det(M) of the precision matrix of Eq. 35.
    const double log_det_precision = -num_cells * std::log(bin_variance) + sum_log_curvature
                                     + std::log(sums.sum_shrinkage / num_cells);
    summary.log_likelihood = log_posterior - 0.5 * log_det_precision; // Eq. 37
    return summary;
}

// P(v_b | n_g), the normalised posterior over the grid (SI Eq. 34), from L(v_b) and its maximum over the grid.
static void normalise_posterior(const std::vector<double>& log_likelihood_b, double max_log_likelihood,
                                std::vector<double>& variance_likelihood_b)
{
    const int num_bins = static_cast<int>(log_likelihood_b.size());
    double likelihood_sum = 0.0; // The normalisation, sum_b exp(L_b).
    for (int b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] = log_likelihood_b[b] - max_log_likelihood;
        variance_likelihood_b[b] = std::exp(variance_likelihood_b[b]);
        likelihood_sum += variance_likelihood_b[b];
    }
    for (int b = 0; b < num_bins; b++)
    {
        variance_likelihood_b[b] /= likelihood_sum;
    }
}

// The grid point whose v_b is reported by the MLE, MAP and EAP methods.
static int select_variance_bin(VarianceMethod variance_method, const std::vector<double>& variance_likelihood_b,
                               const std::vector<double>& variance_grid_b, int max_likelihood_bin)
{
    const int num_bins = static_cast<int>(variance_grid_b.size());
    int selected_bin = 0;
    if (variance_method == VarianceMethod::MLE)
    {
        selected_bin = max_likelihood_bin;
    }
    else if (variance_method == VarianceMethod::MAP)
    {
        // The maximum of the posterior density in v_g under a 1/v_g (Jeffreys) prior: the posterior mass
        // of a grid point divided by its width. The grid is equally spaced in log(v_g), so the width of a grid
        // point is proportional to v_b, and EAP and MARG, which weigh the grid points by their mass, assume the same prior.
        double max_posterior = -1e+100;
        for (int b = 0; b < num_bins; b++)
        {
            const double posterior = variance_likelihood_b[b] / variance_grid_b[b];
            if (posterior > max_posterior)
            {
                max_posterior = posterior;
                selected_bin = b;
            }
        }
    }
    else if (variance_method == VarianceMethod::EAP)
    {
        double posterior_mean_variance = 0.0;
        for (int b = 0; b < num_bins; b++)
        {
            posterior_mean_variance += variance_likelihood_b[b] * variance_grid_b[b];
        }
        double min_distance = 1e+100;
        for (int b = 0; b < num_bins; b++)
        {
            const double distance = std::fabs(variance_grid_b[b] - posterior_mean_variance);
            if (distance < min_distance)
            {
                min_distance = distance;
                selected_bin = b;
            }
        }
    }
    else
    {
        throw std::logic_error("select_variance_bin: the marginalising method does not select a grid point");
    }
    return selected_bin;
}

// MLE, MAP and EAP: report v_b, and every estimate at the selected grid point (note Section 2.10, "once we picked a
// 'best' v_g"). Two passes over the cells: x*_gc and the sums (compute_ltqs), then the error bars and x^m_gc.
static void estimate_at_selected_bin(int selected_bin, const std::vector<double>& umi_counts_c, double gene_total, double log_gene_total,
                                     const std::vector<double>& variance_grid_b, const std::vector<double>& log_cell_totals_c,
                                     GeneWorkspace& workspace, RowComputation& result)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const double gene_variance = variance_grid_b[selected_bin];
    const double log_weighted_total = workspace.log_weighted_total_b[selected_bin];
    const std::vector<double>& fractions_c = workspace.fractions_c;
    result.gene_variance = gene_variance;
    reconstruct_fractions(workspace.fractions_c.data(), umi_counts_c, gene_total, gene_variance, num_cells, log_weighted_total, log_cell_totals_c.data());
    const CellSums sums = compute_ltqs(fractions_c, log_gene_total, gene_total * gene_variance, log_cell_totals_c, result.ltq_c.data());
    result.mean_ltq = sums.sum_ltq / num_cells; // mu_g = <x*> (Eq. 33)
    result.mean_ltq_error_sq = mean_ltq_error_sq_from_sums(sums, gene_variance, num_cells); // Eq. 43
    for (int c = 0; c < num_cells; ++c)
    {
        const CellErrorBars bars = cell_error_bars(gene_total * fractions_c[c], result.ltq_c[c], result.mean_ltq, gene_variance, sums);
        result.ltq_error_sq_c[c] = bars.ltq_error_sq;
        result.measured_ltq_c[c] = bars.measured_ltq;
        result.measured_ltq_error_sq_c[c] = bars.measured_ltq_error_sq;
    }
}

// MARG: make the cell estimates at every grid point and average everything over the posterior P(v_b | n_g)
// (note Eqs. 47-50). The reported v_g is the posterior mean of the variance.
// The averages over the grid are accumulated one grid point at a time with West's weighted update of the mean and
// the spread, so only the x*_gc of the current grid point are held. x^m_gc and its error bar are averaged in the
// same way as x*_gc (Eqs. 47, 48 with 1/(n_g f*_gc(v_b)) as the variance at each grid point); the note does not
// define a marginal x^m_gc.
static void marginalise_over_variance(const std::vector<double>& umi_counts_c, double gene_total, double log_gene_total,
                                      const std::vector<double>& variance_grid_b, const std::vector<double>& log_cell_totals_c,
                                      GeneWorkspace& workspace, RowComputation& result)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const int num_bins = static_cast<int>(variance_grid_b.size());
    const std::vector<double>& variance_likelihood_b = result.variance_likelihood_b;
    const std::vector<double>& fractions_c = workspace.fractions_c;
    double* bin_ltq_c = workspace.bin_ltq_c.data();
    // sum_b P(v_b) (x*_gc(v_b) - <x*_gc>)^2 so far, in West's update, and the same for x^m_gc.
    std::vector<double>& ltq_spread_c = workspace.ltq_spread_c;
    std::vector<double>& measured_ltq_spread_c = workspace.measured_ltq_spread_c;
    std::vector<double>& ltq_c = result.ltq_c;                                       // x^marg_gc = sum_b P(v_b) x*_gc(v_b) (Eq. 47)
    std::vector<double>& ltq_error_sq_c = result.ltq_error_sq_c;                     // Eq. 48
    std::vector<double>& measured_ltq_c = result.measured_ltq_c;                     // sum_b P(v_b) x^m_gc(v_b)
    std::vector<double>& measured_ltq_error_sq_c = result.measured_ltq_error_sq_c;   // as Eq. 48 for x^m_gc
    std::fill(ltq_c.begin(), ltq_c.end(), 0.0);
    std::fill(ltq_spread_c.begin(), ltq_spread_c.end(), 0.0);
    std::fill(ltq_error_sq_c.begin(), ltq_error_sq_c.end(), 0.0);
    std::fill(measured_ltq_c.begin(), measured_ltq_c.end(), 0.0);
    std::fill(measured_ltq_spread_c.begin(), measured_ltq_spread_c.end(), 0.0);
    std::fill(measured_ltq_error_sq_c.begin(), measured_ltq_error_sq_c.end(), 0.0);

    // Second pass over the grid: the cell estimates at every grid point with nonzero weight.
    double weight_so_far = 0.0; // sum of P(v_b) over the grid points seen so far.
    for (int b = 0; b < num_bins; ++b)
    {
        const double weight = variance_likelihood_b[b];
        if (weight == 0.0) continue;
        weight_so_far += weight;
        const double bin_variance = variance_grid_b[b];
        reconstruct_fractions(workspace.fractions_c.data(), umi_counts_c, gene_total, bin_variance, num_cells, workspace.log_weighted_total_b[b], log_cell_totals_c.data());
        const CellSums sums = compute_ltqs(fractions_c, log_gene_total, gene_total * bin_variance, log_cell_totals_c, bin_ltq_c);
        const double bin_mean_ltq = sums.sum_ltq / num_cells; // mu_g(v_b)
        for (int c = 0; c < num_cells; ++c)
        {
            const CellErrorBars bars = cell_error_bars(gene_total * fractions_c[c], bin_ltq_c[c], bin_mean_ltq, bin_variance, sums);
            const double deviation = bin_ltq_c[c] - ltq_c[c];
            ltq_c[c] += (weight / weight_so_far) * deviation;
            ltq_spread_c[c] += weight * deviation * (bin_ltq_c[c] - ltq_c[c]);
            ltq_error_sq_c[c] += weight * bars.ltq_error_sq;
            const double measured_deviation = bars.measured_ltq - measured_ltq_c[c];
            measured_ltq_c[c] += (weight / weight_so_far) * measured_deviation;
            measured_ltq_spread_c[c] += weight * measured_deviation * (bars.measured_ltq - measured_ltq_c[c]);
            measured_ltq_error_sq_c[c] += weight * bars.measured_ltq_error_sq;
        }
    }
    // Eq. 48: <epsilon^2(v_b)> + <(x*(v_b) - x^marg)^2>; P(v_b) sums to one.
    for (int c = 0; c < num_cells; ++c)
    {
        ltq_error_sq_c[c] += ltq_spread_c[c];
        measured_ltq_error_sq_c[c] += measured_ltq_spread_c[c];
    }

    // The mean LTQ averaged over the grid (Eq. 49), and its variance (Eq. 50) with sigma_mu^2(v_b) from stage 1.
    double mean_ltq = 0.0;
    for (int b = 0; b < num_bins; b++)
    {
        mean_ltq += variance_likelihood_b[b] * workspace.mean_ltq_b[b];
    }
    double mean_ltq_error_sq = 0.0;
    for (int b = 0; b < num_bins; b++)
    {
        const double deviation = workspace.mean_ltq_b[b] - mean_ltq;
        mean_ltq_error_sq += variance_likelihood_b[b] * (workspace.mean_ltq_error_sq_b[b] + deviation * deviation);
    }
    result.mean_ltq = mean_ltq;
    result.mean_ltq_error_sq = mean_ltq_error_sq;

    // The reported v_g is the posterior mean of the variance.
    double gene_variance = 0.0;
    for (int b = 0; b < num_bins; ++b)
    {
        gene_variance += variance_grid_b[b] * variance_likelihood_b[b];
    }
    result.gene_variance = gene_variance;
}

RowComputation get_gene_expression_level(const std::vector<double> &umi_counts_c, double gene_total, const std::vector<double> &variance_grid_b,
                                         VarianceMethod variance_method, const std::vector<double>& log_cell_totals_c,
                                         double sum_cell_totals, GeneWorkspace& workspace)
{
    const int num_cells = static_cast<int>(umi_counts_c.size());
    const int num_bins = static_cast<int>(variance_grid_b.size());
    const double log_gene_total = std::log(gene_total); // log(n_g), the constant of Eq. 32.

    // Stage 1: fit every grid point, keeping q_g(v_b), L(v_b), mu_g(v_b) and sigma_mu^2(v_b).
    std::optional<double> previous_log_weighted_total; // q_g of the previous grid point, the starting guess for the next fit.
    double max_log_likelihood = -1e+100; // Maximum of L(v_b) over the grid.
    int max_likelihood_bin = 0; // Grid point at which L(v_b) is maximal (the MLE).
    for (int b = 0; b < num_bins; ++b)
    {
        const double bin_variance = variance_grid_b[b];
        // We fit log_weighted_total = q_g(v_b) (SI Eq. 22; note Eq. 31); fit_fractions also fills fractions_c with f*_gc.
        const double log_weighted_total = fit_fractions(workspace.fractions_c.data(), umi_counts_c, gene_total, bin_variance, num_cells, previous_log_weighted_total, log_cell_totals_c.data(), sum_cell_totals, workspace.lambert_log_arg_base_c.data());

        previous_log_weighted_total = log_weighted_total;
        workspace.log_weighted_total_b[b] = log_weighted_total;

        // Log marginal likelihood L(v_b) at this grid point, stored per grid point, while tracking the maximum for the MLE.
        const BinSummary summary = bin_log_likelihood(umi_counts_c, workspace.fractions_c, log_cell_totals_c, gene_total, log_gene_total, bin_variance, log_weighted_total);
        workspace.log_likelihood_b[b] = summary.log_likelihood;
        workspace.mean_ltq_b[b] = summary.mean_ltq;
        workspace.mean_ltq_error_sq_b[b] = summary.mean_ltq_error_sq;
        if (summary.log_likelihood > max_log_likelihood)
        {
            max_log_likelihood = summary.log_likelihood;
            max_likelihood_bin = b;
        }
    }

    RowComputation result;
    result.variance_likelihood_b.resize(num_bins);
    result.ltq_c.resize(num_cells);
    result.ltq_error_sq_c.resize(num_cells);
    result.measured_ltq_c.resize(num_cells);
    result.measured_ltq_error_sq_c.resize(num_cells);
    normalise_posterior(workspace.log_likelihood_b, max_log_likelihood, result.variance_likelihood_b);

    // Stage 2: the reported estimates, either at one grid point or averaged over the grid.
    if (variance_method == VarianceMethod::MARG)
    {
        marginalise_over_variance(umi_counts_c, gene_total, log_gene_total, variance_grid_b, log_cell_totals_c, workspace, result);
    }
    else
    {
        // The grid point used depends on the variance method: the maximum likelihood (MLE), the posterior maximum (MAP), or the point closest to the posterior mean of v_g (EAP).
        const int selected_bin = select_variance_bin(variance_method, result.variance_likelihood_b, variance_grid_b, max_likelihood_bin);
        estimate_at_selected_bin(selected_bin, umi_counts_c, gene_total, log_gene_total, variance_grid_b, log_cell_totals_c, workspace, result);
    }
    return result;
}
