#ifndef SANITY_OPTIONS_HPP
#define SANITY_OPTIONS_HPP

#include <string>

#include "gene_fit.hpp"

// Which LTQ estimates the default output holds (--ltq-estimate); -e 1 always writes both kinds.
enum class LtqEstimate
{
    Posterior = 0, // x*_gc with the posterior error bars epsilon_gc (note Eqs. 32, 41); the default.
    Measured = 1,  // x^m_gc with the measurement error bars 1/sqrt(n_g f*_gc) (note Eqs. 44, 46).
    Both = 2       // Both pairs of files.
};

// Everything the command line of Sanity can set, with its default.
struct SanityOptions
{
    std::string in_file;                    // The count matrix (-f).
    std::string in_file_extension;          // "mtx" or "tsv", derived from in_file.
    std::string gene_name_file = "none";    // Gene names of an MTX input (-mtx_genes).
    std::string cell_name_file = "none";    // Cell names of an MTX input (-mtx_cells).
    std::string out_folder = "./";          // Output folder, always ending in '/' (-d).
    int num_threads = 4;                    // -n
    bool extended_output = false;           // -e
    double variance_min = 0.001;            // v_min, the smallest candidate gene variance (SI Section S1.1.4).
    double variance_max = 50.0;             // v_max, the largest candidate gene variance.
    int num_bins = 160;                     // B, the number of grid points of the gene variance between v_min and v_max.
    VarianceMethod variance_method = VarianceMethod::MAP; // How to estimate v_g (see gene_fit.hpp).
    bool gzip_output = false;               // -gz: compress the numerical text files; overridden by -npy.
    bool npy_output = false;                // -npy: write the numerical output files as NPY; overrides -gz.
    bool bonsai_output_mode = false;        // --bonsai-output-mode: the Bonsai file set; overrides -e, -gz and -npy.
    LtqEstimate ltq_estimate = LtqEstimate::Posterior; // --ltq-estimate: what the default files hold; Bonsai mode defaults to Measured.
    MeasuredExpansion measured_expansion = MeasuredExpansion::PosteriorMode; // --measured-expansion: the formula of the measurement LTQs (gene_fit.hpp).
};

enum class ParseResult
{
    Continue,
    VersionRequested,
    HelpRequested,
    Error
};

// Fills `options` from the command line. Errors are logged before Error is returned.
ParseResult parse_options(int argc, char **argv, SanityOptions &options);

void show_usage();

// The name of a variance method as written to sanity_command.txt ("MARG", "MLE", "MAP" or "EAP").
const char *variance_method_name(VarianceMethod variance_method);

// The name of an LTQ estimate as written to sanity_command.txt ("posterior", "measured" or "both").
const char *ltq_estimate_name(LtqEstimate ltq_estimate);

// The name of a measured expansion as written to sanity_command.txt ("posterior-mode" or "likelihood-maximum").
const char *measured_expansion_name(MeasuredExpansion measured_expansion);

#endif
