# Sanity

**Sa**mpling **N**oise based **I**nference of **T**ranscription Activit**Y** : Filtering of Poisson noise on a single-cell RNA-seq UMI count matrix

Single-cell RNA sequencing normalization algorithm presented in the publication [Bayesian inference of gene expression states from single-cell RNA-seq data -
J Breda, M Zavolan, E van Nimwegen - Nature Biotechnology, 2021](https://www.nature.com/articles/s41587-021-00875-x).

Sanity infers the log expression levels *x<sub>gc</sub>* of gene *g* in cell *c* by filtering out 
the Poisson noise on the UMI count matrix *n<sub>gc</sub>* of gene *g* in cell *c*.

The fit estimates the log transcription quotients (LTQs) *x<sub>gc</sub>* directly, with error bars from the full posterior over the cells, following Section 2 of E. van Nimwegen's note "Getting proper error-bars in Sanity".



### Reproducibility
The raw UMI count and normalized datasets mentioned in benchmarking in the associated [publication](https://www.nature.com/articles/s41587-021-00875-x) are available on [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4009187.svg)](https://zenodo.org/record/4009187). Files are named `[dataset name]_UMI_counts.txt.gz` and `[dataset name]_[tool name]_normalization.txt.gz`.

The scripts used for running the benchmarked normalization methods and for making the figures of the preprint are in the reproducibility folder.

## Input

* `-f`: UMI count matrix: *(N<sub>g</sub> x N<sub>c</sub>)* matrix with *N<sub>g</sub>* the number of genes and *N<sub>c</sub>* the number of cells. (`'path/to/count_matrix.tsv'`)

| GeneID | Cell 1 | Cell 2 | Cell 3 | ... |
|:-------|:------:|:------:|:------:|:---:|
| Gene 1 | 1.0 | 2.0 | 0.0 | ... |
| Gene 2 | 6.0 | 3.0 | 1.0 | ... |
| ... | ... | ... | ... | ... |

**The count matrix must be tab-separated (TSV).**
The tab character is the only field separator Sanity accepts: comma-separated (CSV) and space-separated files are *not* supported and will be misread.
Convert them to TSV before running Sanity, for example with `pandas.DataFrame.to_csv(path, sep='\t')` or `awk 'BEGIN{FS=",";OFS="\t"}{$1=$1;print}' in.csv > out.tsv`.
Because tabs alone delimit fields, gene and cell names may contain spaces and commas.

The first line is the header with the cell names.
It may either start with a label for the gene-ID column (the same number of fields as a data row) or start straight away with the first cell name (one field fewer, as written by `pandas.DataFrame.to_csv`).
Every data row must have exactly one gene ID followed by one value per cell; Sanity stops with an error if a row has too few or too many fields.

The count matrix can be compressed with gzip (with `.tsv.gz` extension).

* `-f`: (Alternatively) Matrix Market File Format: Sparse matrix of UMI counts. Automatically recognized by `.mtx` extension of the input file. Example: `matrix.mtx` by cellranger 2.1.0 and 3.1.0 (10x Genomics). (`'path/to/text_file.mtx'`)
**Important**: MTX file should be sorted by row (gene) indices! If not then you can use `sort_mtx_by_row.py` script in the `Sanity/scripts` folder to sort it. See [Sorting MTX file by row (gene) indices](#sorting-mtx-file-by-row-gene-indices).
The MTX file can be compressed with gzip (with `.gz` extension).
	* `-mtx_genes`: (optional) Gene ID file: text file with one gene ID per line. The order of gene IDs should match the order of genes in the count matrix. Examples: `genes.tsv` by cellranger 2.1.0 and `features.tsv` by cellranger 3.1.0 (10x Genomics). (`'path/to/text_file'`)
	* `-mtx_cells`: (optional) Cell ID file: text file with one cell ID per line. The order of cell IDs should match the order of cells in the count matrix. Examples: `barcodes.tsv` by cellranger 2.1.0 and 3.1.0 (10x Genomics).  (`'path/to/text_file'`)
* `-d`: (optional) Destination folder (`'path/to/output/folder'`, default: `cwd`)
* `-n`: (optional) Number of threads (integer, default: `4`). The same setting controls input indexing and inference.

  Uncompressed, row-sorted MTX files with nonnegative integral counts are indexed in parallel using bounded input buffers; no index cache is used. Counts written as `2.0` or `2e0` also qualify.

  Fractional counts, gzip input, and counts or totals outside the exact integer range of double precision use the existing serial MTX reader. If an unsupported count is encountered during parallel indexing, the input scan restarts with that reader. The indexing log reports the reader, actual thread count and elapsed time. TSV input retains its existing reader.
  During inference, each worker reuses its count buffer. Plain MTX rows are read through a bounded buffer with direct numeric parsing; fractional counts are supported. TSV and gzip retain their existing numeric conversion.

* `-e`: (optional) Write every output file (Boolean, `'true', 'false', '1'` or `'0'`, default: `false`). See [Output](#output).
* `-ltq`/`--ltq-estimate`: (optional) Which per-cell LTQ estimates are written: `posterior` (default; *x\*<sub>gc</sub>* with the posterior error bars), `measured` (*x<sup>m</sup><sub>gc</sub>* with the measurement error bars, the default in Bonsai output mode) or `both`. With `-e 1` both kinds are always written. See [Output](#output).
* `-gz`: (optional flag, no value) Compress the numerical text output files with gzip (`.txt.gz`). Overridden by `-npy`. See [Output formats](#output-formats-optional).
* `-npy`: (optional flag, no value) Write the numerical output files as NumPy arrays (`.npy`; float32 matrices, float64 per-gene vectors) instead of text. Overrides `-gz`. See [Output formats](#output-formats-optional).
* `--bonsai-output-mode`: (optional flag, no value) The preset for Bonsai: `-npy -v_m EAP --ltq-estimate measured`. See [Bonsai output mode](#bonsai-output-mode-optional).
* `-v_m`: (optional, expert-user-only) Choose the method to estimate gene-variances $v_g$. In the MAP, EAP, MLE-options, one value for $v_g$ is fixed, and the corresponding gene expression estimates are returned, in the MARG-option, the gene expression estimates are obtained by marginalizing $v_g$. The options are:
    * MAP (**default**): Use the maximum a posteriori estimate for $v_g$.
    * EAP: expected value of $v_g$ over the posterior.
    * MLE: maximum likelihood estimate for $v_g$.
    * MARG: (original method) The reported $v_g$ is the EAP estimate, but unlike the EAP option above, this value is not fixed when computing the gene expression estimates — instead, the gene expression estimates are obtained by marginalizing over the posterior of $v_g$.
* `-vmin/-vmax`: (optional, expert-user-only) Minimal and maximal considered values of the variance in log transcription quotients (double, default: *v<sub>min</sub>=*`0.001` *v<sub>max</sub>=*`50`)
* `-nbin`: (optional, expert-user-only) Number of bins for the variance in log transcription quotients (integer, default: `160`)

### How the fit works

For every gene, Sanity fits the optimal cell fractions *f\*<sub>gc</sub>* and the log marginal likelihood *L(v<sub>g</sub>)*
at every point of the variance grid, keeping only four numbers per grid point. For MAP, MLE and EAP one grid point is
then selected and all cell estimates are made there, in two passes over the cells and without per-cell scratch arrays:
the posterior LTQs *x\*<sub>gc</sub> = log(f\*<sub>gc</sub>) + log(n<sub>g</sub>) − log(N<sub>c</sub>)*, their posterior
error bars *&epsilon;<sub>gc</sub>* (the diagonal of the posterior covariance of the Laplace approximation, which includes
the uncertainty on the gene mean), the mean LTQ *&mu;<sub>g</sub> = &lang;x\*&rang;* with its error bar, and the
measurement LTQs *x<sup>m</sup><sub>gc</sub>* (the posterior with the prior divided out) with error bars
*1/sqrt(n<sub>g</sub> f\*<sub>gc</sub>)*. MARG makes the same estimates at every grid point and averages them over the
posterior of *v<sub>g</sub>*, bin by bin. The formulas are those of Section 2 of the note "Getting proper error-bars in
Sanity" (Eqs. 31–50), and the selected *v<sub>g</sub>* is the same as in version 2.0. Zero-count cells get the same posterior error bar as every other cell; the earlier asymmetric
zero-count error bar is gone. Runs with different `-n` give bit-identical results.

Text output uses fixed six-decimal formatting. Each inference worker formats
one gene row, and complete rows are written in gene order. Formatting uses the
same `-n` worker count as inference; there is no separate writer thread setting.
The number of formatted rows held in memory is bounded by the worker count.

## Output

Every file name says what the file holds, whatever the options of the run. Sanity reports two kinds of LTQ estimates
per cell. The **posterior** estimates *x\*<sub>gc</sub>* with error bars *&epsilon;<sub>gc</sub>* are the best
estimates of the LTQs given the counts and the fitted prior of the gene (mean *&mu;<sub>g</sub>*, variance
*v<sub>g</sub>*). The **measurement** estimates *x<sup>m</sup><sub>gc</sub>* with error bars
*1/sqrt(n<sub>g</sub> f\*<sub>gc</sub>)* are the posterior estimates with that prior divided out, for downstream
methods such as Bonsai that combine the measurements with a prior of their own; their error bars are the pure sampling
noise and can be very large for cells with few expected counts. `--ltq-estimate` selects which kind is written
(`posterior` unless `--bonsai-output-mode` is given, `measured`, or `both`); `-e 1` writes both kinds and every other
file.

The measurement estimates are the Gaussian that a second-order expansion of the Poisson log-likelihood of each cell
around *x\*<sub>gc</sub>* makes of it (Eqs. 44 and 46 of the note). They are a local expansion, not a point estimate:
for a count far above the count expected at the gene mean, *x<sup>m</sup><sub>gc</sub>* lies far from any plausible
LTQ (it can even be positive) with a correspondingly wide error bar. Read as a precision and a slope they are what a
downstream method with a prior of its own needs, and that is how Bonsai uses them; do not plot or threshold
`ltqs_measurement_approx` as a value.

The per-cell matrices are TSV files with a header line of the cell IDs and one line per gene that starts with the
gene ID:

  | GeneID | Cell 1 | Cell 2 | Cell 3 | ... |
  |:-------|:------:|:------:|:------:|:---:|
  | Gene 1 | -13.7227 | -13.722 | -13.729 | ... |
  | Gene 2 |  -9.96744 | -10.2522 | -10.1453 | ... |
  | ... | ... | ... | ... | ... |

* `ltqs_posterior.txt`: the posterior estimates *x\*<sub>gc</sub>* of the log-transcription quotients (LTQs). The LTQ *x<sub>gc</sub>* of gene *g* in cell *c* is the logarithm of the fraction of mRNAs in cell *c* that belong to gene *g*, so that *&Sigma;<sub>g</sub> exp(x<sub>gc</sub>) = 1* for each cell *c*; the estimated number of mRNAs of gene *g* in cell *c* is *exp(x<sub>gc</sub>)* times the estimated total number of mRNAs in the cell. Written with `--ltq-estimate posterior` (the default), `both` or `-e 1`.
* `ltqs_posterior_error_bars.txt`: the error bars *&epsilon;<sub>gc</sub>* (standard deviations) on the posterior estimates. They are the diagonal of the posterior covariance over the cells and include the uncertainty on the gene mean *&mu;<sub>g</sub>*, which is shared by all cells of the gene. Written together with `ltqs_posterior.txt`.
* `ltqs_measurement_approx.txt` and `ltqs_measurement_approx_error_bars.txt`: the measurement estimates *x<sup>m</sup><sub>gc</sub>* and their error bars *1/sqrt(n<sub>g</sub> f\*<sub>gc</sub>)*, in the same layout. Written with `--ltq-estimate measured`, `both`, `-e 1` or `--bonsai-output-mode`.

The per-gene vectors, one value per line in the order of `geneID.txt`, are always written:

* `gene_means.txt`: the mean LTQ *&mu;<sub>g</sub> = &lang;x\*<sub>gc</sub>&rang;* of each gene over all cells.
* `gene_means_error_bars.txt`: the error bars (standard deviations) on *&mu;<sub>g</sub>*.
* `gene_variances.txt`: the estimated variance *v<sub>g</sub>* of the LTQs *x<sub>gc</sub>* across cells. These variances are different, and generally larger, than the variance of the estimates in `ltqs_posterior.txt`, because they take the uncertainty on each *x<sub>gc</sub>* into account; when estimates of true gene expression variability are needed, use this file.

Always written as well: `geneID.txt` and `cellID.txt`, the gene and cell IDs in the order of the rows and columns,
and `sanity_command.txt`, which records the timestamp, the Sanity version, `# Method:` (the `-v_m` method),
`# LTQ estimate:` (`posterior`, `measured` or `both`), `# Output mode: bonsai` when `--bonsai-output-mode` was given,
and the command line.

## Extended output (optional)

`-e 1` writes both kinds of per-cell estimates, every per-gene vector, and in addition:

* `gene_variances_likelihood.txt`: the posterior distribution of each gene's true variance in log-expression. For the numerical calculation of this distribution, the variance is a priori assumed to lie in the range *[v<sub>min</sub>,v<sub>max</sub>]*, and the distribution is evaluated on a grid of *N<sub>b</sub>* values spaced uniformly on a logarithmic scale. The file contains a matrix with likelihoods for each gene *g* and each grid point *b*. The first line lists the grid values *v<sub>b</sub> = v<sub>min</sub> exp(b &middot; &Delta;v)* for *b = 0,...,N<sub>b</sub>-1* with *&Delta;v = log(v<sub>max</sub>/v<sub>min</sub>)/(N<sub>b</sub>-1)*, so the first and last values are *v<sub>min</sub>* and *v<sub>max</sub>* themselves; the following lines correspond to the likelihood values of genes.

  | | | | | |
  |:-------|:------:|:------:|:------:|:---:|
  | Variance | 0.001000 | 0.001070 | 0.001146 | ... |
  | Gene 1 | 0.018 | 0.019 | 0.020 | ... |
  | Gene 2 | 0.0006 | 0.0051 | 0.0031 | ... |
  | ... | ... | ... | ... | ... |

## Output formats (optional)

`-e` and `--ltq-estimate` choose which results are written; `-gz` and `-npy` only change how they are written. Both combine with `-e`; `-npy` overrides `-gz`.

* `-gz`: every numerical text file gets a `.gz` extension and is gzip-compressed; `geneID.txt` and `cellID.txt` stay plain text. `Sanity_distance` and `Sanity_gene_correlation` read uncompressed text only.
* `-npy`: every numerical file is written as an NPY array with the same base name, e.g. `ltqs_posterior.npy`. The per-cell matrices, of shape `(genes, cells)`, and `gene_variances_likelihood.npy`, of shape `(genes, bins)`, are float32, which carries the six decimals of the text output; the per-gene vectors, of shape `(genes, 1)`, are float64. With `-e`, `gene_variances_grid.npy` (float64, shape `(bins, 1)`) holds the variance grid that heads the text likelihood table. Error bars are standard deviations, as in the text files.
* With both flags, `-npy` wins: the NPY files are not compressed, so that they can still be memory-mapped.

| Files | default (`-ltq posterior`) | `-ltq measured` | `-ltq both` | `-e 1` (any `-ltq`) |
|:--|:--|:--|:--|:--|
| per-cell matrices | `ltqs_posterior`, `ltqs_posterior_error_bars` | `ltqs_measurement_approx`, `ltqs_measurement_approx_error_bars` | both pairs | both pairs + `gene_variances_likelihood` |
| always | `gene_means`, `gene_means_error_bars`, `gene_variances`, `geneID.txt`, `cellID.txt`, `sanity_command.txt` | the same | the same | the same + `gene_variances_grid.npy` with `-npy` |

The numerical files are `.txt`, `.txt.gz` with `-gz` or `.npy` with `-npy`.

## Bonsai output mode (optional)

`--bonsai-output-mode` is the preset `-npy -v_m EAP --ltq-estimate measured`: it writes the measurement estimates
as NPY, with the EAP gene variances that Bonsai needs, so a Bonsai run needs no other output option.

```bash
./bin/Sanity -f counts.tsv -d bonsai_output -n 12 --bonsai-output-mode
```

The folder then holds `ltqs_measurement_approx.npy` and `ltqs_measurement_approx_error_bars.npy` (float32,
`(genes, cells)`), `gene_means.npy`, `gene_means_error_bars.npy` and `gene_variances.npy` (float64, `(genes, 1)`),
`geneID.txt`, `cellID.txt` and `sanity_command.txt`, which records `# Output mode: bonsai` and
`# LTQ estimate: measured`. Bonsai reads the files by name and checks the two lines. Any other output option
combines with the preset as usual (`-e 1` adds the posterior files and the likelihood table); an explicit `-v_m`
other than `EAP` is an error, since Bonsai's gene variances are an inference choice, not an output choice.

Arrays use NPY version 1.0, little-endian values and C order. Row order matches
`geneID.txt`; matrix columns match `cellID.txt`. Inference and uncertainty
calculations remain double precision. In Python, arrays can be opened without loading the whole
matrix using `numpy.load(path, mmap_mode="r", allow_pickle=False)`.

## Usage
```
  ./Sanity <option(s)>
  Options:
    -h, --help                          Show this help message
    -v, --version                       Show the current version
    -f, --file                          Specify the input transcript count text file (.mtx for Matrix Market File Format)
    -mtx_genes, --mtx_gene_name_file    Specify the gene name text file (only needed if .mtx input file)
    -mtx_cells, --mtx_cell_name_file    Specify the cell name text file (only needed if .mtx input file)
    -d, --destination                   Specify the destination path (default: pwd)
    -n, --n_threads                     Specify the number of threads to be used (default: 4)
    -e, --extended_output               Option to print every output file (default: false, choice: false,0,true,1)
    -gz, --gzip_output                  Compress the numerical text output files with gzip (.gz); overridden by -npy
    -npy, --npy_output                  Write the numerical output files as NumPy arrays (.npy; float32 matrices, float64 per-gene vectors) instead of text; overrides -gz
    --bonsai-output-mode                The preset for Bonsai: -npy -v_m EAP --ltq-estimate measured
    -ltq, --ltq-estimate                Which per-cell LTQ estimates are written: posterior (x* with posterior error bars, default), measured (x^m with measurement error bars; the default in Bonsai output mode) or both; -e 1 writes both kinds
    -v_m, --v_method                    Option to specify the method for variance estimation (default: MAP, choice: MAP, EAP, MLE, MARG)
    -vmin, --variance_min               Minimal value of variance in log transcription quotient (default: 0.001)
    -vmax, --variance_max               Maximal value of variance in log transcription quotient (default: 50)
    -nbin, --number_of_variance_bins    Number of bins for the variance in log transcription quotient  (default: 160)
```

## Installation

Text formatting requires floating-point `std::to_chars` support in the C++17
standard library.

* Clone the GitHub repository
```
git clone https://github.com/jmbreda/Sanity.git
```
* Install OpenMP library
	* On Linux  
	If not already installed (Check with `ldconfig -p | grep libgomp`, no output if not installed), do
	```
	sudo apt-get update
	sudo apt-get install libgomp1
	```
	* On mac OS  
	Apple's own `g++`/`clang++` does not support OpenMP, so you need a real GCC.

		* Using brew (recommended)  
		```
		brew install gcc
		```
		The Makefile detects this automatically and picks the right `g++` binary, no further changes needed.

		* Using macports  
		```
		port install gcc13
		```
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Tell the Makefile to use it by compiling with `make CC=g++-mp-13` (adjust the version number to the package you installed).
* Move to the source code directory and compile.
```
cd Sanity/src
make
```
* The binary file is located in
```
Sanity/bin/Sanity
```
* Alternatively, pre-built binaries are published on the repository's [Releases](https://github.com/jmbreda/Sanity/releases) page.

## Sanity_distance
Compute cell-cell distances from Sanity output files. Needs the text output of `--ltq-estimate both` or `-e 1`: the posterior estimates for the signal-to-noise cut-off and the Euclidean distance, the measurement estimates (the prior divided out, SI Eqs. 62 and 63) for the distance with error bars.
### Input
* `-f`: The output folder of the Sanity run, specified with the `-d` option in Sanity (`'path/to/folder'`)
* `-s2n`: (optional) The gene signal to noise ratio used as gene cut-off (double, default: `1.0`)
* `-err`: (optional) Compute cell-cell distance with or without error bars (boolean, default: `1` or `true`)
* `-n`: (optional) Number of threads (integer, default: `4`)
### Output
* Cell-cell distance: *(N<sub>c</sub>(N<sub>c</sub>-1)/2)* vector of cell to cell distances *dist(cell<sub>i</sub>,cell<sub>j</sub>), i=1,...,N<sub>c</sub>-1, j=i+1,...,N<sub>c</sub>*, with *N<sub>c</sub>* the number of cells. 
  ||
  |:------:|
  |*dist(cell<sub>1</sub>,cell<sub>2</sub>)*|
  |*dist(cell<sub>1</sub>,cell<sub>3</sub>)*|
  |*dist(cell<sub>1</sub>,cell<sub>4</sub>)*|
  |...|
  |*dist(cell<sub>N<sub>c</sub>-2</sub>,cell<sub>N<sub>c</sub>-1</sub>)*|
  |*dist(cell<sub>N<sub>c</sub>-2</sub>,cell<sub>N<sub>c</sub></sub>)*|
  |*dist(cell<sub>N<sub>c</sub>-1</sub>,cell<sub>N<sub>c</sub></sub>)*|
located in the Sanity output folder (specified with `-f` option), named `cell_cell_distance_[...].txt`, depending on the `-err` and `-s2n` options.
### Usage
```
./Sanity_distance <option(s)>
Options:
	-h,--help		Show this help message
	-v,--version		Show the current version
	-f,--folder		Specify the input folder with extended output from Sanity
	-s2n,--signal_to_noise_cutoff	Minimal signal/noise of genes to include in the distance calculation (default: 1.0)
	-err,--with_error_bars	Compute cell-cell distance taking the errorbar epsilon into account (default: true)
	-n,--n_threads		Specify the number of threads to be used (default: 4)
```

### Installation
Same dependencies as Sanity (see above).

* Move to the source code directory and compile.
```
cd Sanity/src
make Sanity_distance
```
* The binary file is located in
```
Sanity/bin/Sanity_distance
```

<a name="sorting-mtx-file-by-row-gene-indices"></a>
### Sorting MTX file by row (gene) indices

Use script `sort_mtx_by_row.py` in the `Sanity/scripts` folder to sort the MTX file by row (gene) indices. The script takes as input the MTX file. The output is a sorted MTX file by default saved as `sorted_[original filename].gz` in the same folder as the input file. User can specify the output file name and path with the `-o` option. The saved MTX file is compressed with gzip.

Requirements: Python3, sort

Example usage:

```bash
python3 sort_mtx_by_row.py -i input.mtx
```

## Tests
A simple test comparing Sanity's output across the `-v_m` methods (MAP, EAP, MLE, MARG) on a small count matrix is available in the `tests` directory, together with an independent NumPy implementation of the fit's equations (`tests/compare_section2.py`). See [tests/README.md](tests/README.md) for details on running them.

## Help
For any questions or assistance regarding Sanity, please post your question in the issues section.
