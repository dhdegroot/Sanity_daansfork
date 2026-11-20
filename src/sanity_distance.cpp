#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <omp.h>

#include "Version.h"

using namespace std;

// "SI" in the comments is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021 (Section S1.3).
// Index letters at the end of a name give what an array runs over: g = genes, c = cells, a = grid points of alpha.
// Two letters mean a matrix stored in that order, so log_fold_change_cg is indexed [cell][gene].

// Distances functions :
// with error bars:
void get_distance_errorbar(double **log_fold_change_cg, double **log_fold_change_error_sq_cg, double *gene_variance_g, double **distance_cc, int num_genes, int num_cells, int num_threads);
void distance_row_errorbar(double **unshrunk_log_fold_change_cg, double **unshrunk_log_fold_change_error_sq_cg, double **alpha_times_variance_ga, double *distance_row_c, int num_genes, int num_cells, int i, int num_alpha_bins);
// without errorbars:
void get_distance_euclidean(double **log_fold_change_cg, double **distance_cc, int num_genes, int num_cells, int num_threads);
void distance_row_euclidean(double **log_fold_change_cg, double *distance_row_c, int num_genes, int num_cells, int i);

// I/O
void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, bool &with_error_bar, int &num_threads, int &num_genes, int &num_cells);
static void show_usage(void);
void load_matrix(string sanity_folder, double **log_fold_change_cg, double **log_fold_change_error_sq_cg, double *gene_variance_g, double **measurement_log_fold_change_cg, double **measurement_error_sq_cg, int num_genes, int num_cells, bool with_error_bar);
void print_distance(string out_file, double **distance_cc, int num_cells);

int main(int argc, char** argv){


	// Declare arguments :
	string sanity_folder("");
	double s2n_cutoff = 1.0;
	bool with_error_bar(true);
	int num_threads(4);
	int num_genes(0);
	int num_cells(0);
	parse_argv(argc, argv, sanity_folder, s2n_cutoff, with_error_bar, num_threads, num_genes, num_cells);

	cout << "Sanity folder :\n\t" << sanity_folder << "\n";
	cout << "Signal/noise cut-off :" << s2n_cutoff << "\n";
	if (with_error_bar)
		cout << "With error bars";
	else
		cout << "Without error bars";
	cout << "\n";
	cout << "Nr. of gene :" << num_genes << "\n";
	cout << "Nr. of cell :" << num_cells << "\n";
	cout << "Using " << num_threads << " cpus\n";


	// Declare temporary variables
	int g, c;
	double **log_fold_change_all_cg = new double *[num_cells]; // <delta_gc> (SI Eq. 39), indexed [cell][gene], from ltqs_posterior.txt minus gene_means.txt; before the signal-to-noise filter.
	for(c=0;c<num_cells;c++)
		log_fold_change_all_cg[c] = new double[num_genes];

	double **log_fold_change_error_sq_all_cg = new double *[num_cells]; // epsilon_gc^2 (SI Eq. 42), indexed [cell][gene], from the error bars in ltqs_posterior_error_bars.txt.
	for(c=0;c<num_cells;c++)
		log_fold_change_error_sq_all_cg[c] = new double[num_genes];
	
	double *gene_variance_all_g = new double [num_genes]; // v_g (SI Eq. 17), read from gene_variances.txt.

	// delta~_gc and eta_gc^2 (SI Eqs. 62, 63): the measurement LTQs of Sanity minus the gene mean, with their error
	// bars squared, from ltqs_measurement_approx.txt and ltqs_measurement_approx_error_bars.txt; only with error bars.
	double **measurement_log_fold_change_all_cg = NULL;
	double **measurement_error_sq_all_cg = NULL;
	if(with_error_bar){
		measurement_log_fold_change_all_cg = new double *[num_cells];
		measurement_error_sq_all_cg = new double *[num_cells];
		for(c=0;c<num_cells;c++){
			measurement_log_fold_change_all_cg[c] = new double[num_genes];
			measurement_error_sq_all_cg[c] = new double[num_genes];
		}
	}

	// Get data
	cout << "Get data... ";
	load_matrix(sanity_folder, log_fold_change_all_cg, log_fold_change_error_sq_all_cg, gene_variance_all_g, measurement_log_fold_change_all_cg, measurement_error_sq_all_cg, num_genes, num_cells, with_error_bar);	

	// Compute signal 2 noise
	bool keep_gene_g[num_genes]; // True for the genes whose signal-to-noise reaches the cut-off.
	int num_genes_kept(0);
	if (s2n_cutoff > 0.0){
	
		double mean_log_fold_change; // Mean of <delta_gc> over the cells.
		double log_fold_change_sample_var; // Sample variance of <delta_gc> over the cells.
		double mean_log_fold_change_error_sq; // Mean of epsilon_gc^2 over the cells.
		for(g=0;g<num_genes;g++){
			
			mean_log_fold_change = 0.0;
			for(c=0;c<num_cells;c++){
				mean_log_fold_change += log_fold_change_all_cg[c][g];
			}
			mean_log_fold_change /= (double)num_cells;

			log_fold_change_sample_var = 0.0;
			for(c=0;c<num_cells;c++){
				log_fold_change_sample_var += (log_fold_change_all_cg[c][g] - mean_log_fold_change)*(log_fold_change_all_cg[c][g] - mean_log_fold_change);
			}
			log_fold_change_sample_var /= (double)num_cells - 1.0 ;
	
			mean_log_fold_change_error_sq = 0.0;
			for(c=0;c<num_cells;c++){
				mean_log_fold_change_error_sq += log_fold_change_error_sq_all_cg[c][g];
			}
			mean_log_fold_change_error_sq /= (double)num_cells;

			if(log_fold_change_sample_var/mean_log_fold_change_error_sq >= s2n_cutoff){
				keep_gene_g[g] = true;
				num_genes_kept += 1;
			}
			else{
				keep_gene_g[g] = false;
			}
		}
	}else{
		num_genes_kept = num_genes;
	}
	
	// The values the distance is computed from, for the genes above the cut-off, indexed [cell][gene]: with error
	// bars delta~_gc and eta_gc^2 (SI Eqs. 62, 63), without them <delta_gc> alone.
	double **distance_input_all_cg = with_error_bar ? measurement_log_fold_change_all_cg : log_fold_change_all_cg;
	double **distance_input_error_sq_all_cg = with_error_bar ? measurement_error_sq_all_cg : log_fold_change_error_sq_all_cg;
	double **log_fold_change_cg = new double *[num_cells];
	for(c=0;c<num_cells;c++)
		log_fold_change_cg[c] = new double[num_genes_kept];

	double **log_fold_change_error_sq_cg = new double *[num_cells];
	for(c=0;c<num_cells;c++)
		log_fold_change_error_sq_cg[c] = new double[num_genes_kept];

	double *gene_variance_g = new double [num_genes_kept]; // v_g of the genes above the cut-off.

	if (s2n_cutoff > 0.0){
		// Copy genes with signal to noise above cut-off
		int k = 0;
		for(g=0;g<num_genes;g++){
			if(keep_gene_g[g]){
				for(c=0;c<num_cells;c++){
					log_fold_change_cg[c][k] = distance_input_all_cg[c][g];
					log_fold_change_error_sq_cg[c][k] = distance_input_error_sq_all_cg[c][g];
				}
				gene_variance_g[k] = gene_variance_all_g[g];
				k++;
			}
		}
	}else{
		for(g=0;g<num_genes;g++){
            for(c=0;c<num_cells;c++){
                log_fold_change_cg[c][g] = distance_input_all_cg[c][g];
                log_fold_change_error_sq_cg[c][g] = distance_input_error_sq_all_cg[c][g];
            }
            gene_variance_g[g] = gene_variance_all_g[g];
        }
	}

	num_genes = num_genes_kept;
	cout << "N_gene after cut-off :" << num_genes << "\n"; 

	// delete temporary variables
	for(c=0;c<num_cells;c++){
		delete[] log_fold_change_all_cg[c];
		delete[] log_fold_change_error_sq_all_cg[c];
		if(with_error_bar){
			delete[] measurement_log_fold_change_all_cg[c];
			delete[] measurement_error_sq_all_cg[c];
		}
	}
	delete[] log_fold_change_all_cg;
	delete[] log_fold_change_error_sq_all_cg;
	delete[] measurement_log_fold_change_all_cg;
	delete[] measurement_error_sq_all_cg;
	delete[] gene_variance_all_g;

	// Get distance
	cout << "Get distance... ";
	double **distance_cc = new double *[num_cells]; // d_cc', the distance between every pair of cells (SI Eqs. 50, 72).
	for(c=0; c<num_cells; c++)
		distance_cc[c] = new double[num_cells];

	int i,j;
	for(i=0;i<num_cells;i++){
		for(j=0;j<num_cells;j++){
			distance_cc[i][j] = 0.0;
		}
	}
	if(with_error_bar)
		get_distance_errorbar(log_fold_change_cg,log_fold_change_error_sq_cg,gene_variance_g,distance_cc,num_genes,num_cells,num_threads);
	else
		get_distance_euclidean(log_fold_change_cg,distance_cc,num_genes,num_cells,num_threads);
  	cout << "\n";

	// Print distance
	string out_file = sanity_folder + "cell_cell_distance";
	if(with_error_bar)
		out_file += "_with_errorbar";
	else
		out_file += "_euclidean";
	if(s2n_cutoff > 0.0){
		string s2n_str = to_string(s2n_cutoff);
		// Strip trailing zeros from the fractional part only, then the '.' itself.
		// Stopping at the '.' is what keeps 10 -> "10" rather than "1".
		if(s2n_str.find('.') != string::npos){
			while(s2n_str.back() == '0')
				s2n_str = s2n_str.substr(0, s2n_str.size()-1);
			if(s2n_str.back() == '.')
				s2n_str = s2n_str.substr(0, s2n_str.size()-1);
		}
		out_file += "_s2n_gt_" + s2n_str;
	}
	out_file += ".txt";

	cout << "Print distance in " << out_file;
	print_distance(out_file, distance_cc, num_cells);
  	cout << "\n";

	return 0;
}

void get_distance_errorbar(double **log_fold_change_cg, double **log_fold_change_error_sq_cg, double *gene_variance_g, double **distance_cc, int num_genes, int num_cells, int num_threads){

	// log_fold_change_cg and log_fold_change_error_sq_cg hold delta~_gc and eta_gc^2 (SI Eqs. 62, 63), the
	// measurement LTQs of Sanity minus the gene mean: the prior N(0, v_g) is already divided out.
	int g;

	// alpha_times_variance_ga = alpha*v_g for every gene and grid point of alpha (SI Eq. 56)
	int a;
	int num_alpha_bins(401); // Number of grid points of the hyperparameter alpha (SI Eq. 56).
	double alpha_step(0.005); // Spacing of the grid of alpha, which runs from 0 to 2.
	double alpha; // Current value of alpha: 0 means identical cells and 2 means as distant as a random pair (SI Eq. 56).

	double **alpha_times_variance_ga = new double *[num_genes];
	for(g=0;g<num_genes;g++)
		alpha_times_variance_ga[g] = new double [num_alpha_bins];

	for(g=0;g<num_genes;g++){
		alpha = 0.0;
		for(a=0;a<num_alpha_bins;a++){
			alpha_times_variance_ga[g][a] = alpha*gene_variance_g[g];
			alpha += alpha_step;
		}
	}

	int i, num_row_pairs;
	if ( (num_cells % 2) == 0 ){
		num_row_pairs = num_cells/2;
	}else{
		num_row_pairs =(num_cells-1)/2;
		i = num_row_pairs;
		distance_row_errorbar(log_fold_change_cg,log_fold_change_error_sq_cg,alpha_times_variance_ga,distance_cc[i],num_genes,num_cells,i,num_alpha_bins);
	}
	#pragma omp parallel for num_threads(num_threads)
	for (i=0; i<num_row_pairs; i++){
		distance_row_errorbar(log_fold_change_cg,log_fold_change_error_sq_cg,alpha_times_variance_ga,distance_cc[i],num_genes,num_cells,i,num_alpha_bins);
		distance_row_errorbar(log_fold_change_cg,log_fold_change_error_sq_cg,alpha_times_variance_ga,distance_cc[num_cells-1-i],num_genes,num_cells,num_cells-i-1,num_alpha_bins);
	}
	return;
}

// Fills distance_row_c[j] = d_ij for all cells j > i.
// unshrunk_log_fold_change_cg is delta~_gc and unshrunk_log_fold_change_error_sq_cg is eta_gc^2 (SI Eqs. 62, 63), after the rescaling in get_distance_errorbar.
void distance_row_errorbar(double **unshrunk_log_fold_change_cg, double **unshrunk_log_fold_change_error_sq_cg, double **alpha_times_variance_ga, double *distance_row_c, int num_genes, int num_cells, int i, int num_alpha_bins){
	
	int j, g, a;
	double pair_total_variance, shrinkage; // pair_total_variance = alpha v_g + eta_g^2, the variance of x_g (SI Eq. 67); shrinkage = f_g(alpha) (SI Eq. 71).
	double log_fold_change_diff_sq_g[num_genes]; // x_g^2, the squared difference of the estimated log fold changes of the two cells (SI Eq. 52).
	double pair_error_sq_g[num_genes]; // eta_g^2 = eta_gc^2 + eta_gc'^2, the noise on x_g (SI text after Eq. 63).
	double alpha_likelihood_a[num_alpha_bins]; // First log P(D_cc' | alpha, v) (SI Eq. 67), then the posterior P(alpha | D_cc', v) (SI Eq. 68).
	double max_log_likelihood, likelihood_sum; // Maximum of the log-likelihood over alpha, and the sum used to normalise it.
	for (j=i+1; j<num_cells; j++){
		
		// Compute
		// x_g^2 = (delta_gi - delta_gj)^2
		// pair_error_sq_g = eta_gc^2 + eta_gc'^2
		// shrinkage = f_g(alpha) = alpha_times_variance_ga/(pair_error_sq_g + alpha_times_variance_ga) (SI Eq. 71)
		for(g=0;g<num_genes;g++){
			log_fold_change_diff_sq_g[g] = (unshrunk_log_fold_change_cg[i][g] - unshrunk_log_fold_change_cg[j][g])*(unshrunk_log_fold_change_cg[i][g] - unshrunk_log_fold_change_cg[j][g]);
			pair_error_sq_g[g] = unshrunk_log_fold_change_error_sq_cg[i][g] + unshrunk_log_fold_change_error_sq_cg[j][g];
		}

		// get likelihood
		max_log_likelihood = -1000000000.0;
		for(a=0;a<num_alpha_bins;a++){
			// log-likelihood over all genes
			alpha_likelihood_a[a] = 0.0;
			for(g=0;g<num_genes;g++){
				pair_total_variance= pair_error_sq_g[g] + alpha_times_variance_ga[g][a];
				alpha_likelihood_a[a] -= .5*log_fold_change_diff_sq_g[g]/pair_total_variance;
				alpha_likelihood_a[a] -= .5*log(pair_total_variance);
			}
			// keep  track of the maximum log-likelihood
			if(alpha_likelihood_a[a] > max_log_likelihood)
				max_log_likelihood = alpha_likelihood_a[a];
		}

		// normalize
		likelihood_sum = 0.0;
		for(a=0;a<num_alpha_bins;a++){
			alpha_likelihood_a[a] -= max_log_likelihood;
			// turn into likelihood
			alpha_likelihood_a[a] = exp(alpha_likelihood_a[a]);
			// keep track of normalization
			likelihood_sum += alpha_likelihood_a[a];
		}

		// calculate expectation values
		distance_row_c[j] = 0.0;
		for(a=0;a<num_alpha_bins;a++){
			// posterior of this alpha
			alpha_likelihood_a[a] /= likelihood_sum;
			for(g=0;g<num_genes;g++){
				// This is f_g(alpha), SI Eq. 71
				shrinkage = alpha_times_variance_ga[g][a]/(alpha_times_variance_ga[g][a] + pair_error_sq_g[g]);
				// this is the contribution of gene g to the squared distance (SI Eq. 72)
				distance_row_c[j] += alpha_likelihood_a[a]*(shrinkage*shrinkage*log_fold_change_diff_sq_g[g] + shrinkage*pair_error_sq_g[g]);
			}
		}
		distance_row_c[j] = sqrt(distance_row_c[j]);
	}
	return;
}


void get_distance_euclidean(double **log_fold_change_cg, double **distance_cc, int num_genes, int num_cells, int num_threads){

    int i, num_row_pairs;
    if ( (num_cells % 2) == 0 ){
        cout << "N_cell even\n";
        num_row_pairs = num_cells/2;
    }else{
        cout << "N_cell odd\n";
        num_row_pairs =(num_cells-1)/2;
        i = num_row_pairs;
        distance_row_euclidean(log_fold_change_cg,distance_cc[i],num_genes,num_cells,i);
    }
    #pragma omp parallel for num_threads(num_threads)
    for (i=0; i<num_row_pairs; i++){
        distance_row_euclidean(log_fold_change_cg,distance_cc[i],num_genes,num_cells,i);
        distance_row_euclidean(log_fold_change_cg,distance_cc[num_cells-1-i],num_genes,num_cells,num_cells-i-1);
    }
    return;
}

void distance_row_euclidean(double **log_fold_change_cg, double *distance_row_c, int num_genes, int num_cells, int i){

    int j, g;
    double log_fold_change_diff; // Difference of the log fold changes of the two cells for one gene.
    for (j=i+1; j<num_cells; j++){
        distance_row_c[j] = 0.0;
        for (g=0; g<num_genes; g++){
            log_fold_change_diff = log_fold_change_cg[i][g] - log_fold_change_cg[j][g];
            distance_row_c[j] += log_fold_change_diff*log_fold_change_diff;
        }
        distance_row_c[j] = sqrt(distance_row_c[j]);
    }

    return;
}



// Reads a per-gene table of Sanity 3.0: a header line of cell names, then one line per gene that starts with the
// gene name. store(g, c, value) receives every number. Exits on a short or unreadable line.
template <class Store>
static void read_matrix_table(const string &path, int num_genes, int num_cells, Store store){
	ifstream in(path);
	if(!in){
		cerr << "Cannot open input file " << path << "\n";
		exit(EXIT_FAILURE);
	}
	string line;
	getline(in, line); // The header of cell names.
	for(int g=0; g<num_genes; g++){
		if(!getline(in, line)){
			cerr << "Error: Couldn't read a line at row number " << g << " of " << path << "\n";
			exit(EXIT_FAILURE);
		}
		const char *pos = line.c_str() + strcspn(line.c_str(), " \t"); // Past the gene name.
		char *end;
		for(int c=0; c<num_cells; c++){
			const double value = strtod(pos, &end);
			if(end == pos){
				cerr << "Error: not enough fields on line number " << g << " of " << path << "\n";
				exit(EXIT_FAILURE);
			}
			store(g, c, value);
			pos = end;
		}
		strtod(pos, &end);
		if(end != pos)
			cerr << "Error: too many fields on line number " << g << " of " << path << "\n";
	}
}

// Reads a per-gene vector of Sanity: one number per line.
static void read_vector_file(const string &path, int num_genes, double *values_g){
	ifstream in(path);
	if(!in){
		cerr << "Cannot open input file " << path << "\n";
		exit(EXIT_FAILURE);
	}
	string line;
	for(int g=0; g<num_genes; g++){
		if(!getline(in, line)){
			cerr << "Error: Couldn't read a line at row number " << g << " of " << path << "\n";
			exit(EXIT_FAILURE);
		}
		char *end;
		values_g[g] = strtod(line.c_str(), &end);
		if(end == line.c_str()){
			cerr << "Error: not enough fields on line number " << g << " of " << path << "\n";
			exit(EXIT_FAILURE);
		}
	}
}

// Reads the extended output of Sanity (-e 1): the posterior LTQs x*_gc minus the gene mean mu_g, their error bars
// squared, and with error bars also the gene variances v_g and the measurement LTQs x^m_gc minus mu_g with their
// error bars squared. The measurement LTQs are the posterior with the prior divided out, which this tool used to do
// itself (SI Eqs. 62, 63); Sanity now does it exactly.
void load_matrix(string sanity_folder, double **log_fold_change_cg, double **log_fold_change_error_sq_cg, double *gene_variance_g, double **measurement_log_fold_change_cg, double **measurement_error_sq_cg, int num_genes, int num_cells, bool with_error_bar){

	double *gene_mean_g = new double[num_genes]; // mu_g, read from gene_means.txt.
	read_vector_file(sanity_folder + "gene_means.txt", num_genes, gene_mean_g);
	read_matrix_table(sanity_folder + "ltqs_posterior.txt", num_genes, num_cells,
		[&](int g, int c, double ltq){ log_fold_change_cg[c][g] = ltq - gene_mean_g[g]; });
	read_matrix_table(sanity_folder + "ltqs_posterior_error_bars.txt", num_genes, num_cells,
		[&](int g, int c, double error){ log_fold_change_error_sq_cg[c][g] = error*error; });
	if(with_error_bar){
		read_vector_file(sanity_folder + "gene_variances.txt", num_genes, gene_variance_g);
		read_matrix_table(sanity_folder + "ltqs_measurement_approx.txt", num_genes, num_cells,
			[&](int g, int c, double ltq){ measurement_log_fold_change_cg[c][g] = ltq - gene_mean_g[g]; });
		read_matrix_table(sanity_folder + "ltqs_measurement_approx_error_bars.txt", num_genes, num_cells,
			[&](int g, int c, double error){ measurement_error_sq_cg[c][g] = error*error; });
	}
	delete[] gene_mean_g;
}


void print_distance(string out_file, double **distance_cc, int num_cells){
	
	int i, j;
	string file_name(out_file);
	ofstream outfile(file_name);
	for (i=0; i<(num_cells-1); i++){
		for (j=i+1; j<num_cells; j++)
			outfile << distance_cc[i][j] << "\n";
	}
	outfile.close();
}


void parse_argv(int argc,char** argv, string &sanity_folder, double &s2n_cutoff, bool &with_error_bar, int &num_threads, int &num_genes, int &num_cells){

    if (argc<2)
        show_usage();

    int i;

    string get_help [2];
    get_help[0] = "-h";
    get_help[1] = "--help";

    for(i=1;i<argc;i++){
        if (argv[i] == get_help[0] || argv[i] == get_help[1])
            show_usage();
    }

    string get_version [2];
    get_version[0] = "-v";
    get_version[1] = "--version";
    for(i=1;i<argc;i++){
        if (argv[i] == get_version[0] || argv[i] == get_version[1]){
            cout << "Sanity_distance version " << SANITY_VERSION << "\n";
            exit(0);
        }
    }

    int N_param(4);
    string error_bar("true");
    string to_find[4][2] = {{"-f", "--folder"},
                            {"-s2n", "--signal_to_noise_cutoff"},
                            {"-err","--with_error_bars"},
                            {"-n", "--n_threads"}};

    int j;
    int idx;
    for(j=0;j<N_param;j++){
        idx = 0;
        for(i=1;i<argc;i++){
            if (argv[i] == to_find[j][0] || argv[i] == to_find[j][1]){
                idx = i;
                if ( idx+1 > argc-1 ){
                    cerr << "Error in argument parsing :\n"
                         << argv[i] << " option missing\n";
                    show_usage();
                }
                if(j==0) sanity_folder = argv[idx+1];
                if(j==1) s2n_cutoff = atof(argv[idx+1]);
				if(j==2) error_bar = argv[idx+1];
                if(j==3) num_threads = atoi(argv[idx+1]);

                // add '/' to out_folder if not already. `-f ""` gives an empty string, on
                // which back() would be undefined behaviour, so treat it as the cwd.
                if( j == 0 ){
                    if( sanity_folder.empty() )
                        sanity_folder = "./";
                    else if( sanity_folder.back() != '/' )
                        sanity_folder = sanity_folder + '/';
                }
            }
        }
        if (idx == 0 && j == 0){
            cerr << "Error in argument parsing :\n"
            << "missing input file name\n";
            show_usage();
        }
    }
    
	if ( error_bar == "false" || error_bar == "0" )
        with_error_bar = false;

    // Get number of gene and number of cells
	num_genes = -1;
	num_cells = -1;
	ifstream in(sanity_folder + "ltqs_posterior.txt");
	if (in.is_open()){
		num_genes = 0;
		num_cells = 0;
		string line;
		
		// The header line holds "GeneID" and then the cell names; every later line is one gene.
		getline(in, line);
		istringstream iss(line);
		do{
			string sub;
			iss >> sub;
			if (sub.length())
				++num_cells;
		}
		while(iss);
		--num_cells;

		while ( getline(in, line) ){
			   ++num_genes;
		}
	}else{
		cerr << "Unable to open " + sanity_folder + "ltqs_posterior.txt\n";
	}

}

static void show_usage(void)
{
    cerr << "Usage: Sanity_distance <option(s)> SOURCES\n"
         << "Options:\n"
         << "\t-h,--help\t\tShow this help message\n"
         << "\t-v,--version\t\tShow the current version\n"
         << "\t-f,--folder\t\tSpecify the input folder with extended output from Sanity\n"
         << "\t-s2n,--signal_to_noise_cutoff\tMinimal signal/noise of genes to include in the distance calculation (default: 1.0)\n"
         << "\t-err,--with_error_bars\tCompute cell-cell distance taking the errobar epsilon into account (default: 1)\n"
         << "\t-n,--n_threads\t\tSpecify the number of threads to be used (default: 4)\n";
    exit(0);
}



