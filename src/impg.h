#include <vector>
#include <map>

#include "util.h"
#include "linsubs.h"

const double LAMBDA = 0.1;

// load z-score typed snps file
size_t load_zscore_typed_snps(const char *filename,
	const vector<typed_snp>& typed_snps, vector<zscore_typed_snp> &ztyped_snps);

// load genotypes
size_t load_genotypes(const char *filename, 
	vector<genotype> &genotypes);

// load correlations from file
size_t load_cors(const char *filename, vector<double> &cors);

// load the ld matrix from file to memory
void load_ld_mat(const char *filename, size_t num_typed_snps, double *ld_mat);

// compute allele frequencies for all the snps
void get_all_freqs(const vector<string>& haps, vector<double> &freqs);

// compute the allele frequency of a snp, specified by snp_idx
double get_freq(const vector<string>& haps, size_t snp_idx);

// compute the h frequency between two snps, specified by snp_idx1, snp_idx2
double get_h_freq(const vector<string>& haps,size_t snp_idx1,size_t snp_idx2);

// get variance
double get_var(const vector<string>& haps,
		const vector<double>& freqs, int idx);

// get the command line input for gen_beta program
int get_gen_beta_cmd_line(int argc, char **argv, char **IN_HAP_FILE,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE,
		char **OUT_FILE_PREFIX);
		
// get command line input for impute z-scores
int get_imp_cmd_line(int argc, char **argv, char **PREFIX,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE, char **OUT_FILE);

// get command line input for gen beta with user defined cor
int get_gen_beta_cmd_line_cor(int argc, char **argv, char **IN_HAP_FILE,
                char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE,
                char **OUT_FILE_PREFIX, char **IN_COR_FILE);

// get the command line input for gen_beta with ld
int get_gen_beta_ld_cmd_line(int argc, char **argv, char **IN_HAP_FILE,
		char **IN_ALL_SNP_FILE, char **IN_TYPED_SNP_FILE, char **IN_LD_FILE,
		char **OUT_FILE_PREFIX);

// compute beta and var
void get_beta_var(const string& prefix, const vector<vector<double> >& haps,
	const vector<ref_snp>& all_snps, const vector<typed_snp>& typed_snps,
	const vector<double>& freqs, const vector<char>& impute_flags,
	double *sigma_t_tmp, double lambda);

// compute beta and var using self computed cor
void get_beta_var_cor(const string& prefix, const vector<ref_snp>& all_snps,
        const vector<typed_snp>& typed_snps, const vector<char>& impute_flags,
        const vector<char>& convert_flags, const vector<int>& typed_idx_map,
        double *sigma_t_tmp, double lambda, const vector<double> &cors, size_t sigma_dim);

void get_all_freqs_tmp(const vector<vector<double> >& haps, vector<double> &freqs);

double get_freq_tmp(const vector<vector<double> >& haps, size_t snp_idx);

// compute beta and var
void get_beta_var_raw(const string& prefix, const vector<string>& haps,
	const vector<genotype>& all_snps, const vector<genotype>& typed_snps,
	const vector<double>& freqs, double *sigma_t_tmp, double lambda);
