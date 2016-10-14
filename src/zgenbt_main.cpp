#include <iostream>
#include <cmath>
#include <stdio.h>

#include "impg.h"
#include "linsubs.h"

using namespace std;

int main(int argc, char **argv) {
	string filename, prefix;
	
    // get command line input
	char *IN_HAP_FILE = NULL, *IN_ALL_SNP_FILE = NULL;
	char *IN_TYPED_SNP_FILE = NULL, *OUT_FILE_PREFIX = NULL;
	get_gen_beta_cmd_line(argc, argv, &IN_HAP_FILE, &IN_ALL_SNP_FILE,
		&IN_TYPED_SNP_FILE, &OUT_FILE_PREFIX);
	prefix = string(OUT_FILE_PREFIX);

	// load typed snps
	vector<typed_snp> typed_snps;
	size_t num_typed_snps = load_typed_snps(IN_TYPED_SNP_FILE, typed_snps);
	
	// load all snps
	vector<ref_snp> all_snps;
	size_t num_total_snps = load_all_snps(IN_ALL_SNP_FILE, all_snps);

	// mark snps whose z-scores need to be converted in the imputation step
	vector<char> convert_flags;
	convert_flags.resize(num_typed_snps, 0);
	vector<char> impute_flags;
	impute_flags.resize(num_total_snps, 1);
	mark_snps(typed_snps, all_snps, convert_flags, impute_flags);

	// load haplotypes
	vector<string> haps;
    vector<vector<double> > haps_float;
	//load_haplotypes(IN_HAP_FILE, haps, num_total_snps);
	load_haplotypes_tmp(IN_HAP_FILE, haps_float, num_total_snps);

	// compute allele frequencies for all snps
	vector<double> freqs;
	freqs.resize(num_total_snps, 0.0);
	get_all_freqs_tmp(haps_float, freqs);
	
	// estimate sigma_t matrix using maf filtered typed snps
	size_t sigma_t_tmp_size = num_typed_snps*num_typed_snps;
	double *sigma_t_tmp = (double*)safe_calloc(sigma_t_tmp_size,
		sizeof(double));
	
	for(size_t i = 0; i < num_typed_snps; i++) {
		size_t idxi = typed_snps[i].idx;
		for(size_t j = i+1; j < num_typed_snps; j++) {
			size_t idxj = typed_snps[j].idx;
			double r = corr(haps_float[idxi], haps_float[idxj]);  
			sigma_t_tmp[i*num_typed_snps+j] = r;
			sigma_t_tmp[j*num_typed_snps+i] = sigma_t_tmp[i*num_typed_snps+j];
		}
		sigma_t_tmp[i*num_typed_snps+i] = 1.0;
	}

	// output the betas and vars
	get_beta_var(prefix, haps_float, all_snps, typed_snps,
		freqs, impute_flags, sigma_t_tmp, LAMBDA);

	// clean up
	free(sigma_t_tmp);	
	free(IN_HAP_FILE);
	free(IN_ALL_SNP_FILE);
	free(IN_TYPED_SNP_FILE);
	free(OUT_FILE_PREFIX);
}
