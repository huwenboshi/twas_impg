#include <iostream>
#include <cmath>
#include <stdio.h>

#include "impg.h"
#include "linsubs.h"

using namespace std;

int main(int argc, char **argv) {
	string filename, prefix;
	
    // get command line input
	char *IN_LD_FILE = NULL, *IN_ALL_SNP_FILE = NULL;
	char *IN_TYPED_SNP_FILE = NULL, *OUT_FILE_PREFIX = NULL;
	char *IN_COR_FILE = NULL;
    get_gen_beta_cmd_line_cor(argc, argv, &IN_LD_FILE, &IN_ALL_SNP_FILE,
		&IN_TYPED_SNP_FILE, &OUT_FILE_PREFIX, &IN_COR_FILE);
	prefix = string(OUT_FILE_PREFIX);

    // load cor file
    vector<double> cors;
    load_cors(IN_COR_FILE, cors);

	// load typed snps
	vector<typed_snp> typed_snps;
	size_t num_typed_snps = load_typed_snps(IN_TYPED_SNP_FILE, typed_snps);
	
	// load all snps
	vector<ref_snp> all_snps;
	size_t num_total_snps = load_all_snps(IN_ALL_SNP_FILE, all_snps);

    // append gene expression
    ref_snp rs;
    rs.idx = num_total_snps;
    rs.snp_name = "gene_exp";
    rs.snp_pos = -1;
    rs.ref_allele = 'X';
    rs.alt_allele = 'Y'; 
    all_snps.push_back(rs);

	// mark snps whose z-scores need to be converted in the imputation step
	vector<char> convert_flags;
	convert_flags.resize(num_typed_snps, 0);
	vector<char> impute_flags;
	vector<int> typed_idx_map;
    impute_flags.resize(num_total_snps+1, 1);
	intersect_snps(typed_snps, all_snps, convert_flags, impute_flags, typed_idx_map);

    int num_inter = 0;
    for(size_t i=0; i<convert_flags.size(); i++) {
        if(convert_flags[i] != 0) {
            num_inter++;
        }
    }
    if(num_inter == 0) {
        fprintf(stderr, "Error: Empty intersection\n");
        exit(1);
    }

	// load LD
    size_t sigma_size = num_total_snps*num_total_snps;
	double *sigma = (double*)safe_calloc(sigma_size, sizeof(double));
	load_ld_mat(IN_LD_FILE, num_total_snps, sigma);

	// output the betas and vars
	get_beta_var_cor(prefix, all_snps, typed_snps, impute_flags,
        convert_flags, typed_idx_map, sigma, LAMBDA, cors, num_total_snps);

	// clean up
	free(sigma);	
	free(IN_LD_FILE);
	free(IN_ALL_SNP_FILE);
	free(IN_TYPED_SNP_FILE);
	free(OUT_FILE_PREFIX);
}
