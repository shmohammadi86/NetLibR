#ifndef GRAPHLIBR_H
#define GRAPHLIBR_H

#include <omp.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>

#include <steiner/pcst_fast.h>

#include <bmatching_seq/bMatching.h>
#include <bmatching_seq/mtxReader.h>
#include <bmatching_seq/pcl_stack.h>

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

#define AND_SYM 0
#define OR_SYM 1

#define COR_DIST 0
#define ACTIONet_DIST 1
#define JSD_DIST 2
#define NEG_LOG 3


namespace NetLibR {
	// Graph matching algorithms
	mat MWM_hungarian(mat &G);
	umat MWM_rank1(vec u, vec v, double u_threshold, double v_threshold);

	// [Approx] Nearest-neighbor graph construction

	// K*-NN
	field<sp_mat> buildAdaptiveNet_fromDist(mat &D, int max_kNN, double LC, int sym_method, int thread_no, int dist_type, double sigma); // k*-NN
	field<sp_mat> buildAdaptiveNet_fromEdgeList(mat &dist, mat &idx, int max_kNN, double LC, int sym_method, int thread_no, int dist_type, double sigma); // k*-NN

	// PCSF
	sp_mat PCSF(sp_mat &Adj, vec node_scores, int root, int clusters, int convert_similarity_to_distance);

	// Graph spanner
	sp_mat make_spanner(sp_mat &G_adj, int k);

	// Graph alignment
	sp_mat netAlign_arma(sp_mat A_mat, sp_mat B_mat, sp_mat L_mat, double alpha, double beta, double gamma, int maxiter, bool finalize);
	
	// Graph clustering
	vec unsigned_cluster(sp_mat A, double resolution_parameter, uvec initial_clusters, int seed);
	vec signed_cluster(sp_mat A, double resolution_parameter, uvec initial_clusters, int seed);


	// Graph diffusion
	mat batchPR(sp_mat &G, mat &U, double alpha, int thread_no, double tol);
	mat batch_zoned_diffusion(sp_mat &G, uvec &zones, mat &U, double alpha, int thread_no, double tol);
	
	vec sweepcut(sp_mat &A, vec s);

	// Graph autocorrelation
	field<vec> computeAutocorrelation_Geary(sp_mat &G, mat &scores, int rand_perm_no, int num_shuffles);

	// Graph embedding
}

#endif
