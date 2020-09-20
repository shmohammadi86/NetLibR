#define ARMA_USE_CXX11_RNG
#define ARMA_64BIT_WORD

#include <NetLibR.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat MWM_hungarian(mat &G) {	
	mat G_matched = NetLibR::MWM_hungarian(G);

    return G_matched;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
umat MWM_rank1(vec u, vec v, double u_threshold = 0, double v_threshold = 0) {
	
	umat pairs = NetLibR::MWM_rank1(u, v, u_threshold, v_threshold);
	
	pairs = pairs + 1;
	
	return(pairs);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec signed_cluster(sp_mat A, double resolution_parameter = 1.0, int seed = 0, Nullable<IntegerVector> initial_clusters_ = R_NilValue) {
    set_seed(seed);

	uvec initial_clusters_uvec(A.n_rows);
	if ( initial_clusters_.isNotNull() ) {
        NumericVector initial_clusters(initial_clusters_);
		
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = initial_clusters(i);
	} else {
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
	}

	
	vec clusters = NetLibR::signed_cluster(A, resolution_parameter, initial_clusters_uvec, seed);

    return clusters;	
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec unsigned_cluster(sp_mat A, double resolution_parameter = 1.0, int seed = 0, Nullable<IntegerVector> initial_clusters_ = R_NilValue) {
    set_seed(seed);


	uvec initial_clusters_uvec(A.n_rows);
	if ( initial_clusters_.isNotNull() ) {
        NumericVector initial_clusters(initial_clusters_);
		
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = initial_clusters(i);
	} else {
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
	}

	vec clusters = NetLibR::unsigned_cluster(A, resolution_parameter, initial_clusters_uvec, seed);

    return clusters;	
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat bmatching(sp_mat A, vec b_limit, int seed = 0) {
    set_seed(seed);
	register int i, j, k;
	
    CSR G;
	G.readSpMat(A);
	
    int* nlocks=(int*)malloc(G.nVer*sizeof(int));    //required for all schemes
    int* start=(int*)malloc(G.nVer*sizeof(int));     //required for sorted and part. sorted
    int* end=(int*)malloc(G.nVer*sizeof(int));       //required for part sorted

    int *b=(int*)malloc(G.nVer*sizeof(int));
    for(i = 0; i < (unsigned int)G.nVer; i++) {
		int card = 0;
		for(j = G.verPtr[i]; j < (unsigned int)G.verPtr[i+1]; j++)
			if(0 < G.verInd[j].weight)
				card++;

		if(b_limit[i] < card)
			b[i]=b_limit[i];
		else b[i]=max(1, card);
    }

    Node* S=(Node*)malloc(G.nVer*sizeof(Node));      //Heap data structure
	for(i = 0; i < (unsigned int)G.nVer; i++) {
		S[i].heap=(Info*)malloc(max(1, b[i])*sizeof(Info));      //Each heap of size b
		S[i].curSize = 0;
	}
	
	localDom(&G, S, b, start, end, nlocks, 3, 2, 1);	
	
	
	sp_mat A_matched(size(A));
	
	/*
	for (i = 0; i < A.n_rows; i++) {
		for(int k = 0; k < (unsigned int) S[i].curSize; k++) {
			j = (S[i].heap[k].id) - A.n_rows;
			A_matched(i, j) = A(i, j);
		}
	}	
	*/
	return(A_matched);		
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat make_spanner(sp_mat &G_adj, int k) {
	sp_mat G_sparse = NetLibR::make_spanner(G_adj, k);

    return G_sparse;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat netAlign(sp_mat A, sp_mat B, sp_mat L,
							double alpha = 1.0,
							double beta = 1.0,
							double gamma = 0.99,
							int maxiter = 100,
							bool finalize = false) {

	sp_mat L_matched = NetLibR::netAlign_arma(A, B, L, alpha, beta, gamma, maxiter, finalize);
	
	return(L_matched);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat PCSF(sp_mat &Adj, vec node_scores, double kappa = 1.0, int root = -1, int clusters = 20, int convert_similarity_to_distance = 0) {
	
	sp_mat G = NetLibR::PCSF(Adj, node_scores, root, clusters, convert_similarity_to_distance);

	return(G);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List computeAutocorrelation_Geary (sp_mat &G, mat &scores, int rand_perm = 100, int num_shuffles = 10000) {
	
	field<vec> results = NetLibR::computeAutocorrelation_Geary (G, scores, rand_perm, num_shuffles);
	

	List out_list;		
	out_list["C"] = results(0);
	out_list["mu"] = results(1);
	out_list["sigma"] = results(2);
	out_list["Z"] = results(3);
	
    return out_list;
}






/*

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List buildAdaptiveNet_fromDist(mat &D, int max_kNN = -1, double LC = 1.0, int sym_method = 0, int thread_no=4, int dist_type = -1, double sigma = 0.5) {	
	
	field<sp_mat> NN = NetLibR::buildAdaptiveNet_fromDist(D, max_kNN, LC, sym_method, thread_no, dist_type, sigma);

	List out_list;		
	out_list["sym"] = NN(0);
	out_list["org"] = NN(1);

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List buildAdaptiveNet_fromEdgeList(mat &dist, mat &idx, int max_kNN = -1, double LC = 1.0, int sym_method = 0, int thread_no=4, int dist_type = -1, double sigma = 0.5) {	
	
	idx -= 1;
	
	field<sp_mat> NN = NetLibR::buildAdaptiveNet_fromEdgeList(dist, idx, max_kNN, LC, sym_method, thread_no, dist_type, sigma);

	List out_list;		
	out_list["sym"] = NN(0);
	out_list["org"] = NN(1);

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat make_spanner(sp_mat &G_adj, int k) {
	sp_mat G_sparse = NetLibR::make_spanner(G_adj, k);

    return G_sparse;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat netAlign(sp_mat A, sp_mat B, sp_mat L,
							double alpha = 1.0,
							double beta = 1.0,
							double gamma = 0.99,
							int maxiter = 100,
							bool finalize = false) {

	sp_mat L_matched = NetLibR::netAlign_arma(A, B, L, alpha, beta, gamma, maxiter, finalize);
	
	return(L_matched);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat PCSF(sp_mat &Adj, vec node_scores, double kappa = 1.0, int root = -1, int clusters = 20, int convert_similarity_to_distance = 0) {
	
	sp_mat G = NetLibR::PCSF(Adj, node_scores, root, clusters, convert_similarity_to_distance);

	return(G);
}
	
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat batchPR(sp_mat &G, mat &U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {	
	mat U_smoothed = NetLibR::batchPR(G, U, alpha, thread_no, tol);
	
    return U_smoothed;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat zoned_diffusion(sp_mat &G, uvec& zones, mat &U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {	
	mat U_smoothed = NetLibR::batch_zoned_diffusion(G, zones, U, alpha, thread_no, tol);
	
    return U_smoothed;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec sweepcut(sp_mat &A, vec s) {
    vec conductance = NetLibR::sweepcut(A, s);

    return conductance;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List computeAutocorrelation_Geary (sp_mat &G, mat &scores, int rand_perm = 100, int num_shuffles = 10000) {
	
	field<vec> results = NetLibR::computeAutocorrelation_Geary (G, scores, rand_perm, num_shuffles);
	

	List out_list;		
	out_list["C"] = results(0);
	out_list["mu"] = results(1);
	out_list["sigma"] = results(2);
	out_list["Z"] = results(3);
	
    return out_list;
}
*/
