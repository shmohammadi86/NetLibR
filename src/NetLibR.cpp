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
vec signed_cluster(sp_mat A, double resolution_parameter = 1.0, int seed = 0, Nullable<IntegerVector> initial_clusters_ = R_NilValue) {
    set_seed(seed);

	uvec initial_clusters_uvec(A.n_rows);
	if ( initial_clusters_.isNotNull() ) {
        NumericVector initial_clusters(initial_clusters_);
		
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = initial_clusters(i);
	} else {
		for (int i = 0; i < A.n_rows; i++) initial_clusters_uvec(i) = i;
	}

	
	vec clusters = NetLibR::signed_cluster(A, resolution_parameter, seed, initial_clusters_uvec);

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

	vec clusters = NetLibR::unsigned_cluster(A, resolution_parameter, seed, initial_clusters_uvec);

    return clusters;	
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
// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat computeFullDist(mat &H_stacked, int thread_no = 8, int verbose = 1) {  	
    mat D = NetLibR::computeFullDist(H_stacked, thread_no, verbose);
    
    return D;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
sp_mat computeNearestDist(mat &H_stacked, double kNN, int thread_no = 8) {  
    sp_mat D = NetLibR::computeNearestDist(H_stacked, kNN, thread_no);
    
    return D;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List computeNearestDist_edgeList(mat &H_stacked, double kNN, int thread_no = 8) {  	
    field<mat> NN = NetLibR::computeNearestDist_edgeList(H_stacked, kNN, thread_no);
    
	List out_list;		
	out_list["idx"] = NN(0);
	out_list["dist"] = NN(1);

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List buildAdaptiveACTIONet(mat &H_stacked, double LC = 1.0, double epsilon = 0.0, int thread_no = 8, bool auto_adjust_LC = false, string sym_method = "OR") {  
	
	int sm = sym_method == "AND"? ACTIONet_AND:ACTIONet_OR;
	
	printf("Sym method: %s (%d)\n", sym_method.c_str(), sm);
		
    field<sp_mat> res = NetLibR::buildAdaptiveACTIONet(H_stacked, LC, epsilon, thread_no, auto_adjust_LC, sm);

	List out_list;		
	out_list["ACTIONet"] = res(0);
	out_list["ACTIONet_asym"] = res(1);

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List layoutACTIONet(sp_mat G,
	mat S_r,
	int compactness_level= 50,
	unsigned int n_epochs = 100,
	int thread_no = 8) {

	field<mat> res = NetLibR::layoutACTIONet(G, S_r, compactness_level, n_epochs, thread_no);
    
	List out_list;		
	out_list["coordinates"] = res(0);
	out_list["coordinates_3D"] = res(1);
	out_list["colors"] = res(2);

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List reconstructArchetypes(sp_mat &S, const List& C_trace, const List& H_trace, double z_threshold = 1.0) {
	
	int n_list = H_trace.size();
	vector<mat> C_trace_vec(n_list+1);
	vector<mat> H_trace_vec(n_list+1);
	for (int i = 0; i < n_list; i++) {
		if(Rf_isNull(H_trace[i])) {
			continue;
		}
		C_trace_vec[i+1] = (as<mat>(C_trace[i]));
		H_trace_vec[i+1] = (as<mat>(H_trace[i]));
	}

	NetLibR::multilevel_archetypal_decomposition results = NetLibR::reconstructArchetypes(S, C_trace_vec, H_trace_vec, z_threshold = -1.0);

	List out_list;		
	
	
	for(int i = 0; i < results.selected_archs.n_elem; i++) results.selected_archs[i]++;
	out_list["selected_archs"] = results.selected_archs; 

	for(int i = 0; i < results.landmark_cells.n_elem; i++) results.landmark_cells[i]++;
	out_list["landmark_cells"] = results.landmark_cells; 

	out_list["C_stacked"] = results.C_stacked;
	out_list["H_stacked"] = results.H_stacked;
	
	out_list["archetype_profile"] = results.archetype_profile;

	out_list["backbone"] = results.backbone;

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat extractArchetypeAssociatedSamples(sp_mat &G, mat &H_stacked, double alpha = 0.85) {

	mat scores = NetLibR::extractArchetypeAssociatedSamples(G, H_stacked, alpha);
	
	return scores;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets(sp_mat &S, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = NetLibR::assessFeatureSets(S, feature_sets, rand_perm);		

    return scores;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets_archs(mat &archetype_profile, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = NetLibR::assessFeatureSets_archs(archetype_profile, feature_sets, rand_perm);		

    return scores;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets_decoupled(mat archetype_profile, mat H_stacked, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = NetLibR::assessFeatureSets_decoupled(archetype_profile, H_stacked, feature_sets, rand_perm);	

    return scores;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List computeAutocorrelation (sp_mat &G, mat &scores, int rand_perm = 100, int num_shuffles = 10000) {
	
	field<vec> results = NetLibR::computeAutocorrelation (G, scores, rand_perm, num_shuffles);
	

	List out_list;		
	out_list["C"] = results(0);
	out_list["mu"] = results(1);
	out_list["sigma"] = results(2);
	out_list["Z"] = results(3);
	
    return out_list;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat phenotypeEnrichment(mat &H_stacked, mat &phenotype_associations, int rand_perm_no) {	
	mat Z = NetLibR::phenotypeEnrichment (H_stacked, phenotype_associations, rand_perm_no);

    return Z;
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
sp_mat mergeArchetypes(mat C_stacked, mat H_stacked) {
	sp_mat results = NetLibR::mergeArchetypes(C_stacked, H_stacked);

    return results;	
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

	
	vec clusters = NetLibR::signed_cluster(A, resolution_parameter, seed, initial_clusters_uvec);

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

	vec clusters = NetLibR::unsigned_cluster(A, resolution_parameter, seed, initial_clusters_uvec);

    return clusters;	
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List constructBackbone(mat arch_profile_reduced, double weight_threshold = 0, double pval_threshold = 0.05, double lambda = -1, int thread_no = 8) {
	field<mat> results = NetLibR::FastGGM(arch_profile_reduced, lambda, thread_no);		
	
	mat weights = results(0);
	mat pvals = results(1);

	mat backbone = weights;
	uvec idx = find(weights < weight_threshold);
	backbone(idx).zeros();
	
	idx = find(pvals > pval_threshold);
	backbone(idx).zeros();


	List out_list;		
	
	
	out_list["backbone"] = backbone; 
	out_list["weights"] = weights; 
	out_list["pvals"] = pvals; 

    return out_list;	
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat update_layout_2D(mat coors,
	int compactness_level= 50,
	unsigned int n_epochs = 100,
	int thread_no = 8) {
	
	if(coors.n_cols == 2)
		coors = trans(coors);
	
		
	mat updated_coors = NetLibR::update_layout_2D(coors, compactness_level, n_epochs, thread_no);

    return updated_coors;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat constructKstarNN_fromDist(mat &D, double L_C = 1) {
		
	mat G = NetLibR::constructKstarNN_fromDist(D, L_C);

    return G;
}

*/
