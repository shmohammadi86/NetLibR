// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/NetLibR.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// set_seed
void set_seed(double seed);
RcppExport SEXP _NetLibR_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// MWM_hungarian
mat MWM_hungarian(mat& G);
RcppExport SEXP _NetLibR_MWM_hungarian(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(MWM_hungarian(G));
    return rcpp_result_gen;
END_RCPP
}
// MWM_rank1
umat MWM_rank1(vec u, vec v, double u_threshold, double v_threshold);
RcppExport SEXP _NetLibR_MWM_rank1(SEXP uSEXP, SEXP vSEXP, SEXP u_thresholdSEXP, SEXP v_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type u_threshold(u_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type v_threshold(v_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(MWM_rank1(u, v, u_threshold, v_threshold));
    return rcpp_result_gen;
END_RCPP
}
// signed_cluster
vec signed_cluster(sp_mat A, double resolution_parameter, int seed, Nullable<IntegerVector> initial_clusters_);
RcppExport SEXP _NetLibR_signed_cluster(SEXP ASEXP, SEXP resolution_parameterSEXP, SEXP seedSEXP, SEXP initial_clusters_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type resolution_parameter(resolution_parameterSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type initial_clusters_(initial_clusters_SEXP);
    rcpp_result_gen = Rcpp::wrap(signed_cluster(A, resolution_parameter, seed, initial_clusters_));
    return rcpp_result_gen;
END_RCPP
}
// unsigned_cluster
vec unsigned_cluster(sp_mat A, double resolution_parameter, int seed, Nullable<IntegerVector> initial_clusters_);
RcppExport SEXP _NetLibR_unsigned_cluster(SEXP ASEXP, SEXP resolution_parameterSEXP, SEXP seedSEXP, SEXP initial_clusters_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type resolution_parameter(resolution_parameterSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type initial_clusters_(initial_clusters_SEXP);
    rcpp_result_gen = Rcpp::wrap(unsigned_cluster(A, resolution_parameter, seed, initial_clusters_));
    return rcpp_result_gen;
END_RCPP
}
// bmatching
sp_mat bmatching(sp_mat A, vec b_limit, int seed);
RcppExport SEXP _NetLibR_bmatching(SEXP ASEXP, SEXP b_limitSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< vec >::type b_limit(b_limitSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(bmatching(A, b_limit, seed));
    return rcpp_result_gen;
END_RCPP
}
// netAlign
sp_mat netAlign(sp_mat A, sp_mat B, sp_mat L, double alpha, double beta, double gamma, int maxiter, bool finalize);
RcppExport SEXP _NetLibR_netAlign(SEXP ASEXP, SEXP BSEXP, SEXP LSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP maxiterSEXP, SEXP finalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< sp_mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< sp_mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type finalize(finalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(netAlign(A, B, L, alpha, beta, gamma, maxiter, finalize));
    return rcpp_result_gen;
END_RCPP
}
// PCSF
sp_mat PCSF(sp_mat& Adj, vec node_scores, int root, int clusters, int convert_similarity_to_distance);
RcppExport SEXP _NetLibR_PCSF(SEXP AdjSEXP, SEXP node_scoresSEXP, SEXP rootSEXP, SEXP clustersSEXP, SEXP convert_similarity_to_distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< vec >::type node_scores(node_scoresSEXP);
    Rcpp::traits::input_parameter< int >::type root(rootSEXP);
    Rcpp::traits::input_parameter< int >::type clusters(clustersSEXP);
    Rcpp::traits::input_parameter< int >::type convert_similarity_to_distance(convert_similarity_to_distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(PCSF(Adj, node_scores, root, clusters, convert_similarity_to_distance));
    return rcpp_result_gen;
END_RCPP
}
// computeAutocorrelation_Geary
List computeAutocorrelation_Geary(sp_mat& G, mat& scores, int rand_perm, int num_shuffles);
RcppExport SEXP _NetLibR_computeAutocorrelation_Geary(SEXP GSEXP, SEXP scoresSEXP, SEXP rand_permSEXP, SEXP num_shufflesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< sp_mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< mat& >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< int >::type rand_perm(rand_permSEXP);
    Rcpp::traits::input_parameter< int >::type num_shuffles(num_shufflesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeAutocorrelation_Geary(G, scores, rand_perm, num_shuffles));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NetLibR_set_seed", (DL_FUNC) &_NetLibR_set_seed, 1},
    {"_NetLibR_MWM_hungarian", (DL_FUNC) &_NetLibR_MWM_hungarian, 1},
    {"_NetLibR_MWM_rank1", (DL_FUNC) &_NetLibR_MWM_rank1, 4},
    {"_NetLibR_signed_cluster", (DL_FUNC) &_NetLibR_signed_cluster, 4},
    {"_NetLibR_unsigned_cluster", (DL_FUNC) &_NetLibR_unsigned_cluster, 4},
    {"_NetLibR_bmatching", (DL_FUNC) &_NetLibR_bmatching, 3},
    {"_NetLibR_netAlign", (DL_FUNC) &_NetLibR_netAlign, 8},
    {"_NetLibR_PCSF", (DL_FUNC) &_NetLibR_PCSF, 5},
    {"_NetLibR_computeAutocorrelation_Geary", (DL_FUNC) &_NetLibR_computeAutocorrelation_Geary, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_NetLibR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
