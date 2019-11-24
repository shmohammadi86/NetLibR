#include <NetLibR.h>

namespace NetLibR {
	// G is the cell-cell network, scores is a cell x geneset matrix
	field<vec> computeAutocorrelation_Geary (sp_mat &G, mat &scores, int rand_perm_no, int num_shuffles = 1000) {
		int nV = G.n_rows;
		int feature_set_no = scores.n_cols;
		
		printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, feature_set_no);
		int nnz = G.n_nonzero;
		int idx = 0;		

		vec vals(nnz);
		umat subs(2, nnz);
		for(sp_mat::iterator it= G.begin(); it != G.end(); ++it) {
			vals(idx) = *it;
			subs(0, idx) = it.row();
			subs(1, idx) = it.col();
			idx++;
		}

		//double scale_factor = sum(sum(G)) / (sample_no-1); // 2W / (N-1)
		double total_weight = sum(vals);

		// Compute graph Laplacian
		vec d = vec(trans(sum(G)));
		
		sp_mat L(nV, nV);
		L.diag() = d;
		L -= G;
		
		// Summary stats
		vec Cstat = zeros(feature_set_no);


		
		printf("Computing autocorrelations ..."); fflush(stdout);
		for(int i = 0; i < feature_set_no; i++) {
			vec x = scores.col(i);
			double stat = dot(x, L*x); // x'Lx = sum_{i,j} (w_{ij}(x_i - x_j)^2)						
			double norm_fact = var(x)*total_weight;
			Cstat(i) = 1 - (stat / norm_fact);
		}
		printf("done\n"); fflush(stdout);
		
		
		vec mu, sigma, Cstat_Z;
		if(rand_perm_no > 0) {
			int perc = 0;
			mat rand_stat = zeros(rand_perm_no, feature_set_no);
			printf("Computing significance of autocorrelations ...\n"); fflush(stdout);		
			for(int k = 0; k < rand_perm_no; k++) {
				if(round(100.0*k / rand_perm_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}
				
				
				/*
				sp_mat G_rand = shuffleEdges(subs, vals, nV, num_shuffles, 1000*k);			
								
				sp_mat L_rand(nV, nV);	
				L_rand.diag() = d;
				L_rand -= G_rand;
				*/
				
				for(int i = 0; i < feature_set_no; i++) {
					vec x = scores.col(i);
					uvec rand_perm = sort_index(randu(size(x)));
					vec x_permuted = x(rand_perm);					
					
					double stat = dot(x_permuted, L*x_permuted); // x'Lx = sum_{i,j} (w_{ij}(x_i - x_j)^2)						
					double norm_fact = var(x_permuted)*total_weight;
					rand_stat(k, i) = 1 - (stat / norm_fact);
				}
			}
			printf("done\n"); fflush(stdout);
			
			mu = trans(mean(rand_stat));
			sigma = trans(stddev(rand_stat));
			Cstat_Z = (Cstat - mu) / sigma;
		}
		else {
			mu = zeros(feature_set_no);
			sigma = zeros(feature_set_no);
			Cstat_Z = zeros(feature_set_no);
		}
		
		field<vec> results(4);
		results(0) = Cstat;		
		results(1) = mu;
		results(2) = sigma;
		results(3) = Cstat_Z;
		
		return(results);
	}
}
