#include <NetLibR.h>

// Algorithm adopted and reimplemented from:
// Anava, O., & Levy, K. (2016). k\ast -Nearest Neighbors: From Global to Local. In D. D. Lee, M. Sugiyama, U. V Luxburg, I. Guyon, & R. Garnett (Eds.), Advances in Neural Information Processing Systems 29 (pp. 4916–4924).
// https://github.com/kfirkfir/k-Star-Nearest-Neighbors
namespace NetLibR {
	field<sp_mat> buildAdaptiveNet_fromDist(mat &D, int max_kNN = -1, double LC = 1.0, int sym_method = AND_SYM, int thread_no=4, int dist_type = -1, double sigma = 0.5) {
		int nV = D.n_rows;
		printf("buildAdaptiveNet_fromD:: Running k*-nearest neighbors algorithm with given distance matrix D over %d vertices (LS = %.2f)\n", nV, LC);		

		if(max_kNN == -1) {
			double kappa = 5; // Random, really!
			max_kNN = min(nV-2, (int)(kappa*round(sqrt(nV))));
		}
		
		D.diag().zeros();
		
		printf("\tConverting from distance matrix to neighborhood list ... ");
		mat dist = zeros(max_kNN+1, nV);
		mat idx  = zeros(max_kNN+1, nV);	
		
		for(int i = 0; i < nV; i++) {
			vec d = D.col(i);
			uvec indices = sort_index(d);
			
			idx.col(i) = conv_to< vec >::from(indices(span(0, max_kNN))); 			
			dist.col(i) = d(indices(span(0, max_kNN)));
		}
		dist = trans(dist);
		idx = trans(idx);
		printf("done\n");
		
		printf("\tComputing Lambda ...");
		mat beta = LC*dist;
		vec beta_sum = zeros(nV);
		vec beta_sq_sum = zeros(nV);
		
		register int k;
		mat lambda = zeros(size(beta));
		for(k = 1; k <= max_kNN; k++) {
			beta_sum += beta.col(k);
			beta_sq_sum += square(beta.col(k));
			
			lambda.col(k) = (1.0/(double)k) * ( beta_sum + sqrt(k + square(beta_sum) - k*beta_sq_sum) );
		}
		lambda.replace(datum::nan, 0); 			
		
		lambda = trans(lambda);
		beta = trans(beta);
		
		printf("done\n");
		

		printf("\tComputing the number of nearest neighbors ... ");
		mat Delta = lambda - beta;		
		Delta.shed_row(0);

		vec NN_count = zeros(nV);
		vec mean_dist = zeros(nV); // for generic distance <-> matrix based on the Network Similarity Fusion (NSF) algorithm
		for(int v = 0; v < nV; v++) {				
			vec delta = Delta.col(v);		
					
			uvec rows = find(delta < 0, 1, "first");
			NN_count(v) = rows.n_elem == 0?max_kNN:(rows(0));
			
			rowvec d = dist.row(v);
			mean_dist(v) = mean(d(span(1, NN_count(v))));			
		}
		printf("done\n");
		
		printf("\tConstructing k*-NN graph ... ");
//		# pragma omp parallel for shared(G) num_threads(thread_no)

//		# pragma omp parallel num_threads(thread_no) {
			sp_mat G(nV, nV);
	//		# pragma omp for shared(G)
			for(int v = 0; v < nV; v++) {				
				int neighbor_no = NN_count(v);
				
				int dst = v;								
				rowvec v_dist = dist.row(v);
				rowvec v_idx = idx.row(v);			
				
				for (int i = 1; i <= neighbor_no; i++) {				
					int src = v_idx(i);
						
					if( (dist_type == COR_DIST) || (dist_type == ACTIONet_DIST) || (dist_type == JSD_DIST)) {
						G(src, dst) = 1.0 - v_dist(i);
					} 
					else {
						double denom = sigma*(mean_dist(src) + mean_dist(dst) + v_dist(i)) / 3;
						G(src, dst) = exp( -(v_dist(i)*v_dist(i)) / denom);
					}
				} 
			}
		//}
		printf("done\n");
		
		
		printf("\tSymmetrizing network ...  (sym method = %d)", sym_method);
		G.replace(datum::nan, 0);  // replace each NaN with 0

		sp_mat Gt = trans(G);	
		
		sp_mat G_sym;
		if(sym_method == AND_SYM) {
			G_sym = sqrt(G % Gt);
		} else if(sym_method == OR_SYM) {
			G_sym = (G + Gt);
			G_sym.for_each( [](sp_mat::elem_type& val) { val /= 2.0; } );			
		} else { // Default to MNN
			sp_mat G_sym = sqrt(G % Gt);
		}
		sp_mat G_asym = normalise(G, 1, 0);
		G_sym = G_sym / (double)max(max(G_sym));
		
		field<sp_mat> output(2);
		output(0) =	G_sym;
		output(1) = G_asym;
		
		printf("done\n");
		
		return(output);		
	}



	field<sp_mat> buildAdaptiveNet_fromEdgeList(mat &dist, mat &idx, int max_kNN = -1, double LC = 1.0, int sym_method = AND_SYM, int thread_no=4, int dist_type = -1, double sigma = 0.5) {
		int nV = dist.n_rows;
		printf("buildAdaptiveNet:: Running k*-nearest neighbors algorithm with given distance matrix D over %d vertices (LS = %.2f)\n", nV, LC);		

		if(max_kNN == -1) {
			double kappa = 5; // Random, really!
			int max_kNN = min(nV-1, (int)(kappa*round(sqrt(nV))));
		}		

		mat beta = LC*dist;
		vec beta_sum = zeros(nV);
		vec beta_sq_sum = zeros(nV);
		
		register int k;
		mat lambda = zeros(size(beta));
		for(k = 1; k <= max_kNN; k++) {
			beta_sum += beta.col(k);
			beta_sq_sum += square(beta.col(k));
			
			lambda.col(k) = (1.0/(double)k) * ( beta_sum + sqrt(k + square(beta_sum) - k*beta_sq_sum) );
		}
		lambda.replace(datum::nan, 0); 			
		
		lambda = trans(lambda);
		beta = trans(beta);
		
		printf("done\n");
		

		printf("\tComputing the number of nearest neighbors ... ");
		mat Delta = lambda - beta;		
		Delta.shed_row(0);

		vec NN_count = zeros(nV);
		vec mean_dist = zeros(nV); // for generic distance <-> matrix based on the Network Similarity Fusion (NSF) algorithm
		for(int v = 0; v < nV; v++) {				
			vec delta = Delta.col(v);		
					
			uvec rows = find(delta < 0, 1, "first");
			NN_count(v) = rows.n_elem == 0?max_kNN:(rows(0));
			
			rowvec d = dist.row(v);
			mean_dist(v) = mean(d(span(1, NN_count(v))));			
		}
		printf("done\n");
		
		printf("\tConstructing k*-NN graph ... ");
//		# pragma omp parallel for shared(G) num_threads(thread_no)

//		# pragma omp parallel num_threads(thread_no) {
			sp_mat G(nV, nV);
			
	//		# pragma omp for shared(G)
			for(int v = 0; v < nV; v++) {				
				int neighbor_no = NN_count(v);
				
				int dst = v;								
				rowvec v_dist = dist.row(v);
				rowvec v_idx = idx.row(v);			
				
				for (int i = 1; i <= neighbor_no; i++) {				
					int src = v_idx(i);
					
					if( (dist_type == COR_DIST) || (dist_type == ACTIONet_DIST) || (dist_type == JSD_DIST)) {
						G(src, dst) = 1.0 - v_dist(i);
					} 
					else {
						double denom = sigma*(mean_dist(src) + mean_dist(dst) + v_dist(i)) / 3;
						G(src, dst) = exp( -(v_dist(i)*v_dist(i)) / denom);
					}
				} 
			}
		//}
		printf("done\n");
		
		
		printf("\tSymmetrizing network ...  (sym method = %d)", sym_method);
		G.replace(datum::nan, 0);  // replace each NaN with 0

		sp_mat Gt = trans(G);	
		
		sp_mat G_sym;
		if(sym_method == AND_SYM) {
			G_sym = sqrt(G % Gt);
		} else if(sym_method == OR_SYM) {
			G_sym = (G + Gt);
			G_sym.for_each( [](sp_mat::elem_type& val) { val /= 2.0; } );			
		} else { // Default to MNN
			sp_mat G_sym = sqrt(G % Gt);
		}
		sp_mat G_asym = normalise(G, 1, 0);
		G_sym = G_sym / (double)max(max(G_sym));
		
		field<sp_mat> output(2);
		output(0) =	G_sym;
		output(1) = G_asym;
		
		printf("done\n");
		
		return(output);		
	}

}






