#include <NetLibR.h>

#include <prpack_utils.h>
#include <prpack_result.h>
#include <prpack_solver.h>

vec pagerank(sp_mat &G, vec u_vec, double alpha = 0.85, double tol = 1e-6) {

	prpack::prpack_base_graph g(&G);

	g.normalize_weights(); 
	
	
    prpack::prpack_solver solver(&g, false);

	u_vec = normalise(u_vec, 1);


	double* u = new double[u_vec.size()];
	memcpy(u, u_vec.memptr(), u_vec.size()*sizeof(double));
    double* v = u;
    
    
    const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");
    

    vec pr(res->x, res->num_vs);
        
    return pr;
}


vec pagerank_scaled(sp_mat &G, vec u_vec, double alpha = 0.85, double tol = 1e-6) {
	
	prpack::prpack_base_graph g(&G);
	g.normalize_weights(); 
    prpack::prpack_solver solver(&g, false);

	u_vec = normalise(u_vec, 1);
	
	double* u = new double[u_vec.size()];
	memcpy(u, u_vec.memptr(), u_vec.size()*sizeof(double));
    double* v = u;
    
    const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");    
    vec pr(res->x, res->num_vs);
    
    
    vec o = ones(size(u_vec)) / (double)u_vec.n_elem;	
	memcpy(u, o.memptr(), o.size()*sizeof(double));
    v = u;
        
    res = solver.solve(alpha, tol, u, v, "");    
    vec pr0(res->x, res->num_vs);
    
    vec log_ratio = log2(pr / pr0);
    
    log_ratio.replace(datum::nan, 0); 

    return log_ratio;
}

namespace NetLibR {
	mat batch_zoned_diffusion(sp_mat &G, uvec &zones, mat &U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {
		mat ZonedDiff = zeros(size(U));
		
		if(U.n_rows != G.n_rows) {
			fprintf(stderr, "batch_zoned_diffusion:: Number of rows in U doesn't match with the number of vertices in G\n");
			return(ZonedDiff);
		}
		U.transform( [](double val) { return (val < 0? 0:val); } );
		U = normalise(U, 1);
		
		G.diag().zeros();
		
		uvec uz = unique(zones);					
		printf("Zoned diffusion (%d zones)\n", uz.n_elem);
		
		uvec vertex_id(size(zones));	
		for(int i = 0; i < uz.n_elem; i++) {
			uvec idx = find(zones == uz(i));

			for(int j = 0; j < idx.n_elem; j++) {
				vertex_id(idx(j)) = j;
			}
					
			sp_mat subG(idx.n_elem, idx.n_elem);
			
			sp_mat::iterator it     = G.begin();
			sp_mat::iterator it_end = G.end();
			for(; it != it_end; ++it) {
				if( (zones(it.row()) != uz(i)) || (zones(it.col()) != uz(i)) )
					continue;
				
				int ii = vertex_id(it.row());
				int jj = vertex_id(it.col());

				subG(ii, jj) = subG(jj, ii) = (*it);
			}
			
			prpack::prpack_base_graph g(&subG);		
			g.normalize_weights(); 
			prpack::prpack_solver solver(&g, false);

			
			#pragma omp parallel for num_threads(thread_no) 			
			for(int k = 0; k < U.n_cols; k++) {
				vec init_pr = U.col(k);
				//double total_max = max(init_pr);
				double weight = sum(init_pr(idx));				
				
				
				vec sub_init_pr = normalise(init_pr(idx), 1);
				double* u = new double[sub_init_pr.size()];
				memcpy(u, sub_init_pr.memptr(), sub_init_pr.size()*sizeof(double));
				double* v = u;
			
				const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");
											
				vec pr(res->x, res->num_vs);
				pr.replace(datum::nan, 0);

				vec x = ZonedDiff.col(k);
				x(idx) = weight*pr;
				ZonedDiff.col(k) = x;
			}			
		}

		ZonedDiff = normalise(ZonedDiff, 1);
		return(ZonedDiff);	
	}


	
	
	
	
	mat batchPR(sp_mat &G, mat &U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {
		
		prpack::prpack_base_graph g(&G);		
		
		g.normalize_weights(); 
		prpack::prpack_solver solver(&g, false);
		
		U = normalise(U, 1);
		
		mat PR = zeros(size(U));
		
		
		#pragma omp parallel for num_threads(thread_no) 			
		for(int i = 0; i < U.n_cols; i++) {			
			vec u_vec = U.col(i);
			double* u = new double[u_vec.size()];
			memcpy(u, u_vec.memptr(), u_vec.size()*sizeof(double));
			double* v = u;
			
			const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");
			
			vec pr(res->x, res->num_vs);
			PR.col(i) = pr;
		}	
		
		PR.replace(datum::nan, 0);
		uvec isolated_vertices = find(sum(PR) < 1e-6);
		PR.cols(isolated_vertices) = normalise(ones(PR.n_rows, isolated_vertices.n_elem), 1);
		
		return PR;
	}

	


	vec sweepcut(sp_mat &A, vec s) {
		int top_ignore = 5;
		
		A.diag().zeros();		
		int nV = A.n_rows;
		
		vec w = vec(sum(A, 1));
		double total_vol = sum(w);

		vec conductance = datum::inf*ones(w.n_elem);				
		
		uvec perm = sort_index( s, "descend" );
		vec x = zeros(nV);
		x(perm(span(0, top_ignore-1))).ones();
		double vol = sum(w(perm(span(0, top_ignore-1))));

		double cut_size = vol;
		for(int i = 0; i < top_ignore; i++) {
			for(int j = 0; j < top_ignore; j++) {
				cut_size -= A(i, j);
			}
		}
		
		for(register int i = top_ignore; i < nV-top_ignore-1; i++) {
			int u = perm(i);
			vol += w[u];
			
			x(u) = 1;			

			sp_mat::col_iterator it = A.begin_col( u );							
			for(; it != A.end_col( u ); it++) {
				int v = it.row();
				if(x[v] == 0) { // v is in S_prime (not yet selected)
					cut_size += (*it);
				}	
				else {
					cut_size -= (*it);
				}			
			}
						
			double vol_prime = total_vol - vol;
			conductance(i) = cut_size / min(vol, vol_prime);
		}

		
		return(conductance);
	}
	
	vec nonlinear_diffusion(sp_mat A, vec s, double p = 0.5, double h = 0.001, int k = 10) {
		sp_mat A_norm = normalise(A, 1, 0);
		sp_mat I = speye(size(A));
		sp_mat L = I - A_norm;
		
		printf("h = %e, p = %e\n", h, p);
		
		vec u = s;		
		for(register int i = 0; i < k; i++) {
			vec delta = h * L * pow(u, p);
			u = u - delta;
			u = clamp( u, 2.0e-8, 1 );
		}
		
		return(u);
	}



}
