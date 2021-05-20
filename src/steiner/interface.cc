#include <NetLibR.h>

#include <pcst_fast.h>

using namespace cluster_approx;

void WriteToStderr(const char* s) {
  fprintf(stderr, "%s", s);
  fflush(stderr);
}

namespace NetLibR {
	sp_mat PCSF(sp_mat &Adj, vec node_scores, int root = -1, int clusters = 20, int convert_similarity_to_distance = 0) {
		int V = Adj.n_rows;

		vector <double> prizes(node_scores.memptr(), node_scores.memptr()+node_scores.n_elem);
		
		
		// Extract network edges and convert edge weight to distances/costs
		sp_mat G = trimatl(Adj);
		sp_mat::iterator it     = G.begin();
		sp_mat::iterator it_end = G.end();

		vector<double> costs(G.n_nonzero);
		vector< pair<int, int> > edges(G.n_nonzero);
		for(; it != it_end; ++it) {
			//if( node_scores(it.row()) == 0 || node_scores(it.col()) == 0 )
				//continue;
				
			edges.push_back(make_pair(it.row(), it.col()));
			
			switch(convert_similarity_to_distance) {
				case 0:
					costs.push_back(*it);
					break;
					
				case 1: // Usefull for ACTIONet weights
					costs.push_back(1.0 - (*it));
					break;
					
				case 2: // Useful for general distance
					costs.push_back(1.0 / (*it));
					break;
				
				case 3: // Useful for p-values
					costs.push_back( -log(*it) );
			}
		}
		
		PCSTFast *algo;
		if( (root == -1) || (root >= Adj.n_rows) ) {
			printf("PCST: unrooted (%d clusters)\n", clusters);
			algo = new PCSTFast(edges, prizes, costs, -1, clusters, PCSTFast::kStrongPruning, 0, WriteToStderr);
		} else {
			printf("PCST: rooted at v%d\n", root+1);
			algo = new PCSTFast(edges, prizes, costs, root, 0, PCSTFast::kStrongPruning, 0, WriteToStderr);
		}
		
		vector<int> node_result;
		vector<int> edge_result;		
		algo->run(&node_result, &edge_result);
		
		sp_mat out(size(Adj));
		for(int i = 0; i < edge_result.size(); i++) {
			int src = edges[edge_result[i]].first;
			int dst = edges[edge_result[i]].second;
			
			out(src, dst) = Adj(src, dst);			
		}
		
		return(out);
	}
}
