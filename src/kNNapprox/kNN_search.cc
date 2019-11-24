#include <NetLibR.h>
#include <atria.h>
#include <vptree.h>

namespace NetLibR {
	sp_mat computeNearestDist(mat &X, int kNN, double (*distance)( const T&, const T& ), int thread_no=4) {	
		printf("Building distance matrix of the %d nearest neighbors of each node (returns sparse distance matrix)\n", kNN);
		
		double epsilon = 1e-10;
		int sample_no = X.n_rows;
		
		if(kNN >= sample_no || kNN <= 1)
			kNN = min(30, sample_no-1);
			
		umat subs(2, kNN*sample_no);
		vec vv(kNN*sample_no);

		// build ball tree on set
		VpTree<DataPoint, distance>* tree = new VpTree<DataPoint, distance>();
		std::vector<DataPoint> samples(sample_no); //, DataPoint(archetype_no, -1, data));
		for (int i = 0; i < sample_no; i++) {
			samples[i] = DataPoint(sample_no, i, X.rowptr(i));
		}
		tree -> create(samples, 0);
		
		
		int perc = 1;
		int total_counts = 1;
		#pragma omp parallel num_threads(thread_no) 
		{
			vector< vector<int> > ind_arr(sample_no, std::vector<int>(kNN+1));
			vector< vector<double> > dist_arr(sample_no, std::vector<double>(kNN+1));
			#pragma omp for
			for (int v = 0; v < sample_no; v++) {
				total_counts ++;
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}
				
				tree -> search(samples[v], kNN+1, &ind_arr[v], &dist_arr[v]);
				
				int base = v*kNN;
				for(int i = 1; i <= kNN; i++) {
					double d = dist_arr[v][i]; // To match with old scores
					d = d < epsilon?epsilon:d;
						
					subs(0, base + i-1) = ind_arr[v][i]-1;
					subs(1, base + i-1) = v;
					vv(base + i-1) = d;			
				}				
			}
		}
		samples.clear();
		delete tree;
			
		sp_mat D(subs, vv, sample_no, sample_no);	

		return(D);
	}
	

	field<mat> computeNearestDist_edgeList(mmat &X, int kNN, double (*distance)( const T&, const T& ), int thread_no=4) {	
		printf("Building distance matrix of the %d nearest neighbors of each node (returns edge list)\n", kNN);
		
		double epsilon = 1e-10;
		int sample_no = X.n_rows;
		
		if(kNN >= sample_no || kNN <= 1)
			kNN = min(30, sample_no-1);
			
		umat subs(2, kNN*sample_no);
		vec vv(kNN*sample_no);

		// build ball tree on set
		VpTree<DataPoint, distance>* tree = new VpTree<DataPoint, distance>();
		std::vector<DataPoint> samples(sample_no); //, DataPoint(archetype_no, -1, data));
		for (int i = 0; i < sample_no; i++) {
			samples[i] = DataPoint(sample_no, i, X.rowptr(i));
		}
		tree -> create(samples, 0);
		
		
		int perc = 1;
		int total_counts = 1;
		#pragma omp parallel num_threads(thread_no) 
		{
			vector< vector<int> > ind_arr(sample_no, std::vector<int>(kNN+1));
			vector< vector<double> > dist_arr(sample_no, std::vector<double>(kNN+1));
			#pragma omp for
			for (int v = 0; v < sample_no; v++) {
				total_counts ++;
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}
				
				tree -> search(samples[v], kNN+1, &ind_arr[v], &dist_arr[v]);
				
				int base = v*kNN;
				for(int i = 1; i <= kNN; i++) {
					double d = dist_arr[v][i]; // To match with old scores
					d = d < epsilon?epsilon:d;
						
					subs(0, base + i-1) = ind_arr[v][i]-1;
					subs(1, base + i-1) = v;
					vv(base + i-1) = d;			
				}				
			}
		}
		samples.clear();
		delete tree;
			
		field<mat> output(2);
		output(0) = idx;
		output(1) = dist;
		
		return(output);
	}
}






