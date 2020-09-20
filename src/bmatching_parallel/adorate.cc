#include "utils.h"
#include "graph.h"
#include "algo.h"


int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << "usage: "<<argv[0]<<" thread-count inputfile b-limit"<< std::endl;
        return 1;
    }

    int thread_count = std::stoi(argv[1]);
    int b_limit = std::stoi(argv[3]);
    std::string input_filename{argv[2]};
    std::ifstream infile(input_filename);
    graph current_graph;
    std::set<int> nodes;
    construct_graph_from_file(infile, current_graph, nodes);

    for (int b_method = 0; b_method < b_limit + 1; b_method++) {
        int weights = compute_for_method(current_graph, nodes, b_method, thread_count);
        std::cout << weights << std::endl;
    }
}
