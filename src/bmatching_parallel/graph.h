#ifndef MIMUW_ADORATE_42_GRAPH_H
#define MIMUW_ADORATE_42_GRAPH_H

#include "constants.h"
#include "constants.h"

#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <mutex>


struct edge {
    int to;
    int weight;
    bool operator < (const edge& rhs) const {
        return (weight >= rhs.weight);
    }
    inline bool bigger_than(const edge& rhs) const {
//        if (weight == rhs.weight)
//            return to > rhs.to;
        return (weight > rhs.weight);
    }
};
typedef std::map<int, std::vector<edge>> graph;
void add_edge(int node1, int node2, int weight, graph &graph);
void construct_graph_from_file(std::ifstream &infile, graph &graph, std::set<int> &nodes);

void construct_graph_from_file(std::ifstream& infile, graph& current_graph, std::set<int>& nodes) {
    int node1, node2, weight;
    std::string line;
    while (getline(infile, line)) {
        if (line[0] == '#') { continue; }
        std::istringstream iss(line);
        if (!(iss >> node1 >> node2 >> weight)) {
            std::cout << "ERROR" << std::endl;
            break;
        }
        nodes.insert(node1);
        nodes.insert(node2);
        add_edge(node1, node2, weight, current_graph);
        add_edge(node2, node1, weight, current_graph);

    }
}

void add_edge(int node1, int node2, int weight, graph& graph) {
    edge current_edge{node2, weight};
    graph::iterator it = graph.find(node1);
    if (it == graph.end()) { // node1 doesnt exist in graph
        std::vector<edge> adjList;
        adjList.push_back(current_edge);
        graph[node1] = adjList;
    } else { // node1 exists
        graph[node1].push_back(current_edge);
    }
}

#endif //MIMUW_ADORATE_42_GRAPH_H
