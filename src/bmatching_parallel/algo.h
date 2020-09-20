#ifndef MIMUW_ADORATE_42_ALGO_H
#define MIMUW_ADORATE_42_ALGO_H

#include "graph.h"
#include "blimit.hpp"
#include "constants.h"
#include "thread_pool.h"

#include <mutex>
#include <queue>

struct Mutex {
    std::mutex mux;
};

int compute_for_method(graph& current_graph, std::set<int>& nodes, int b_method, int thread_count);
void find_suitors_for_node(const int current_node,
                           const int b_method,
                           std::map<int, std::vector<int>>& T,
                           std::vector<int>& R,
                           std::map<int, std::priority_queue<edge>>& S,
                           graph& current_graph,
                           std::mutex& mutex,
                           std::map<int, Mutex*>& mutexes);
bool is_eligible_candidate_node(std::map<int, std::priority_queue<edge>>& S,
                                edge& candidate_edge,
                                int current_node_b);
int compute_results(std::map<int, std::priority_queue<edge>>& S);


int compute_for_method(graph& current_graph, std::set<int>& nodes, int b_method, int thread_count) {
    std::map<int, std::vector<int>> T;
    std::vector<int> R;
    std::map<int, std::priority_queue<edge>> S;
    std::map<int, Mutex*> mutexes;

    for(int node : nodes) {
        S[node] = std::priority_queue<edge>();
        T[node] = std::vector<int>();
        mutexes[node] = new Mutex;
    }
    std::vector<int> Q(nodes.begin(), nodes.end());
    std::mutex mutex;
    while (!Q.empty()) {
        if (thread_count > 1 ) {
            ThreadPool pool(thread_count-1);
            for (int current_node : Q) {
                pool.enqueue([current_node, b_method, &T, &R, &S, &current_graph, &mutex, &mutexes]
                             {find_suitors_for_node(current_node, b_method, T, R, S, current_graph, mutex, mutexes);});
            }
        } else
        {
            for (int current_node : Q) {
                find_suitors_for_node(current_node, b_method, T, R, S, current_graph, mutex, mutexes);
            }
        }
        Q = R;
        R.clear();
    }
    return compute_results(S);
}


void find_suitors_for_node(const int current_node,
                           const int b_method,
                           std::map<int, std::vector<int>>& T,
                           std::vector<int>& R,
                           std::map<int, std::priority_queue<edge>>& S,
                           graph& current_graph,
                           std::mutex& mutex,
                           std::map<int, Mutex*>& mutexes) {
    int b = bvalue(b_method, current_node);
    while (T[current_node].size() < b) {
        // Find the best candidate
        bool found_candidate = false;
        edge current_best_candidate_edge{-1, 0};
        for (edge candidate_edge : current_graph[current_node]) {
            if (!(contains(T[current_node], candidate_edge.to))) {
                if (is_eligible_candidate_node(S, candidate_edge, bvalue(b_method, candidate_edge.to))) {
                    found_candidate = true;
                    if (candidate_edge.bigger_than(current_best_candidate_edge)) {
                        current_best_candidate_edge = candidate_edge;
                    }
                }
            }
        }
        if (!(found_candidate)) {
            break;
        }
        else { // current_node will adorate current_best_candidate_edge
            int candidate_node = current_best_candidate_edge.to;
            mutex.lock();
//            mutexes[candidate_node]->mux.lock();
//            mutexes[current_node]->mux.lock();

            int canddidate_b_value = bvalue(b_method, candidate_node);
            if (!is_eligible_candidate_node(S, current_best_candidate_edge, canddidate_b_value)) {
                mutex.unlock();
                continue;
            }
            if (S[candidate_node].size() == canddidate_b_value) {
                // Annuling
                edge annulled_edge = S[candidate_node].top();
                mutexes[annulled_edge.to]->mux.lock();
                S[candidate_node].pop();
                std::vector<int>::iterator position = find(T[annulled_edge.to].begin(),
                                                           T[annulled_edge.to].end(),
                                                           candidate_node);
                if (position != T[annulled_edge.to].end())
                    T[annulled_edge.to].erase(position);
                R.push_back(annulled_edge.to);
                mutexes[annulled_edge.to]->mux.unlock();
            }
            T[current_node].push_back(candidate_node);
            // Find corresponding edge in the opposite direction
            std::vector<edge> edges_of_candidate = current_graph[candidate_node];
            for (edge e : edges_of_candidate) {
                if (e.to == current_node) {
                    S[candidate_node].push(e);
                    break;
                }
            }
            mutex.unlock();
        }
    }
}

bool is_eligible_candidate_node(std::map<int, std::priority_queue<edge>>& S, edge& candidate_edge, int current_node_b) {
    bool candidate_has_suitors = !S[candidate_edge.to].empty();
    return ((candidate_has_suitors && candidate_edge.bigger_than(S[candidate_edge.to].top()))
            || (candidate_has_suitors && S[candidate_edge.to].size() < current_node_b)
            || !candidate_has_suitors);
}

int compute_results(std::map<int, std::priority_queue<edge>>& S) {
    std::set<int> already_printed;
    int counter = 0;
    for (std::pair<const int, std::priority_queue<edge>> & e : S) {
        if (!(contains(already_printed, e.first))) {
            already_printed.insert(e.first);

            while(!e.second.empty()) {
                edge top = e.second.top();
                if (!(contains(already_printed, top.to))) {
                    counter += top.weight;
                }

                e.second.pop();
            }
        }
    }
    return counter;
}

#endif //MIMUW_ADORATE_42_ALGO_H
