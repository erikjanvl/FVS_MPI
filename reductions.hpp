#ifndef FEEDBACKVERTEXSET_REDUCTIONS_H
#define FEEDBACKVERTEXSET_REDUCTIONS_H

#include <queue>
#include "inplace.hpp"
#include "newgraph.hpp"

bool remove_nodes_selfloops(Graph& g, inplace_arraylist<Graph::Node>& res);

bool remove_nodes_degree_0_1(Graph& graph, bool store = false);
bool remove_nodes_degree_0_1(Graph& graph, std::queue<Graph::Node>& potential_nodes, bool store = false);

bool remove_degree_two_two_neighbors(Graph& graph, bool store = false);
//bool remove_leaves_nodes_wo_edges(Graph& graph);

//template <class ResultContainer>
bool remove_degree_two_one_neighbor(Graph& graph, inplace_arraylist<Graph::Node>& res, bool store = false );
bool remove_degree_two_one_neighbor(Graph& graph, std::list<Graph::Node>& res, bool store = false );

bool initial_kernelize(Graph& graph, inplace_arraylist<Graph::Node>& res);
//bool kernelize(Graph& g, inplace_arraylist<Graph::Node>& F, std::vector<int>& weights, std::list<Graph::Node>& bafnaSOL);

#endif //FEEDBACKVERTEXSET_REDUCTIONS_H
