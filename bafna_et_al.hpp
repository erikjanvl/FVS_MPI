#include <math.h>
#include <map>
#pragma once

#ifdef NEWCEIL
#define REALCEIL ceil
#else
#define REALCEIL std::ceil
#endif

#include "newgraph.hpp"
#include "reductions.hpp"
#include "fvs_optimization.hpp"

class BafnaEtAl {
public:
    std::tuple<int, std::list<Graph::Node>> approximate_fvs(Graph& graph, FvsOptimizer& optimizer) {
        std::map<Graph::Node, float> w;
        for (Graph::NodeIt node(graph); node != lemon::INVALID; ++node) {
            w[node] = 1;
        }
        return approximate_fvs(graph, optimizer, w);
    }

    std::tuple<int, std::list<Graph::Node>> approximate_fvs(Graph& graph, FvsOptimizer& optimizer, std::map<Graph::Node, float>& w) {
        double lower_bound = 0;
        std::list<Graph::Node> F;
        graph.disableMark();    // make sure we have a marker
        remove_nodes_degree_0_1(graph, true);
        while (graph.nodeCount() != 0) {
            bool has_degree_two_neighbours;
            int u_id, v_id;
            std::tie(has_degree_two_neighbours, u_id, v_id) = has_get_degree_two_neighbours(graph);
            if (has_degree_two_neighbours) {
                Graph::Node u = graph.nodeFromId(u_id), v = graph.nodeFromId(v_id);
                int u__id, v__id;
                std::tie(u__id, v__id) = get_next_neighbours(graph, u, v);
                if (u__id != v__id) {
                    graph.contractTwoStore(v);
                    w[u] = std::min(w[u], w[v]);
                } else {
                    Graph::Node u_ = graph.nodeFromId(u__id);
                    Graph::Node min_node = u;
                    for (Graph::Node n : (Graph::Node[]){v, u_}) if (w[n] < w[min_node]) min_node = n;
                    double gamma = w[min_node];
                    for (Graph::Node n : (Graph::Node[]){v, u, u_}) w[n] -= gamma;
                    lower_bound += gamma;
                    F.push_back(min_node);
                    graph.disableStore(min_node);
                    std::queue<Graph::Node> min_nodes_neighbours;
                    for (Graph::IncEdgeIt e(graph, min_node); e != lemon::INVALID; ++e) min_nodes_neighbours.push(graph.oppositeNode(min_node, e));
                    remove_nodes_degree_0_1(graph, min_nodes_neighbours, true);
                }
            } else {
                // Compute gamma from smallest node
                Graph::Node node0 = Graph::NodeIt(graph);
                float gamma = w[node0] / (graph.degree(node0) - 1);
                for (Graph::NodeIt node(graph); node != lemon::INVALID; ++node) {
                    float val = w[node] / (float) (graph.degree(node) - 1);
                    if (val < gamma) {
                        gamma = val;
                    }
                }
                // Subtract gamma from all nodes
                for (Graph::NodeIt node(graph); node != lemon::INVALID; ++node) {
                    w[node] -= gamma * (graph.degree(node) - 1);
                }
                lower_bound += gamma * (graph.edgeCount() - graph.nodeCount() + 1);

                // Remove nodes with w == 0
                std::vector<Graph::Node> remove_nodes;
                for (Graph::NodeIt node(graph); node != lemon::INVALID; ++node) {
                    if (w[node] <= 0) {
                        graph.disableStore(node);
                        F.push_back(node);
                    }
                }
                remove_nodes_degree_0_1(graph, true);
            }
        }
        //  Undo changes to graph up to the marker
        graph.restoreDisabled(lemon::INVALID,true);
        
        optimizer.optimizeList(graph, F);

        return std::make_pair(REALCEIL(lower_bound), F);
    }

    std::tuple<bool, int, int> has_get_degree_two_neighbours(Graph &graph) {
//        for (Graph::EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
//            Graph::Node u = graph.u(edge);
//            Graph::Node v = graph.v(edge);
//            if (graph.degree(u) == 2 && graph.degreeUnique(u) == 2 && graph.degree(v) == 2 && graph.degreeUnique(v) == 2) {
//                return std::make_tuple(true, graph.id(u), graph.id(v));
//            }
//        }
        for( Graph::NodeIt node(graph); node != lemon::INVALID; ++node ) {
            if( graph.degree(node) == 2 && graph.degreeUnique(node) == 2 ) {
                for( Graph::IncEdgeIt edge(graph,node); edge != lemon::INVALID; ++edge ) {
                    Graph::Node opp = graph.oppositeNode(node, edge);
                    if( graph.degree(opp) == 2 && graph.degreeUnique(opp) == 2 ) {
                        return std::make_tuple(true, graph.id(node), graph.id(opp));
                    }
                }
            }
        }
        return std::make_tuple(false, -1, -1);
    }

    std::tuple<int, int> get_next_neighbours(Graph& graph, Graph::Node& u, Graph::Node& v) {
        Graph::Node u_, v_;
        for (Graph::IncEdgeIt e(graph, u); e != lemon::INVALID; ++e) {
            Graph::Node opposite = graph.oppositeNode(u, e);
            if (opposite != v) u_ = opposite;
        }
        for (Graph::IncEdgeIt e(graph, v); e != lemon::INVALID; ++e) {
            Graph::Node opposite = graph.oppositeNode(v, e);
            if (opposite != u) v_ = opposite;
        }
        return std::make_tuple(graph.id(u_), graph.id(v_));
    }
};