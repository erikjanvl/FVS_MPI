#include <queue>

#include "newgraph.hpp"
#include "reductions.hpp"
#include <lemon/core.h>
#include <iostream>

bool remove_nodes_selfloops(Graph& g, inplace_arraylist<Graph::Node>& res) {
    bool didSomething = false;
    for( Graph::NodeIt n(g); n != lemon::INVALID; ++n ) {
        if( g.hasSelfLoop(n) ) {
            didSomething = true;
            res.push_back(n);
            g.disable(n);
        }
    }
    
    return didSomething;
}

bool remove_nodes_degree_0_1(Graph& graph, std::queue<Graph::Node>& potential_nodes, bool store) {
    bool didSomething = false;
    while(!potential_nodes.empty()) {
        Graph::Node node = potential_nodes.front();
        potential_nodes.pop();
        if( !graph.isEnabled(node) )
            continue;
        if (graph.degree(node) == 0) {
            didSomething = true;
            if( store ) {
                graph.disableStore(node);
            } else {
                graph.disable(node);
            }
        } else if (graph.degree(node) == 1) {
            didSomething = true;
            potential_nodes.push(graph.oppositeNode(node, Graph::IncEdgeIt(graph, node)));
            if( store ) {
                graph.disableStore(node);
            } else {
                graph.disable(node);
            }
//        } else if (graph.degree(node) == 2 && graph.degreeUnique(node) == 2) {
//            Graph::Node u, v;
//            bool wasAlreadyDouble;
////            Graph::Node u = Graph::NeighborIt(graph,node);
////            Graph::NeighborIt v(graph,node);
////            ++v;
////            graph.disable(node);
////            graph.increaseAdjacency(u,v);
////            std::tie(u,v,wasAlreadyDouble) = graph.contractTwo(node);
////            potential_nodes.push(u);
////            potential_nodes.push(v);
//            //if( u == v ) throw 0;
        }
    }
    
    return didSomething;
}

bool remove_nodes_degree_0_1(Graph& graph, bool store ) {
    std::queue<Graph::Node> potential_nodes;
    if (potential_nodes.empty()) {
        for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
            if( graph.degree(n) <= 1 ) {
                potential_nodes.push(n);
            }
        }
    }
    return remove_nodes_degree_0_1(graph, potential_nodes, store);
}

bool remove_degree_two_two_neighbors(Graph& graph, bool store) {
    bool changed = false;
    std::list<Graph::Node> bridges;
    for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
        if (graph.degree(n) == 2 && graph.degreeUnique(n) == 2 ) {
            bridges.push_back(n);
        }
    }
    for (Graph::Node node : bridges) {
        if( store ) {
            graph.contractTwoStore(node);
        } else {
            graph.contractTwo(node);
        }
    }
    if (bridges.size() > 0) {
        changed = true;
    }
    return changed;
}

//bool remove_leaves_nodes_wo_edges(Graph& graph) {
//    bool changed_overall = false;
//    bool changed = false;
//    do {
//        changed = false;
//        std::list<Graph::Node> remove_nodes;
//        for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
//            if (graph.degree(n) <= 1) {
//                remove_nodes.push_back(n);
//            }
//        }
//        for (Graph::Node node : remove_nodes) {
//            graph.disableStore(node);
//        }
//        if (remove_nodes.size() > 0) {
//            changed = true;
//            changed_overall = true;
//        }
//    } while (changed);
//    return changed_overall;
//}

//template <class ResultContainer>
bool remove_degree_two_one_neighbor(Graph& graph, inplace_arraylist<Graph::Node>& res, bool store  )
{
    bool changed = false;
    std::list<Graph::Node> nodes;
    for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
        if(graph.degree(n) == 2 && graph.degreeUnique(n) == 1 ) {
            Graph::Node opp = Graph::NeighborIt(graph,n);
            nodes.push_back(opp);
            changed = true;
            break;
        }
    }
    for (Graph::Node node : nodes) {
        if( store ) {
            graph.disableStore(node);
        } else {
            graph.disable(node);
        }
        res.push_back(node);
    }
    return changed;
}

bool remove_degree_two_one_neighbor(Graph& graph, std::list<Graph::Node>& res, bool store  )
{
    bool changed = false;
    std::list<Graph::Node> nodes;
    for (Graph::NodeIt n(graph); n != lemon::INVALID; ++n) {
        if(graph.degree(n) == 2 && graph.degreeUnique(n) == 1 ) {
            Graph::Node opp = Graph::NeighborIt(graph,n);
            nodes.push_back(opp);
            changed = true;
            break;
        }
    }
    for (Graph::Node node : nodes) {
        if( store ) {
            graph.disableStore(node);
        } else {
            graph.disable(node);
        }
        res.push_back(node);
    }
    return changed;
}

bool initial_kernelize(Graph& g, inplace_arraylist<Graph::Node>& res) {
    remove_nodes_selfloops(g,res);
    
    bool didSomething = true;
    while( didSomething ) {
        didSomething = false;
        didSomething = remove_nodes_degree_0_1(g) || didSomething;
        while( remove_degree_two_two_neighbors(g) ) {
            didSomething = true;
        }
        didSomething = remove_degree_two_one_neighbor(g, res, false) || didSomething;
    }
    
    return didSomething;
}