#ifndef LEMON_ONLY_TEMPLATES
#define LEMON_ONLY_TEMPLATES
#endif

#include <assert.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>

#include <lemon/core.h>
#include <lemon/error.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/bits/map_extender.h>
#include <lemon/bits/default_map.h>
#include <ctime>

#include "newgraph.hpp"
#include "inplace.hpp"
#include "becker_geiger.hpp"
#include "componentFinder.hpp"
#include "fvs_optimization.hpp"
#include "reductions.hpp"
#include "bafna_et_al.hpp"
#include "branch_algo.hpp"


void split_string_on_white_space(const std::string &s, char delim,
                                 std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
  //while (std::getline(ss, item, delim)) {
  //  elems.push_back(item);
  //}
    while( ss >> item ) {
        elems.push_back(item);
    }
}

void read_lines(
    std::set<std::pair<std::string, std::string>> &edge_string_pairs,
    std::set<std::string> &node_strings) {
    std::string line;
    while (std::getline(std::cin, line)) {
        std::string prefix = "#";
        if (std::equal(prefix.begin(), prefix.end(), line.begin())){
            continue;
        }
        std::vector<std::string> two_nodes;
        split_string_on_white_space(line, ' ', two_nodes);
      
        node_strings.insert(two_nodes[0]);
        node_strings.insert(two_nodes[1]);
        edge_string_pairs.insert(std::make_pair(two_nodes[0], two_nodes[1]));
    }
}

void fill_graph(
    Graph &G,
    const std::set<std::pair<std::string, std::string>> &edge_string_pairs,
    const std::set<std::string> &node_strings,
    std::map<Graph::Node, std::string> &node_to_string_map) {

  std::vector<Graph::Node> nodes = G.buildEnableNodes();
  assert(nodes.size() == node_strings.size());

  std::map<std::string, Graph::Node> string_to_node_map;
  int i = 0;
  for (auto str : node_strings) {
    string_to_node_map[str] = nodes[i];
    node_to_string_map[nodes[i]] = str;
    i++;
  }

  for (auto edge_pair : edge_string_pairs) {
    G.increaseAdjacency(string_to_node_map[edge_pair.first],
                        string_to_node_map[edge_pair.second]);
  }
}

template <class Container>
void write_removed_nodes(const Container& removed_nodes, std::map<Graph::Node, std::string>& node_to_string_map) {
    for( int i = 0; i < removed_nodes.size(); i++) {
        std::cout << node_to_string_map[removed_nodes.at(i)] << std::endl;
        //std::cout << removed_nodes.at(i).id << std::endl;
    }
}

template <>
void write_removed_nodes(const std::list<Graph::Node>& removed_nodes, std::map<Graph::Node, std::string>& node_to_string_map) {
    for (auto node : removed_nodes) {
        std::cout << node_to_string_map[node] << std::endl;
    }
}


int main() {
    time_t start;
    time(&start);
    std::set<std::pair<std::string, std::string>> edge_string_pairs,edge_string_pair;
    std::set<std::string> node_strings;

    read_lines(edge_string_pairs, node_strings);

    Graph G(node_strings.size());
    Graph G1(node_strings.size());
    std::map<Graph::Node, std::string> node_to_string_map,node_to_string_map1;

    fill_graph(G, edge_string_pairs, node_strings, node_to_string_map);
    fill_graph(G1, edge_string_pairs, node_strings, node_to_string_map1);
    
    ComponentFinder c( G.nodeCountOriginal() );
    FvsOptimizer opt(G.nodeCountOriginal());


//    int i = 1;
//    for (Graph::EdgeIt edge(G); edge != lemon::Invalid(); ++edge) {
//        Graph::Node u = G.u(edge);
//        Graph::Node v = G.v(edge);
//        std::cout << i << " " << node_to_string_map[u] << ", " << node_to_string_map[v] << std::endl;
//        i++;
//    }
//    i = 1;
//    Graph::Node n = Graph::NodeIt();
//    for (Graph::NeighborIt node(G,n); node != lemon::Invalid(); ++node) {
//        Graph::Node v = node;
//        std::cout << i << " " << node_to_string_map[n] << ", " << node_to_string_map[v] << std::endl;
//        i++;
//    }
//    int i = 1;
//    for( Graph::NodeIt node(G); node != lemon::Invalid(); ++node ) {
//        for (Graph::IncEdgeIt edge(G,node); edge != lemon::Invalid(); ++edge) {
//            Graph::Node u = node;
//            Graph::Node v = G.oppositeNode(node, edge);
//            std::cout << i << " " << u.id << ", " << v.id << std::endl;
//            i++;
//        }
//    }

    inplace_arraylist<Graph::Node> removed_nodes( G.nodeCountOriginal() );
    
    //std::cerr << "Starting Pre " << G.nodeCount() << " " << G.edgeCount() << std::endl;
    initial_kernelize(G, removed_nodes);
    //std::cerr << "Mid Pre " << G.nodeCount() << " " << G.edgeCount() << std::endl;
    //while( remove_bridges(G) ) ;
    //std::cerr << "End Pre " << G.nodeCount() << " " << G.edgeCount() << std::endl;
    //std::cerr << removed_nodes.size() << std::endl;
    
    inplace_arraylist<Graph::Node> bg_nodes( G.nodeCountOriginal() );
    BeckerGeiger bg( G.nodeCountOriginal() );
    bg.init(G);
    bg.execute(G, opt, bg_nodes, false);
    
    int bafnaLB;
    std::list<Graph::Node> bafnaSOL;
    std::tie(bafnaLB, bafnaSOL) = BafnaEtAl().approximate_fvs(G, opt);
    
    if( bafnaLB == bg_nodes.size() ) {
        for( int i = 0; i < bg_nodes.size(); i++ ) {
            removed_nodes.push_back(bg_nodes.at(i));
        }
    } else if( bafnaLB == bafnaSOL.size() ) {
        for(std::list<Graph::Node>::const_iterator n = bafnaSOL.cbegin(); n != bafnaSOL.cend(); n++) {
            removed_nodes.push_back(static_cast<Graph::Node>(*n));
        }
    }
    else {
        //std::cerr << "Starting Branch" << std::endl;
        RecursiveFVS rec( G.nodeCountOriginal() );
        //std::cerr << "Init Branch" << std::endl;
        rec.init(G);
        //std::cerr << "Run Branch" << std::endl;
        try {
            rec.execute(G, c, opt, start);
        } catch(...) {}
        inplace_arraylist<Graph::Node>& sol = rec.getReturn();
        for( int i = 0; i < sol.size(); i++ ) {
            removed_nodes.push_back(sol.at(i));
        }
        //std::cerr << "After Branch" << std::endl;
    }
//    else {
//        int upp = bg_nodes.size();
//        if( bafnaSOL.size() < upp ) {
//            upp = bafnaSOL.size();
//        }
//        std::list<Graph::Node> beckerSOL = BeckerEtAl().randomized_fvs_once(G,100,upp);
//        std::cerr << beckerSOL.size() << std::endl;
//        for(std::list<Graph::Node>::const_iterator n = beckerSOL.cbegin(); n != beckerSOL.cend(); n++) {
//            removed_nodes.push_back(static_cast<Graph::Node>(*n));
//        }
//    }
    
//    i = 1;
//    for( Graph::NodeIt node(G); node != lemon::Invalid(); ++node ) {
//        for (Graph::IncEdgeIt edge(G,node); edge != lemon::Invalid(); ++edge) {
//            Graph::Node u = node;
//            Graph::Node v = G.oppositeNode(node, edge);
//            std::cout << i << " " << u.id << ", " << v.id << std::endl;
//            i++;
//        }
//    }


    //std::list<Graph::Node> removed_nodes;
    // int lower_bound;
    // std::tie(lower_bound, removed_nodes) = BafnaEtAl().approximate_fvs(graph);
    // std::cout << "lower bound: " << lower_bound << std::endl;
    // Graph g;
    // read_graph(g);
    // FvsOptimizer().optimize(graph, removed_nodes);
    
    //std::cerr << "Init opt" << std::endl;
    //opt.optimize(G, removed_nodes);
    //std::cerr << "After Opt" << std::endl;
    
    write_removed_nodes(removed_nodes, node_to_string_map);
    //std::cerr << "Size: " << removed_nodes.size() << std::endl;
    
    if( !c.isForestExclude(G1,removed_nodes) ) {
        ;//std::cerr << "WRONGa " << G.nodeCount() << " " << G.edgeCountUnique() << std::endl;
    }
    
    for( int i = 0; i < removed_nodes.size(); i++) {
        G.disableStore(removed_nodes.at(i));
    }
    if( !c.isForest(G) ) {
        ;//std::cerr << "WRONGb " << G.nodeCount() << " " << G.edgeCountUnique() << std::endl;
    }
    G.restoreDisabled();
    
    return 0;
}
