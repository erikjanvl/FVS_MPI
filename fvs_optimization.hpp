// #include "graph_algorithm.hpp"

#ifndef FEEDBACKVERTEXSET_FVS_OPTIMIZATION_H
#define FEEDBACKVERTEXSET_FVS_OPTIMIZATION_H

#include <vector>

#include "newgraph.hpp"
#include "uf.hpp"
#include "componentFinder.hpp"

class FvsOptimizer {
    UnionFind<Graph::Node> components;
    std::vector<unsigned int> seen;
    unsigned int currentSeen = 0;
public:
    FvsOptimizer(int c) : seen(c,0), components(c) {}
    
    template <class Container>
    void optimize( Graph& g, Container& fvs ) {
        components.clear();
        currentSeen++;
        unsigned int baseSeen = currentSeen;
        
        
        // initialize the union find datastructure
        g.disableMark();
        for( int i = 0; i < fvs.size(); i++ ) {
            if( g.isEnabled(fvs.at(i)) ) {
                g.disableStore(fvs.at(i));
            }
            seen.at(fvs.at(i).id) = baseSeen;  // this ensures that the node is not visible later
        }
        for(Graph::EdgeIt e(g); e != lemon::INVALID; ++e) {
            components.join(g.u(e), g.v(e));
        }
        g.restoreDisabled(lemon::INVALID, true);
        
        for( int i = fvs.size() - 1; i >= 0; i-- ) {
            if( !g.isEnabled(fvs.at(i)) ) {
                continue;
            }
            currentSeen++;
            Graph::Node node = fvs.at(i);
            bool necessary = false;
            bool noneighbors = true;
            
            // look for two neighbors in the same component of G - fvs
            for( Graph::NeighborIt nn(g, node); nn != lemon::INVALID; ++nn) {
                if( seen.at(nn.id) != baseSeen ) {
                    // seen != 0 -> nn is not in the FVS
                    noneighbors = false;
                    if( seen.at(components.find(nn)) == currentSeen || g.isAdjacentDouble(node, nn) ) {
                        // we have seen the component that contains nn before this iteration
                        // means that node has two neighbors in that component
                        necessary = true;
                        break;
                    } else {
                        seen.at(components.find(nn)) = currentSeen;
                    }
                }
            }
            
            if( !necessary || noneighbors ) {
                // node is not necessary -> remove from list and mark as visible
                fvs.clear(i);
                seen.at( node.id ) = baseSeen+1;
                
                //now merge the node with the rest
                for( Graph::NeighborIt nn(g, node); nn != lemon::INVALID; ++nn) {
                    if( seen.at(nn.id) != baseSeen ) {
                        // seen != 0 -> nn is not in the FVS
                        components.join(node,nn);
                    }
                }
            }
        }
        
    }
    
    
    void optimizeList( Graph& g, std::list<Graph::Node>& fvs ) {
        components.clear();
        currentSeen++;
        unsigned int baseSeen = currentSeen;
        
        
        // initialize the union find datastructure
        g.disableMark();
        for(std::list<Graph::Node>::const_iterator itn = fvs.cbegin(); itn != fvs.cend(); itn++) {
            Graph::Node node = static_cast<Graph::Node>(*itn);
            
            if( g.isEnabled(node) ) {
                g.disableStore(node);
            }
            seen.at(node.id) = baseSeen;  // this ensures that the node is not visible later
        }
        for(Graph::EdgeIt e(g); e != lemon::INVALID; ++e) {
            components.join(g.u(e), g.v(e));
        }
        g.restoreDisabled(lemon::INVALID, true);
        
        
        for(std::list<Graph::Node>::const_reverse_iterator itn = fvs.crbegin(); itn != fvs.crend(); itn++) {
            Graph::Node node = static_cast<Graph::Node>(*itn);
            
            if( !g.isEnabled(node) ) {
                continue;
            }
            currentSeen++;
            //Graph::Node node = fvs.at(i);
            bool necessary = false;
            bool noneighbors = true;
            
            // look for two neighbors in the same component of G - fvs
            for( Graph::NeighborIt nn(g, node); nn != lemon::INVALID; ++nn) {
                if( seen.at(nn.id) != baseSeen ) {
                    // seen != 0 -> nn is not in the FVS
                    noneighbors = false;
                    if( seen.at(components.find(nn)) == currentSeen || g.isAdjacentDouble(node, nn) ) {
                        // we have seen the component that contains nn before this iteration
                        // means that node has two neighbors in that component
                        necessary = true;
                        break;
                    } else {
                        seen.at(components.find(nn)) = currentSeen;
                    }
                }
            }
            
            if( !necessary || noneighbors ) {
                // node is not necessary -> remove from list and mark as visible
                fvs.remove(node);
                seen.at( node.id ) = baseSeen+1;
                
                //now merge the node with the rest
                for( Graph::NeighborIt nn(g, node); nn != lemon::INVALID; ++nn) {
                    if( seen.at(nn.id) != baseSeen ) {
                        // seen != 0 -> nn is not in the FVS
                        components.join(node,nn);
                    }
                }
            }
        }
        
    }
    
    
//    void optimize(BaseGraph& graph, std::list<BaseGraph::Node>& fvs) {
//        std::vector<BaseGraph::Node> remove;
//        BaseGraph graph_wo;
//        BaseGraph::NodeMap<BaseGraph::Node> nm(graph);
//        BaseGraph::EdgeMap<BaseGraph::Edge> em(graph);
//        lemon::graphCopy(graph, graph_wo).nodeRef(nm).edgeCrossRef(em).run();
//        for (BaseGraph::Node& node : fvs) {
//            graph_wo.erase(nm[node]);
//        }
//        std::list<BaseGraph::Node>::iterator i = fvs.begin();
//        while (i != fvs.end()) {
//            BaseGraph::Node node = *i;
//            BaseGraph::Node new_node = graph_wo.addNode();
//            nm.set(node, new_node);
//            for (BaseGraph::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
//                graph_wo.addEdge(new_node, nm[graph.oppositeNode(node, e)]);
//            }
//            if (!graph_wo.is_forrest()) {
//                graph_wo.erase(new_node);
//                ++i;
//            } else {
//                fvs.erase(i++);
//            }
//        }
//    }
};

#endif //FEEDBACKVERTEXSET_FVS_OPTIMIZATION_H
