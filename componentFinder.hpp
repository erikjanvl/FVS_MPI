//
//  ComponentFinder.hpp
//  back
//
//  Created by Erik Jan van Leeuwen on 29/7/16.
//  Copyright © 2016 Erik Jan van Leeuwen. All rights reserved.
//

#ifndef ComponentFinder_h
#define ComponentFinder_h

#include "newgraph.hpp"
#include "inplace.hpp"
#include <vector>
#include <limits>
#include <iostream>

class ComponentFinder {
    std::vector<unsigned int> visited;
    std::vector<unsigned int> visitTypeSupport;
    inplace_stack<Graph::Node> parents;
    unsigned int currentVisit;
    Graph::Node currentStartNode;
    
    enum VisitType {
        NOT_IMPORTANT = 0,
        ONLY = 1,
        EXCLUDED = 2
    };
    VisitType vt;
    
    // this starts the next dfs visit and make sure integers do not overflow
    void nextVisit() {
        unsigned int maxVisit = std::numeric_limits<unsigned int>::max();
        if( currentVisit + 2 == maxVisit ) {
            currentVisit = 0;
            std::fill( visited.begin(), visited.end(), 0 );
        } else {
            currentVisit++;
        }
    }
    
    // this checks whether we are allowed to visit a particular node
    inline bool canVisitNode( const Graph::Node& n ) const {
        if( vt == VisitType::EXCLUDED ) {
            if( visited[n.id] == currentVisit+1 ) {
                return false;
            } else {
                return true;
            }
        }
        else if( vt == VisitType::ONLY ) {
            if( visited[n.id] == currentVisit+1 || visited[n.id] == currentVisit ) {
                // a node that was visited this round or is allowed to be visited anyway
                return true;
            } else {
                return false;
            }
        } else { // if vt == VisitType::NOT_IMPORTANT
            return true;
        }
    }
    
    // this checks whether a particular node has already been visited
    inline bool hasBeenVisited( const Graph::Node& n ) const {
        if( vt == VisitType::EXCLUDED || vt == VisitType::NOT_IMPORTANT ) {
            if( visited[n.id] == currentVisit ) {
                // clearly, this node was visited
                return true;
            } else {
                // might be that visited[n.id] == currentVisit+1, i.e. node is excluded, but regardless it wasn't visited
                return false;
            }
        } else { // if vt == VisitType::ONLY
            if( visited[n.id] == currentVisit+1 ) {
                // node is allowed to be visited, but clearly wasn't visited yet
                return false;
            } else if( visited[n.id] == currentVisit ) {
                // clearly visited
                return true;
            } else {
                // node is not allowed to visited or wasn't visited -> not visited
                return false;
            }
        }
    }

    
public:
    ComponentFinder(int c) : visited(c), parents(c), currentVisit(0), visitTypeSupport(c), vt(VisitType::NOT_IMPORTANT) {}
    
    // checks if the given graph induced by the vertices in only is a forest
    template <class Container>
    bool isForestOnly( const Graph& g, const Container& only ) {
        nextVisit();
        /*for(std::vector<Graph::Node>::const_iterator n = only.cbegin(); n != only.cend(); n++) {
            visited[n->id] = currentVisit+1;
        }*/
        for( int i = 0; i < only.size(); i++ ) {
            visited[only.at(i).id] = currentVisit+1;
        }
        vt = VisitType::ONLY;
        
        bool ret = isForestHelper(g);
        
        // we must do this because we used the +1
        nextVisit();
        
        return ret;
    }
    
    // checks if the given graph minus the vertices in exclude is a forest
    template <class Container>
    bool isForestExclude( const Graph& g, const Container& exclude ) {
        int excludeCount = exclude.size();
        if( g.edgeCountUnique() - (excludeCount * excludeCount/2) - (excludeCount*(g.nodeCount() - excludeCount)) >= g.nodeCount() - excludeCount ) {
            return false;
        }
        nextVisit();
        for( int i = 0; i < exclude.size(); i++ ) {
            visited[exclude.at(i).id] = currentVisit+1;
        }
        vt = VisitType::EXCLUDED;
        
        bool ret = isForestHelper(g);
        
        // we must do this because we used the +1 trick
        nextVisit();
        
        return ret;
    }
    
    // checks if the given graph is a forest
    bool isForest( const Graph& g ) {
        if( g.edgeCount() >= g.nodeCount() ) {
            return false;
        }
        Graph::Edge e;
        g.firstDouble(e);
        if( e.is_valid() ) {
            return false;
        }
        nextVisit();
        vt = VisitType::NOT_IMPORTANT;
        return isForestHelper(g);
    }
    
protected:
    
    bool isForestHelper( const Graph& g ) {
        for( Graph::NodeIt node(g); node != lemon::INVALID; ++node ) {
            if( canVisitNode(node) && !hasBeenVisited(node) ) {
                // start the visit from the current node
                if(!isForestHelper(g,node,node)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    bool isForestHelper( const Graph& g, const Graph::Node& node, const Graph::Node& parent ) {
        visited[node.id] = currentVisit;
        //for( Graph::NeighborIt opp(g, node); opp != lemon::INVALID; ++opp ) {
        for( Graph::IncUniqueEdgeIt e(g, node); e != lemon::INVALID; ++e ) {
            Graph::Node opp = g.oppositeNode(node,e);
            if( opp != parent && canVisitNode(opp) ) {
                if( g.isAdjacentDouble(node, opp) ) {
                    return false;
                } else if( hasBeenVisited(opp) ) {
                    return false;
                } else {
                    if(!isForestHelper(g,opp,node) ) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    
public:
    void initComponentFinder( Graph& g ) {
        nextVisit();
        currentStartNode = Graph::NodeIt(g);
    }
    
    void firstFromComponent( Graph& g, Graph::Node& n ) {
        parents.clear();
        for( Graph::NodeIt currentNodeit(g,currentStartNode) ; currentNodeit != lemon::INVALID; ++currentNodeit ) {
            currentStartNode = currentNodeit;
            if( visited[currentStartNode.id] != currentVisit ) {
                parents.push( currentStartNode );
                n = currentNodeit;
                return;
            }
        }
        n.invalidate();
    }
    
    void nextFromComponent( Graph& g, Graph::Node& n ) {
        while( !parents.empty() ) {
            Graph::Node current = parents.top();
            for( Graph::IncUniqueEdgeIt edge(g, current); edge != lemon::INVALID; ++edge ) {
                Graph::Node opp = g.oppositeNode(current,edge);
                if( visited[opp.id] != currentVisit ) {
                    parents.push( opp );
                    n = opp;
                    return;
                }
            }
            // visited all nodes of the current guy, continue looking at parent
            parents.pop();
        }
        
        n.invalidate();
    }
};

#endif /* ComponentFinder_h */
