//
//  becker_geiger.h
//  back
//
//  Created by Erik Jan van Leeuwen on 29/7/16.
//  Copyright © 2016 Erik Jan van Leeuwen. All rights reserved.
//

#ifndef becker_geiger_h
#define becker_geiger_h

#include <iostream>
#include <map>
#include "newgraph.hpp"
#include "inplace.hpp"
#include "fvs_optimization.hpp"

class BeckerGeiger
{
    inplace_fifo<Graph::Node> deleteQueue;
    std::vector<double> weights;
    
    inline double computeGamma( const Graph& g, const Graph::Node& n ) const {
        return weights[n.id] / g.degree(n);
    }
    
public:
    BeckerGeiger(int c) : deleteQueue(2*c), weights(c) {}

    void init( const Graph& g ) {
        for( Graph::NodeIt node(g); node != lemon::INVALID; ++node ) {
            weights[static_cast<Graph::Node>(node).id] = 1.0;
        }
    }
    
    void init( const Graph& g, const Graph::Node& n, double initialWeight ) {
        init(g);
        weights[n.id] = initialWeight;
    }
    
    template <class ResultContainer>
    void execute( Graph& g, FvsOptimizer& optimizer, ResultContainer& result, bool clearResult = false ) {
        if( clearResult ) {
            result.clear();
        }
        g.disableMark();
        
        while( g.nodeCount() > 0 ) {
            deleteQueue.clear();
            double minGamma = computeGamma(g, g.firstNode());
            
            // compute the element of minimum gamma
            for( Graph::NodeIt node(g); node != lemon::INVALID; ++node ) {
                double gamma = computeGamma(g, node);
                if( gamma < minGamma ) {
                    minGamma = gamma;
                    deleteQueue.clear();
                    deleteQueue.push( node );
                } else if( gamma == minGamma ) {
                    deleteQueue.push( node );
                }
            }
            
            // add the vertices of minimum gamma to the solution
            for( int i = 0; i < deleteQueue.size(); i++ ) {
                result.push_back( deleteQueue.at(i) );
                g.disableStore( deleteQueue.at(i) );
            }
            
            deleteQueue.clear();
            bool didSomething = true;
            while( didSomething ) {
                didSomething = false;
                for( Graph::NodeIt node(g); node != lemon::INVALID; ++node ) {
                    if( g.degree(node) == 1 ) {
                        Graph::Node nn = Graph::NeighborIt(g,node);
                        weights[nn.id] -= minGamma;
                    }
                    if( g.degree(node) <= 1 ) {
                        g.disableStore(node);
                    }
                }
            }
            
//            // now iterate over the vertices on the deleteQueue and remove them
//            // also make sure that if they have vertices that get degree <= 1, then we mark them for deletion as well
//            while( deleteQueue.size() > 0 ) {
//                Graph::Node node = deleteQueue.front();
//                if( g.degree(node) == 1 ) {
//                    //Guse.node[nn]['w'] = Guse.node[nn]['w'] - mind
//                    Graph::Node nn = Graph::NeighborIt(g,node);
//                    weights[nn.id] -= minGamma;
//                }
//                for( Graph::IncSingleEdgeIt edge(g, node); edge != lemon::INVALID; ++edge ) {
//                    Graph::Node opp = g.oppositeNode(node, edge);
//                    if( g.degree(opp) - 1 <= 1 ) {
//                        deleteQueue.push(opp);
//                    }
//                }
//                for( Graph::IncDoubleEdgeIt edge(g, node); edge != lemon::INVALID; ++edge ) {
//                    Graph::Node opp = g.oppositeNode(node, edge);
//                    if( g.degree(opp) - 2 <= 1 ) {
//                        deleteQueue.push(opp);
//                    }
//                }
//                g.disableStore(node);
//                deleteQueue.pop();
//            }
        }
        
        // restore the graph to its original happiness
        g.restoreDisabled(lemon::INVALID, true);
        
        optimizer.optimize( g, result );
    }
};



#endif /* becker_geiger_h */
