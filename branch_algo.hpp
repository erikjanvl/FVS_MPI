//
//  branch_algo.hpp
//  back
//
//  Created by Erik Jan van Leeuwen on 1/8/16.
//  Copyright © 2016 Erik Jan van Leeuwen. All rights reserved.
//

#ifndef branch_algo_h
#define branch_algo_h

#include "inplace.hpp"
#include "newgraph.hpp"
#include "bafna_et_al.hpp"
#include "componentFinder.hpp"
#include "reductions.hpp"
#include "fvs_optimization.hpp"
#include <iostream>
#include <limits>
#include <ctime>

class RecursiveFVS {
    inplace_arraylist<Graph::Node> bestSoFar;
    int bestSoFarSize = 0;
    std::map<Graph::Node, float> bafna_weights;
    std::vector<int> weights;
    inplace_arraylist<Graph::Node> F;
    
    int largeWeight;
    
    void updateBestSoFar() {
        bestSoFar.clear();
        for( int i = 0; i < F.size(); i++) {
            bestSoFar.push_back(F.at(i));
        }
        bestSoFarSize = bestSoFar.size();
//        for( int i = 0; i < bestSoFar.size(); i++ ) {
//            std::cerr << bestSoFar.at(i).id << " " << std::endl;
//        }
//        std::cerr << " best" << std::endl;
    }
    
    void updateBestSoFar(const std::list<Graph::Node>& bafnaSOL) {
        updateBestSoFar();
        for(std::list<Graph::Node>::const_iterator n = bafnaSOL.cbegin(); n != bafnaSOL.cend(); n++) {
            bestSoFar.push_back(static_cast<Graph::Node>(*n));
        }
        bestSoFarSize = bestSoFar.size();
//        for( int i = 0; i < bestSoFar.size(); i++ ) {
//            std::cerr << bestSoFar.at(i).id << " " << std::endl;
//        }
//        std::cerr << " best2" << std::endl;
    }
    
public:
    RecursiveFVS(int c) : bestSoFar(c), bafna_weights(), weights(c,1), F(c) {}
    
    void init( Graph& g ) {
        for( Graph::NodeIt node(g); node != lemon::Invalid(); ++node) {
            weights[node.id] = 1;
        }
        bestSoFar.clear();
        F.clear();
        largeWeight = std::numeric_limits<int>::max();
        bestSoFarSize = largeWeight;
    }
    
    inplace_arraylist<Graph::Node>& getReturn() {
        return bestSoFar;
    }
    
    void execute(Graph &g, ComponentFinder& cf, FvsOptimizer& optimizer, time_t start ) {
        time_t now;
        time(&now);
        if( difftime(now, start) >= 29*60 ) {
            throw 0;
        }
        
        
        if( cf.isForest(g) ) {
            if( F.size() < bestSoFarSize ) {
                updateBestSoFar();
            }
            return;
        }
        
        // set the weights for bafna
        for( Graph::NodeIt node(g); node != lemon::Invalid(); ++node) {
            bafna_weights[static_cast<Graph::Node>(node)] = weights[node.id];
        }
        
//        int i = 1;
//            for (Graph::EdgeIt edge(g); edge != lemon::Invalid(); ++edge) {
//                Graph::Node u = g.u(edge);
//                Graph::Node v = g.v(edge);
//                std::cerr << i << " " << u.id << ", " << v.id << std::endl;
//                i++;
//            }
//        for (Graph::NodeIt edge(g); edge != lemon::Invalid(); ++edge) {
//            std::cerr << edge.id <<std::endl;
//        }
        
        // call bafna
        int bafnaLB;
        std::list<Graph::Node> bafnaSOL;
        std::tie(bafnaLB, bafnaSOL) = BafnaEtAl().approximate_fvs(g, optimizer, bafna_weights);
        
        // check if we have a new best solution
        // or whether we can reject the current branch
        if( F.size() + bafnaLB >= bestSoFarSize ) {
            //std::cerr << "OUT" <<std::endl;
            return;
        } else if( F.size() + bafnaSOL.size() < bestSoFarSize ) {
            //std::cerr << "F+bafna" << std::endl;
            updateBestSoFar(bafnaSOL);
//            if( bestSoFarSize == 6 ) throw 0;
//            if( bestSoFarSize == 6 )
//            {
//                int i = 1;
//                            for (Graph::EdgeIt edge(g); edge != lemon::Invalid(); ++edge) {
//                                Graph::Node u = g.u(edge);
//                                Graph::Node v = g.v(edge);
//                                std::cout << i << " " << u.id << ", " << v.id << std::endl;
//                                i++;
//                            }
//            }
        }
        
        //std::cerr << "Now with bestsf " << bestSoFarSize << " and bafna " << bafnaSOL.size() << std::endl;
        
        // add the marker so that I know until where to restore
        g.disableMark();
        int Fsize = F.size();
        //std::cerr << "F ";
//        for( int i = 0; i < Fsize; i++ ) {
//            std::cerr << F.at(i).id << " ";
//        }
//        std::cerr << std::endl;
        
        // kernelize
        // ret == 0 -> no rule applied
        // ret == 1 -> recurse
        // ret == 2 -> no need to recurse, just return
        short ret = kernelize( g, cf, optimizer, bafnaSOL, start );
        if( ret == 1 ) {
            //std::cerr << "Kernel did" << std::endl;
            execute(g, cf, optimizer, start);
        } else if( ret == 0 ) {
            //   branch on v of highest degree
            Graph::Node branchNode;
            bool hasBranchNode = false;
            
            for( Graph::NodeIt itn(g); itn != lemon::Invalid(); ++itn) {
                Graph::Node node = static_cast<Graph::Node>(itn);
                if( weights[node.id] == 1 ) {
                    if( !hasBranchNode ) {
                        hasBranchNode = true;
                        branchNode = node;
                    } else if( g.degree(node) > g.degree(branchNode) ) {
                        branchNode = node;
                    }
                }
            }
            
            if( hasBranchNode ) {
                //std::cerr << "BRANCH 1 " << branchNode.id << std::endl;
                //   recurse on branchNode
                //   Case 1: pick branchNode
                int Fbefore = F.size();
                F.push_back(branchNode);
                g.disable(branchNode);
                execute(g, cf, optimizer,start);
                g.enable(branchNode);
                F.clear(Fbefore);
                
                //std::cerr << "BRANCH 2 " << branchNode.id << std::endl;
                //   Case 2: do not pick branchNode
                weights[branchNode.id] = largeWeight;
                execute(g, cf, optimizer,start);
                weights[branchNode.id] = 1;
                
                //std::cerr << "OUT " << branchNode.id << std::endl;
            }
        }
        
        // restore until the marker
        g.restoreDisabled( lemon::INVALID, true );
        for( int i = F.size()-1; i >= Fsize; i-- ) {
            F.clear(i);
        }
    }
    
private:
    bool strongly_forced( Graph& g, FvsOptimizer& optimizer, const std::list<Graph::Node>& initialSOL ) {
        //std::cerr << "Forced " << std::endl;
        for(std::list<Graph::Node>::const_iterator itn = initialSOL.cbegin(); itn != initialSOL.cend(); itn++) {
            Graph::Node node = static_cast<Graph::Node>(*itn);
            //std::cerr << "Forced " << node.id << std::endl;
            
            // set the weights for bafna
            for( Graph::NodeIt bafnode(g); bafnode != lemon::Invalid(); ++bafnode) {
                bafna_weights[static_cast<Graph::Node>(bafnode)] = weights[bafnode.id];
            }
            bafna_weights[node] = largeWeight;
            
            // call bafna algorithm
            int bafnaLB;
            std::list<Graph::Node> bafnaSOL;
            std::tie(bafnaLB, bafnaSOL) = BafnaEtAl().approximate_fvs(g, optimizer, bafna_weights);
            
            if( F.size() + bafnaLB >= bestSoFarSize ) {
                // need to pick node
                //std::cerr << "Need to " << std::endl;
                F.push_back(node);
                g.disableStore(node);
                return true;
            }
        }
        
        return false;
    }
    
    bool strongly_forced_pair( Graph& g, ComponentFinder& cf, FvsOptimizer& optimizer, const std::list<Graph::Node>& initialSOL, time_t start ) {
        //std::cerr << "Forced pair" << std::endl;
        for(std::list<Graph::Node>::const_iterator itn = initialSOL.cbegin(); itn != initialSOL.cend(); itn++) {
            Graph::Node node = static_cast<Graph::Node>(*itn);
            
            for( Graph::NodeIt ito(g); ito != lemon::Invalid(); ++ito ) {
                Graph::Node other = static_cast<Graph::Node>(ito);
                if( g.isAdjacentDouble(node, other) ) {
                    continue;
                }
                
                //std::cerr << "Forced " << node.id << " " << other.id << std::endl;
                
                // set the weights for bafna
                for( Graph::NodeIt bafnode(g); bafnode != lemon::Invalid(); ++bafnode) {
                    bafna_weights[static_cast<Graph::Node>(bafnode)] = weights[bafnode.id];
                }
                bafna_weights[node] = largeWeight;
                bafna_weights[other] = largeWeight;
                
                // call bafna algorithm
                int bafnaLB;
                std::list<Graph::Node> bafnaSOL;
                std::tie(bafnaLB, bafnaSOL) = BafnaEtAl().approximate_fvs(g, optimizer, bafna_weights);
                
                if( F.size() + bafnaLB >= bestSoFarSize ) {
                    // make node, other a double edge
                    bool isadj = g.isAdjacent(node,other);
                    if( !isadj ) {
                        g.increaseAdjacency(node, other);
                    }
                    g.increaseAdjacency(node, other);
                    
                    //std::cerr << "Recurse " << node.id << " " << other.id << std::endl;
                    execute(g, cf, optimizer, start);
                    
                    // restore adjacency
                    g.decreaseAdjacency(node, other);
                    if( !isadj ) {
                        g.decreaseAdjacency(node, other);
                    }
                    
                    return true;
                }
            }
        }
        
        return false;
    }
    
    // kernelize
    // ret == 0 -> no rule applied
    // ret == 1 -> recurse
    // ret == 2 -> no need to recurse, just return
    short kernelize(Graph& g, ComponentFinder& cf, FvsOptimizer& optimizer, const std::list<Graph::Node>& bafnaSOL, time_t start) {
        bool didSomething = true;
        bool didSomethingNow = false;
        bool globalDidSomething = false;
        
        while( true ) {
            didSomethingNow = remove_nodes_degree_0_1(g, true);
            globalDidSomething = didSomethingNow || globalDidSomething;
            
            while( remove_degree_two_two_neighbors(g, true) ) {
                didSomethingNow = true;
            }
            globalDidSomething = didSomethingNow || globalDidSomething;
            
            didSomethingNow = remove_degree_two_one_neighbor(g, F, true);
            globalDidSomething = didSomethingNow || globalDidSomething;
            if( didSomethingNow ) continue;
            
            didSomethingNow = strongly_forced(g, optimizer, bafnaSOL );
            globalDidSomething = didSomethingNow || globalDidSomething;
            if( didSomethingNow ) break;    // we need a new initial solution now
            
            break;
            
            didSomethingNow = strongly_forced_pair(g, cf, optimizer, bafnaSOL, start );
            globalDidSomething = didSomethingNow || globalDidSomething;
            if( didSomethingNow ) return 2;    // just go back
            
            if( !didSomethingNow ) {
                break;
            }
        }
        
        if( globalDidSomething ) {
            return 1;
        }
        else {
            return 0;
        }
    }
};


#endif /* branch_algo_h */
