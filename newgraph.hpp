/*
*  newgraph.hpp
*  Created by Erik Jan van Leeuwen on 14/7/16.
*
* This file re-uses a large quantity of code that is a part of LEMON, a generic C++ optimization library.
* Copyright (C) 2003-2013
* Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
* (Egervary Research Group on Combinatorial Optimization, EGRES).
*
*/

//

#ifndef newgraph_h
#define newgraph_h

#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/maps.h>
#include <lemon/concept_check.h>
#include <lemon/core.h>
#include <lemon/error.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/bits/map_extender.h>
#include <lemon/bits/default_map.h>

//#include <assert.h>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <stdexcept>
#include <utility>
#include "inplace.hpp"

//template <class node_id_type = int, class edge_id_type = int>
class BaseGraph {
public:
    typedef int node_id_type;
    typedef int edge_id_type;
    
protected:
    template <class id_type, class place_id_type>
    class ElementBase {
        //template <class A, class B>
        friend class BaseGraph;
        
    public:
        id_type id;         // What is the actual id of this element
        place_id_type place; // where in the list are we now?
    public:
        ElementBase() {}
        ElementBase( lemon::Invalid ) : id(Invalid), place(Invalid) {}
        ElementBase( const ElementBase& n ) : id(n.id), place(n.place) {}
        ElementBase( id_type _id, place_id_type _place = Invalid ) : id(_id), place(_place) {}
        
        void invalidate() {
            id = Invalid;
            place = Invalid;
        }
        
        bool is_invalid() const {
            return id == Invalid;
        }
        
        bool is_valid() const {
            return !is_invalid();
        }
        
        bool operator==(const ElementBase& n) const {
            return n.id == id;
        }
        
        bool operator==(const lemon::Invalid& i) const {
            return is_invalid();
        }
        
        bool operator!=(const ElementBase& n) const {
            return n.id != id;
        }
        
        bool operator!=(const lemon::Invalid& i) const {
            return is_valid();
        }
        
        bool operator<(const ElementBase& n) const {
            return id < n.id;
        }
        
        
    };
    
public:
    
    /// The node type of the graph
    /// This class identifies a node of the graph. It also serves
    /// as a base class of the node iterators,
    /// thus they convert to this type.
    class Node : public ElementBase<node_id_type,node_id_type> {//, public std::_List_const_iterator<::BaseGraph::Node *> {
    public:
        Node() : ElementBase<edge_id_type,node_id_type>() {}
        Node(lemon::Invalid i) : ElementBase<edge_id_type,node_id_type>(i) {}
        Node( const Node& e ) : ElementBase<edge_id_type,node_id_type>(e) {}
        Node( node_id_type _id, node_id_type _place = Invalid ) : ElementBase<edge_id_type,node_id_type>(_id,_place) {}
        
        //the operators are implemented by ElementBase
    };
    
    /// The edge type of the graph
    
    /// This class identifies an edge of the graph. It also serves
    /// as a base class of the edge iterators,
    /// thus they will convert to this type.
    class Edge : public ElementBase<edge_id_type,node_id_type> {
    public:
        bool reverseDirection = false;
        
        Edge() : ElementBase<edge_id_type,node_id_type>() {}
        Edge(lemon::Invalid i) : ElementBase<edge_id_type,node_id_type>(i) {}
        Edge( const Edge& e ) : ElementBase<edge_id_type,node_id_type>(e), reverseDirection(e.reverseDirection) {}
        Edge( edge_id_type _id, node_id_type _place = Invalid, bool _reverseDirection = false ) : ElementBase<edge_id_type,node_id_type>(_id,_place), reverseDirection(_reverseDirection) {}
        
        //the operators are implemented by ElementBase
    };
    
    
    
    /// The arc type of the graph
    
    /// This class identifies a directed arc of the graph. It also serves
    /// as a base class of the arc iterators,
    /// thus they will convert to this type.
    class Arc : public ElementBase<edge_id_type,node_id_type> {
    public:
        Arc() : ElementBase<edge_id_type,node_id_type>() {}
        Arc(lemon::Invalid i) : ElementBase<edge_id_type,node_id_type>(i) {}
        Arc( const Arc& e ) : ElementBase<edge_id_type,node_id_type>(e) {}
        Arc( edge_id_type _id, node_id_type _place = Invalid ) : ElementBase<edge_id_type,node_id_type>(_id,_place) {}
        
        operator Edge() const {
            if( ElementBase<edge_id_type,node_id_type>::is_valid() ) {
                if( INPLACE_HAS_BIT(this->id, 0) ) {
                    return Edge(this->id - 1);
                } else {
                    return Edge(id);
                }
            } else {
                return lemon::INVALID;
            }
        }
        
        //the other operators are implemented by ElementBase
    };
    
private:
    
    //inplace_list<node_id_type> nodes;
    inplace_layered_arraylist<node_id_type,node_id_type> active_nodes;
    std::vector<node_id_type> active_nodes_pointers;
    inplace_layered_arraylist<node_id_type,node_id_type> adj_lst;
    inplace_square_matrix<node_id_type> adj_mat;
    inplace_square_matrix_4 adj_cnt;
    inplace_stack<Node,node_id_type> disabledNodeRevertStack;
    std::vector<bool> disabledNodeRevertStackDegreeTwo;
    
    edge_id_type no_edges = 0;
    edge_id_type no_edges_double = 0;
    node_id_type max_no_nodes = 0;
    
    const static edge_id_type Invalid = -1;
    
protected:
    /// Default constructor.
    BaseGraph() : active_nodes(0,0), active_nodes_pointers(0), adj_lst(0,0), adj_mat(0, -1), adj_cnt(0), disabledNodeRevertStack(0), disabledNodeRevertStackDegreeTwo(0) {}
    
    /// BaseGraphs are \e not copy constructible. Use GraphCopy instead.
    BaseGraph(const BaseGraph&) : active_nodes(0,0), active_nodes_pointers(0), adj_lst(0,0), adj_mat(0,-1), adj_cnt(0), disabledNodeRevertStack(0), disabledNodeRevertStackDegreeTwo(0) {}
    /// \brief Assignment of a graph to another one is \e not allowed.
    /// Use GraphCopy instead.
    void operator=(const BaseGraph&) {}
    
    edge_id_type nodes_to_arc( node_id_type u, node_id_type v, bool normalDirection, bool isDoubleEdge ) const {
        if( u > v ) {
            std::swap(u,v);
        }
        return 4 * (u * max_no_nodes + v) + (normalDirection ? 0 : 1) + (isDoubleEdge ? 2 : 0);
    }
    
    edge_id_type nodes_to_edge( node_id_type u, node_id_type v ) const {
        if( u > v ) {
            std::swap(u,v);
        }
        return 4 * (u * max_no_nodes + v);
    }
    
    edge_id_type nodes_to_double_edge( node_id_type u, node_id_type v ) const {
        if( u > v ) {
            std::swap(u,v);
        }
        return 4 * (u * max_no_nodes + v) + 2;
    }
    
    // the first returned node is always the source
    // the second returned node is always the target
    std::tuple<node_id_type,node_id_type> arc_to_nodes( edge_id_type arc_id ) const {
        edge_id_type edge = arc_id;
        if( INPLACE_HAS_BIT(arc_id,0) ) {
            edge -= 1;
        }
        if( INPLACE_HAS_BIT(arc_id, 1) ) {
            edge -= 2;
        }
        node_id_type row = floor( edge / (4*max_no_nodes) );
        node_id_type column = (edge/4) % max_no_nodes;
        if( INPLACE_HAS_BIT(arc_id,0) ) {
            // higher id should come first
            if( column > row ) {
                std::swap(row, column);
            }
        } else {
            // lower id should come first
            if( column < row ) {
                std::swap(row,column);
            }
        }
        return std::make_tuple( row, column );
    }
    
    std::tuple<node_id_type,node_id_type> edge_to_nodes( edge_id_type edge_id ) const {
        if( INPLACE_HAS_BIT(edge_id,1) ) {
            edge_id -= 2;
        }
        node_id_type row = floor( edge_id / (4*max_no_nodes) );
        node_id_type column = (edge_id/4) % max_no_nodes;
        /*if( row > column ) {
            // lower id always comes first
            std::swap( row, column );
        }*/
        return std::make_tuple( row, column );
    }
    
    
public:
    BaseGraph( node_id_type n ) : max_no_nodes(n), active_nodes(1,n), active_nodes_pointers(n), adj_lst( 2*n, n), adj_mat(n,-1), adj_cnt(n), disabledNodeRevertStack(500*n), disabledNodeRevertStackDegreeTwo(500*n) {}
    
    
    /// \brief Undirected graphs should be tagged with \c UndirectedTag.
    ///
    /// Undirected graphs should be tagged with \c UndirectedTag.
    ///
    /// This tag helps the \c enable_if technics to make compile time
    /// specializations for undirected graphs.
    typedef lemon::True UndirectedTag;
    
    
    Node u(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        return Node((u <= v ? u : v));
    }
    
    Node v(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        return Node((u <= v ? v : u));
    }
    
    Node uf(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        Node n((u <= v ? u : v));
        n.place = (u <= v ? active_nodes_pointers.at(u) : active_nodes_pointers.at(v));
        return n;
    }
    
    Node vf(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        Node n((u <= v ? v : u));
        n.place = (u <= v ? active_nodes_pointers.at(v) : active_nodes_pointers.at(u));
        return n;
    }
    
    Node source(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        return Node((u <= v && !e.reverseDirection ? u : v));
    }
    
    Node target(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        return Node((u <= v && !e.reverseDirection ? v : u));
    }
    
    Node sourcef(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        Node n((u <= v && !e.reverseDirection ? u : v));
        n.place = (u <= v && !e.reverseDirection ? active_nodes_pointers.at(u) : active_nodes_pointers.at(v));
        return n;
    }
    
    Node targetf(const Edge& e) const {
        node_id_type u, v;
        std::tie(u,v) = edge_to_nodes(e.id);
        Node n((u <= v && !e.reverseDirection ? v : u));
        n.place = (u <= v && !e.reverseDirection ? active_nodes_pointers.at(v) : active_nodes_pointers.at(u));
        return n;
    }
    
    Node source(const Arc& a) const {
        node_id_type source, target;
        std::tie(source,target) = arc_to_nodes(a.id);
        return Node(source);
    }
    
    Node target(const Arc& a) const {
        node_id_type source, target;
        std::tie(source,target) = arc_to_nodes(a.id);
        return Node(target);
    }
    
    Node sourcef(const Arc& a) const {
        node_id_type source, target;
        std::tie(source,target) = arc_to_nodes(a.id);
        return Node(source, active_nodes_pointers.at(source));
    }
    
    Node targetf(const Arc& a) const {
        node_id_type source, target;
        std::tie(source,target) = arc_to_nodes(a.id);
        return Node(target, active_nodes_pointers.at(target));
    }
    
    int id(const Node& n) const {
        return n.id;
    }
    
    int id(const Edge& e) const {
        return e.id;
    }
    
    int id(const Arc& a) const {
        return a.id;
    }

    Node nodeFromId(int i) const {
        return Node(i);
    }
    
    Edge edgeFromId(int i) const {
        return Edge(i);
    }
    
    Arc arcFromId(int i) const {
        return Arc(i);
    }
    
    int maxNodeId() const {
        return max_no_nodes;
    }
    
    int maxEdgeId() const {
        return nodes_to_double_edge(max_no_nodes, max_no_nodes);
    }
    
    int maxArcId() const {
        return nodes_to_arc(max_no_nodes, max_no_nodes, false, true);
    }
    
    bool direction(const Arc& a) const {
        return !(INPLACE_HAS_BIT(a.id, 0) );
    }
    
    /// \brief Direct the edge.
    ///
    /// Direct the given edge. The returned arc
    /// represents the given edge and its direction comes
    /// from the bool parameter. If it is \c true, then the direction
    /// of the arc is the same as the inherent orientation of the edge.
    Arc direct(const Edge& e, bool b) const {
        if( e.is_valid() ) {
            bool is_double = INPLACE_HAS_BIT(e.id, 1);
            return Arc( nodes_to_arc(u(e).id, v(e).id, b, is_double) );
        } else {
            return lemon::INVALID;
        }
    }
    
    /// \brief Direct the edge.
    ///
    /// Direct the given edge. The returned arc represents the given
    /// edge and its source node is the given node.
    Arc direct(const Edge& e, const Node& n) const {
        if( e.is_valid() && n.is_valid() ) {
            bool is_double = INPLACE_HAS_BIT(e.id, 1);
            if( u(e) == n ) {
                return Arc( nodes_to_arc(u(e).id, v(e).id, true, is_double) );
            } else if( v(e) == n ) {
                return Arc( nodes_to_arc(u(e).id, v(e).id, false, is_double) );
            } else {
                return lemon::INVALID;
            }
        } else {
            return lemon::INVALID;
        }
    }
    
    Arc oppositeArc(const Arc& a) const {
        if( a.is_valid() ) {
            node_id_type u, v;
            std::tie(u,v) = arc_to_nodes(a.id);
            bool is_double = INPLACE_HAS_BIT(a.id, 1);
            return Arc( nodes_to_arc(u, v, !(u < v), is_double) );
        } else {
            return lemon::INVALID;
        }
    }
    
    Node oppositeNode(const Node& n, const Edge& e) const {
        if( n == u(e) ) {
            return v(e);
        } else if( n == v(e) ) {
            return u(e);
        } else {
            return lemon::INVALID;
        }
    }
    
    Node oppositeNode( const Node& n, const Arc& a ) const {
        if( n == source(a) ) {
            return target(a);
        } else if ( n == target(a) ) {
            return source(a);
        } else {
            return lemon::INVALID;
        }
    }
    
    bool hasSelfLoop( const Node &n ) const {
        return isAdjacent(n,n);
    }
    
    bool isAdjacent( const Node &u, const Node &v ) const {
        if( u.is_invalid() || v.is_invalid() ) {
            return false;
        }
        return adj_cnt.get(u.id, v.id) != 0;
    }
    
    bool isAdjacentSingle( const Node &u, const Node &v ) const {
        if( u.is_invalid() || v.is_invalid() ) {
            return false;
        }
        return adj_cnt.get(u.id, v.id) == 1;
    }
    
    bool isAdjacentDouble( const Node &u, const Node &v ) const {
        if( u.is_invalid() || v.is_invalid() ) {
            return false;
        }
        return adj_cnt.get(u.id, v.id) == 2;
    }

    // this does not really work as advertised
//    Edge getEdge( const Node& u, const Node& v ) const {
//        if( !isAdjacent(u, v) ) {
//            Edge e;
//            e.invalidate();
//            return e;
//        } else {
//            return Edge( nodes_to_edge(u.id, v.id), adj_mat.at(u.id,v.id), u.id > v.id );
//        }
//    }
    
    node_id_type degree( const Node& n) const {
        return degreeSingle(n) + 2 * degreeDouble(n);// + (hasSelfLoop(n) ? 1 : 0);
    }
    
    node_id_type degreeSingle(const Node& n) const {
        return adj_lst.size(2*n.id);
    }
    
    node_id_type degreeDouble(const Node& n) const {
        return adj_lst.size(2*n.id+1);
    }
    
    node_id_type degreeUnique(const Node& n) const {
        return degreeSingle(n) + degreeDouble(n);
    }
    
    edge_id_type edgeCount() const {
        return no_edges;
    }
    
    edge_id_type edgeCountDouble() const {
        return no_edges_double;
    }
    
    edge_id_type edgeCountSingle() const {
        return edgeCount() - 2 * edgeCountDouble();
    }
    
    edge_id_type edgeCountUnique() const {
        return edgeCount() - edgeCountDouble();
    }
    
    node_id_type nodeCount() const {
        return active_nodes.size(0);
    }
    
    // returns the node count in the original graph
    node_id_type nodeCountOriginal() const {
        return max_no_nodes;
    }
    
    // enable and disable
    
    std::vector<Node> buildEnableNodes() {
        std::vector<Node> ret;
        ret.reserve( max_no_nodes );
        
        for( node_id_type id = 0; id < max_no_nodes; id++ ) {
            active_nodes.push_back(0, id);
            active_nodes_pointers[id] = id;
            ret.push_back( Node(id) );
        }
        
        return ret;
    }
    
//    Node buildEnable( node_id_type node_id ) {
//        if( node_id != active_nodes.size() ) {
//            throw std::out_of_range("Node id out of range");
//        }
//        Node n = nodeFromId(node_id);
//        enable(n);
//        return n;
//    }
    
    Edge increaseAdjacency( const Edge& e ) {
        return increaseAdjacency(u(e), v(e));
    }
    
    Edge increaseAdjacency( const Node &u, const Node &v ) {
        return increaseAdjacency(u, v, 1, false);
    }
    
protected:
    // doNotListu specifies that we should not modify the list of u
    // the rest always gets modified
    Edge increaseAdjacency( const Node &u, const Node &v, int by, bool doNotUpdateListMatU ) {
        if( u.is_invalid() || v.is_invalid() ) {
            return lemon::INVALID;
        }
        
        int cnt = adj_cnt.get(u.id, v.id);
        if( cnt > 1 ) {
            return Edge( nodes_to_double_edge(u.id, v.id) );
        }
        
        // make sure we set by appropriately
        if( u == v ) {
            // self-loop -> can get max to 1
            by = 2-cnt;
        } else {
            // max 1 or 2, depending on current cnt
            by = std::min((cnt == 1 ? 1 : 2), by);
        }
        
        if( cnt > 1 && by != 0 ) throw 0;
        
        if( by > 0 ) {
            if( cnt == 1 ) {
                // already a single edge -> remove from single edge list
// code is commented to avoid duplication
//                node_id_type oldu, oldv;
//                if( !doNotUpdateListMatU ) {
//                    oldu = adj_mat.at(u.id, v.id);  //oldu is the position in the adjacency list of u where this edge is stored
//                    if( oldu != adj_lst.size(2*u.id) - 1 ) {
//                        // oldu is not the last position in the list, so the order in the list changes
//                        // the guy who is currently last in the adjacency list will move to oldu
//                        node_id_type lastElem = adj_lst.at( 2*u.id, adj_lst.size(2*u.id) - 1 );
//                        if( adj_mat.at(u.id, lastElem) != Invalid ) {
//                            adj_mat.at(u.id, lastElem) = oldu;
//                        } else {
//                            throw 0;
//                        }
//                    }
//                    adj_lst.clear( 2*u.id, oldu );
//                }
//                oldv = adj_mat.at(v.id, u.id);
//                if( oldv != adj_lst.size(2*v.id) - 1 ) {
//                    // oldv is not the last position in the list, so the order in the list changes
//                    // the guy who is currently last in the adjacency list will move to oldv
//                    node_id_type lastElem = adj_lst.at( 2*v.id, adj_lst.size(2*v.id) - 1 );
//                    if( adj_mat.at(v.id, lastElem) != Invalid ) {
//                        adj_mat.at(v.id, lastElem) = oldv;
//                    } else {
//                        throw 0;
//                    }
//                }
//                adj_lst.clear( 2*v.id, oldv );
                decreaseAdjacency(u, v, true, doNotUpdateListMatU);
                cnt = 0;
                by++;
            }
            
            // get the new position in the adjacency list
            node_id_type posu = adj_lst.size( 2*u.id+cnt+by-1 );
            node_id_type posv = adj_lst.size( 2*v.id+cnt+by-1 );
            
            // add to adjacency list
            if( !doNotUpdateListMatU && u != v ) {
                adj_lst.push_back( 2*u.id + cnt+by-1, v.id );
            }
            adj_lst.push_back( 2*v.id + cnt+by-1, u.id );
            
            // add to adjacency matrix
            if( !doNotUpdateListMatU ) {
                adj_mat.at( u.id, v.id ) = posu;
            }
            adj_mat.at( v.id, u.id ) = posv;
            
            no_edges+=by;
            
            // increase the count (this always gets done)
            if( cnt+by == 2 ) {
                adj_cnt.set(u.id,v.id, adj_cnt.value::TWO);
                adj_cnt.set(v.id,u.id, adj_cnt.value::TWO);
                no_edges_double++;
            } else { // cnt + by == 1
                adj_cnt.set(u.id,v.id, adj_cnt.value::ONE);
                adj_cnt.set(v.id,u.id, adj_cnt.value::ONE);
            }
            
            // return the edge
            if( cnt + by == 1 ) {
                return Edge( nodes_to_edge(u.id, v.id) );
            } else { // cnt+by == 2
                return Edge( nodes_to_double_edge(u.id, v.id) );
            }
        } else {
            return lemon::INVALID;
        }
    }
    
    Edge decreaseAdjacency( const Node &u, const Node &v, bool trueDelete, bool doNotUpdateListMatU ) {
        if( u.is_invalid() || v.is_invalid() ) {
            return lemon::INVALID;
        }
        
        int cnt = adj_cnt.get(u.id, v.id);
        if( cnt == 1 ) {
            trueDelete = true;
        }
        
        if( cnt > 0 ) {
            // decrease the count -> this always gets done
            if( trueDelete ) {
                adj_cnt.set(u.id, v.id, adj_cnt.value::ZERO);
                adj_cnt.set(v.id, u.id, adj_cnt.value::ZERO);
            } else {
                adj_cnt.set(u.id, v.id, adj_cnt.value::ONE);
                adj_cnt.set(v.id, u.id, adj_cnt.value::ONE);
            }
            
            // remove from the current list
            node_id_type oldu, oldv;
            if( !doNotUpdateListMatU ) {
                oldu = adj_mat.at(u.id, v.id);  //oldu is the position in the adjacency list of u where this edge is stored
                if( oldu != adj_lst.size(2*u.id+(cnt-1)) - 1 ) {
                    // oldu is not the last position in the list, so the order in the list changes
                    // the guy who is currently last in the adjacency list will move to oldu
                    node_id_type lastElem = adj_lst.at( 2*u.id+(cnt-1), adj_lst.size(2*u.id+(cnt-1)) - 1 );
                    if( adj_mat.at(u.id, lastElem) != Invalid ) {
                        adj_mat.at(u.id, lastElem) = oldu;
                    }
                }
                adj_lst.clear( 2*u.id+(cnt-1), oldu );
            }
            oldv = adj_mat.at(v.id, u.id);
            if( oldv != adj_lst.size(2*v.id+(cnt-1)) - 1 ) {
                // oldv is not the last position in the list, so the order in the list changes
                // the guy who is currently last in the adjacency list will move to oldv
                node_id_type lastElem = adj_lst.at( 2*v.id+(cnt-1), adj_lst.size(2*v.id+(cnt-1)) - 1 );
                if( adj_mat.at(v.id, lastElem) != Invalid ) {
                    adj_mat.at(v.id, lastElem) = oldv;
                }
            }
            adj_lst.clear( 2*v.id+(cnt-1), oldv );

            
            if( cnt == 2 && !trueDelete ) {
                node_id_type posu, posv;
                
                // insert into single adjacency list
                if( !doNotUpdateListMatU ) {
                    posu = adj_lst.size( 2*u.id );
                }
                posv = adj_lst.size( 2*v.id );
            
                // add to single adjacency list
                if( !doNotUpdateListMatU ) {
                    adj_lst.push_back( 2*u.id, v.id );
                }
                adj_lst.push_back( 2*v.id, u.id );
            
                // add to single adjacency matrix
                if( !doNotUpdateListMatU ) {
                    adj_mat.at( u.id, v.id ) = posu;
                }
                adj_mat.at( v.id, u.id ) = posv;
            } else {
                // we are really deleting, clear the adjacency matrix
                if( !doNotUpdateListMatU ) {
                    adj_mat.at( u.id, v.id ) = Invalid;
                }
                adj_mat.at( v.id, u.id ) = Invalid;
            }
            
            if( trueDelete ) {
                no_edges -= cnt;
            } else {
                no_edges--;
            }
            
            // return the edge
            if( cnt == 2 ) {
                no_edges_double--;
                return Edge( nodes_to_edge(u.id, v.id) );
            } else { // cnt == 1
                return lemon::INVALID;
            }
        } else {
            return lemon::INVALID;
        }
    }
    
public:
    Edge decreaseAdjacency( const Node &u, const Node &v ) {
        return decreaseAdjacency(u, v, false, false);
    }
    
    Edge decreaseAdjacency( const Edge& e ) {
        return decreaseAdjacency(u(e), v(e));
    }
    
protected:
    bool disableHelper( const Node& n ) {
        if( n.is_invalid() ) {
            return false;
        }
        int& activeLocation = active_nodes_pointers.at(n.id);
        if( activeLocation == Invalid ) {
            return false;
        }
        
        //idea: keep the adjacency list of n intact, and just remove n from its neighbors
        for( int i = 0; i < adj_lst.size(2*n.id); i++ ) {
            Node on( adj_lst.at(2*n.id,i) );
            decreaseAdjacency(n, on, true, true);
        }
        for( int i = 0; i < adj_lst.size(2*n.id+1); i++ ) {
            Node on( adj_lst.at(2*n.id+1,i) );
            decreaseAdjacency(n, on, true, true);
        }
        
        // recall that if we remove the node from the active pointers, then the last item in active nodes might get moved
        // hence, we might need to update active_nodes_pointers for the guy at the end
        if( activeLocation != active_nodes.size(0) - 1 ) {
            active_nodes_pointers.at( active_nodes.at(0,active_nodes.size(0)-1) ) = activeLocation;
        }
        active_nodes.clear(0, activeLocation);
        //active_nodes_pointers.at(n.id) = Invalid;
        activeLocation = Invalid;
        
        return true;
    }

public:
    void disable( const Node& n ) {
        disableHelper(n);
    }
    
    void disableStore( const Node& n ) {
        if( disableHelper(n) ) {
            disabledNodeRevertStack.push(n);
            disabledNodeRevertStackDegreeTwo.at(disabledNodeRevertStack.size()-1) = false;
        }
    }
    
    void disableMark() {
        Node n(-2,0);
        disabledNodeRevertStack.push(n);
        disabledNodeRevertStackDegreeTwo.at(disabledNodeRevertStack.size()-1) = false;
    }
    
    // restores disabled nodes
    void restoreDisabled() {
        if( !disabledNodeRevertStack.empty() ) {
            Node n = disabledNodeRevertStack.at(0);
            restoreDisabled(n);
        }
    }
    
    // restores disabled nodes until (and including) node
    void restoreDisabled( const Node& node, bool untilMark = false ) {
        while( !disabledNodeRevertStack.empty() ) {
            Node n = disabledNodeRevertStack.top();
            disabledNodeRevertStack.pop();
            if( n.id == -2 ) {
                if( untilMark ) {
                    break;
                } else {
                    continue;
                }
            }
            enable(n);
            bool shouldRestoreAdjacency = disabledNodeRevertStackDegreeTwo.at(disabledNodeRevertStack.size());
            if( shouldRestoreAdjacency ) {
                Node u( adj_lst.at(2*n.id,0) );
                Node v( adj_lst.at(2*n.id,1) );
                decreaseAdjacency(u, v);
            }
            if( !untilMark && n == node ) {
                break;
            }
        }
    }
    
//    // removes disabled nodes from the restore stack until (but NOT including) node
//    void restoreDisabledRollback( const Node& node ) {
//        while( !disabledNodeRevertStack.empty() ) {
//            Node n = disabledNodeRevertStack.top();
//            if( n == node ) {
//                break;
//            }
//            disabledNodeRevertStack.pop();
//            disabledNodeRevertStackDegreeTwo.pop();
//        }
//    }
    
    void enable( const Node& n ) {
        if( n.is_invalid() ) {
            return;
        }
        
        //idea: the adjacency list of n is still intact
        for( int i = 0; i < adj_lst.size(2*n.id); i++ ) {
            Node on( adj_lst.at(2*n.id,i) );
            if( isEnabled(on) ) {
                increaseAdjacency(n, on, 1, true);
            }
        }
        for( int i = 0; i < adj_lst.size(2*n.id+1); i++ ) {
            Node on( adj_lst.at(2*n.id+1,i) );
            if( isEnabled(on) ) {
                increaseAdjacency(n, on, 2, true);
            }
        }
        
        active_nodes_pointers.at(n.id) = active_nodes.size(0);
        active_nodes.push_back(0, n.id);
    }
    
    bool isEnabled( const Node& n ) const {
        if( n.is_invalid() ) {
            return false;
        } else {
            return active_nodes_pointers.at(n.id) != Invalid;
        }
    }
    
    // returns true iff e is a double edge
    bool disable( const Edge& e, bool deleteBothIfDouble = false ) {
        if( e.is_invalid() ) {
            return false;
        }
        Node n = u(e);
        Node nn = v(e);
        bool is_double = isAdjacentDouble(n, nn);
        decreaseAdjacency(n, nn, is_double && deleteBothIfDouble, false);
        return is_double;
    }
    
    void enable( const Edge& e, bool restoreBothIfDouble = true ) {
        if( e.is_invalid() ) {
            return;
        }
        Node n = u(e);
        Node nn = v(e);
        bool is_double = isAdjacentDouble(n, nn);
        int by = (is_double && restoreBothIfDouble ? 2 : 1);
        increaseAdjacency(n, nn, by, false);
    }
    
    /*// returns the number of parallel edges that were made
    void contract( const Edge& e, bool allowParallelEdges = true ) {
        if( e.is_invalid() ) {
            return;
        }
        contract(u(e),v(e),allowParallelEdges);
    }
    
    void contract( const Node& u, const Node& v, bool allowParallelEdges = true ) {
        if( u.is_invalid() || v.is_invalid() ) {
            return;
        }
        
        // make sure we take the node with the smallest degree
        node_id_type nid = u.id;
        if( degreeUnique(u) > degreeUnique(v) ) {
            nid = v.id;
        }
        Node n(nid);
        Node other(nid == u.id ? v.id : u.id);
        
        // now loop over the neighbors of n and make them adjacent to other
        for( int i = 0; i < adj_lst.size(2*n.id); i++ ) {
            Node on( adj_lst.at(2*n.id,i) );
            if( on != n ) {
                if( allowParallelEdges || !isAdjacent(other, on) ) {
                    increaseAdjacency(other, on, 1, false);
                }
            }
        }
        if( allowParallelEdges ) {
            for( int i = 0; i < adj_lst.size(2*n.id+1); i++ ) {
                Node on( adj_lst.at(2*n.id+1,i) );
                if( on != n ) {
                    increaseAdjacency(other, on, 2, false);
                }
            }
        }
    }*/
    
    std::tuple<Node,Node,bool> contractTwoStore( const Node &n ) {
        Node u,v;
        bool shouldRestoreAdjacency;
        std::tie(u,v,shouldRestoreAdjacency) = contractTwo(n);
        if( u.is_valid() && v.is_valid() && n.is_valid() ) {
            disabledNodeRevertStack.push(n);
            disabledNodeRevertStackDegreeTwo.at(disabledNodeRevertStack.size()-1) = shouldRestoreAdjacency;
        }
        return std::make_tuple(u,v,shouldRestoreAdjacency);
    }
    
    // contracts a vertex of degree 2
    // returns the neighbors of the vertex and a bool to indicate whether the neighbors were already double-adjacent before
    std::tuple<Node,Node,bool> contractTwo( const Node &n ) {
        if( n.is_invalid() || degree(n) != 2 || degreeUnique(n) != 2 ) {
            return std::make_tuple(lemon::INVALID,lemon::INVALID,false);
        }
        
        Node u( adj_lst.at(2*n.id,0) );
        Node v( adj_lst.at(2*n.id,1) );
        bool shouldRestoreAdjacency = !isAdjacentDouble(u, v);
        
        disable(n);
        increaseAdjacency(u, v);
        return std::make_tuple(u,v,shouldRestoreAdjacency);
    }
    
    // undo a contraction of a node of degree 2
    // the boolean indicates whether the neighbors of n were double-adjacent prior to contraction
    void uncontractTwo( const Node &n, bool wereNeighborsAdjacentDouble ) {
        if( n.is_invalid() || degree(n) != 2 || degreeUnique(n) != 2 ) {
            return;
        }
        
        enable(n);
        if( !wereNeighborsAdjacentDouble ) {
            Node u( adj_lst.at(2*n.id,0) );
            Node v( adj_lst.at(2*n.id,1) );
            decreaseAdjacency(u, v);
        }
    }
    
    
    
    // iterators
    
    void first(Node& n) const {
        if( active_nodes.empty(0) ) {
            n.invalidate();
        } else {
            n.id = active_nodes.at(0,0);
            n.place = 0;
        }
    }
    void next(Node& n) const {
        if( n.is_valid() && n.place != Invalid ) {
            if( n.place < active_nodes.size(0)-1 ) {
                n.place++;
                n.id = active_nodes.at(0,n.place);
            } else {
                n.invalidate();
            }
        } else {
            n.invalidate();
        }
    }
    Node firstNode() const {
        Node n;
        first(n);
        return n;
    }
   
protected:
    // tries to find the next edge
    // can skip to the next node if necessary
    // if eandnrelated, then e cannot be invalid
    void nextHelperEdge(Node& n, Edge& e, bool eandnrelated = false, bool allowIterate = true, bool careForDirection = false, bool onlySingle = false, bool onlyDouble = false, bool uniqueNeighbors = false ) const {
        bool allowIterateOnce = !allowIterate;
        while( n.is_valid() && (allowIterate || allowIterateOnce) ) {
            allowIterateOnce = false;
            if( (onlySingle && degreeSingle(n) > 0) || (onlyDouble && degreeDouble(n) > 0) || (!onlySingle && !onlyDouble && degree(n) > 0) ) {
                // first check whether e is the 'first' part of a double edge
                if( eandnrelated && !onlySingle && !onlyDouble ) {
                    Node nn = oppositeNode(n, e);
                    if( !INPLACE_HAS_BIT(e.id, 1) && (e.place >= adj_lst.size(2*n.id) || adj_lst.at(2*n.id,e.place) != nn.id) ) {
                        // marked as a single edge, but e.place does not fit in the single-adjacency list or what is in the list does not match the nn.id
                        // move to 'second' part of double edge
                        e.id += 2;
                        return;
                    }
                }
                // now we just need to proceed in the same list as where e is really from
                bool f_done = false;
                if( !onlyDouble && ( !eandnrelated || !INPLACE_HAS_BIT(e.id, 1) ) ) {
                    // if eandnrelated is false, then we start from the beginning with single-adjacent nodes
                    // if there is a relation, then e is a single edge and we continue with single-adjacent nodes
                    node_id_type start = (eandnrelated ? e.place+1 : 0);
                    for( node_id_type i = start; i < adj_lst.size(2*n.id); i++ ) {
                        node_id_type node_id = adj_lst.at(2*n.id, i);
                        if( n.id <= node_id || !careForDirection ) {
                            // only report the edge if n is the 'source'
                            e.id = nodes_to_edge(n.id, node_id);
                            e.place = i;
                            e.reverseDirection = (n.id > node_id);
                            f_done = true;
                            break;
                        }
                    }
                }
                if( !onlySingle && !f_done ) {
                   //( (!eandnrelated && !f_done) || (eandnrelated && INPLACE_HAS_BIT(e.id, 1)) || (eandnrelated && !INPLACE_HAS_BIT(e.id, 1) && !f_done) ) ) {
                    // if eandnrelated is false and we were not done before, then we start from the beginning with double-adjacent nodes
                    // if eandnrelated is true and e is a double edge, then we find the next double edge
                    node_id_type start = (eandnrelated && INPLACE_HAS_BIT(e.id, 1) ? e.place+1 : 0);
                    for( node_id_type i = start; i < adj_lst.size(2*n.id+1); i++ ) {
                        node_id_type node_id = adj_lst.at(2*n.id+1, i);
                        if( n.id <= node_id || !careForDirection ) {
                            // only report the edge if n is the 'source'
                            // take care to first report the 'first' part of the double edge
                            if( !uniqueNeighbors && !onlyDouble ) {
                                e.id = nodes_to_edge(n.id, node_id);
                            } else {
                                e.id = nodes_to_double_edge(n.id, node_id);
                            }
                            e.place = i;
                            e.reverseDirection = (n.id > node_id);
                            f_done = true;
                            break;
                        }
                    }
                }
                if( f_done ) {
                    break;
                }
            }
            if( allowIterate ) {
                next(n);
                eandnrelated = false;
            } else {
                e.invalidate();
                return;
            }
        }
        if( n.is_invalid() ) {
            e.invalidate();
        }
    }
    
public:
    // just a reminder
    // void nextHelperEdge(Node, Edge, eandnrelated = false, allowIterate = true, careForDirection = false, onlySingle = false, onlyDouble = false )
    
    void first(Edge& e) const {
        Node n;
        first(n);
        if( n.is_valid() ) {
            nextHelperEdge(n,e, false, true, true);
        } else {
            e.invalidate();
        }
    }
    void next(Edge& e) const {
        if( e.is_valid() ) {
            Node n = uf(e);
            if( n.is_valid() ) {
                nextHelperEdge(n, e, true, true, true);
            } else {
                e.invalidate();
            }
        }
    }
    
    // TODO: might not need reverseDirection, due to the boolean d
    void firstInc( Edge& e, bool& d, const Node& n ) const {
        if( n.is_valid() && degree(n) > 0 ) {
            Node nn(n);
            nextHelperEdge(nn, e, false, false);
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    void nextInc(Edge& e, bool& d) const {
        if( e.is_valid() ) {
            if( e.reverseDirection ) {
                Node n = vf(e);
                nextHelperEdge(n, e, true, false);
            } else {
                Node n = uf(e);
                nextHelperEdge(n, e, true, false);
            }
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    
    // just a reminder
    // void nextHelperEdge(Node, Edge, eandnrelated = false, allowIterate = true, careForDirection = false, onlySingle = false, onlyDouble = false )
    
    void firstDouble(Edge& e) const {
        Node n;
        first(n);
        if( n.is_valid() ) {
            nextHelperEdge(n,e,false,true,true,false,true);
        } else {
            e.invalidate();
        }
    }
    void nextDouble(Edge& e) const {
        if( e.is_valid() ) {
            Node n = uf(e);
            if( n.is_valid() ) {
                nextHelperEdge(n, e, true,true,true,false,true);
            } else {
                e.invalidate();
            }
        }
    }
    
    void firstDoubleInc( Edge& e, bool& d, const Node& n ) const {
        if( n.is_valid() && degreeDouble(n) > 0 ) {
            Node nn(n);
            nextHelperEdge(nn, e, false, false, false, false, true);
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    void nextDoubleInc(Edge& e, bool& d) const {
        if( e.is_valid() ) {
            if( e.reverseDirection ) {
                Node n = vf(e);
                nextHelperEdge(n, e, true, false, false, false, true);
            } else {
                Node n = uf(e);
                nextHelperEdge(n, e, true, false, false, false, true);
            }
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    
    // just a reminder
    // void nextHelperEdge(Node, Edge, eandnrelated = false, allowIterate = true, careForDirection = false, onlySingle = false, onlyDouble = false )
    
    void firstSingle(Edge& e) const {
        Node n;
        first(n);
        if( n.is_valid() ) {
            nextHelperEdge(n,e,false,true,true,true,false);
        } else {
            e.invalidate();
        }
    }
    void nextSingle(Edge& e) const {
        if( e.is_valid() ) {
            Node n = uf(e);
            if( n.is_valid() ) {
                nextHelperEdge(n, e, true,true,true,true,false);
            } else {
                e.invalidate();
            }
        }
    }
    
    void firstSingleInc( Edge& e, bool& d, const Node& n ) const {
        if( n.is_valid() && degreeSingle(n) > 0 ) {
            Node nn(n);
            nextHelperEdge(nn, e, false, false, false, true, false);
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    void nextSingleInc(Edge& e, bool& d) const {
        if( e.is_valid() ) {
            if( e.reverseDirection ) {
                Node n = vf(e);
                nextHelperEdge(n, e, true, false, false, true, false);
            } else {
                Node n = uf(e);
                nextHelperEdge(n, e, true, false, false, true, false);
            }
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    
    // just a reminder
    // void nextHelperEdge(Node, Edge, eandnrelated = false, allowIterate = true, careForDirection = false, onlySingle = false, onlyDouble = false )
    
    void firstUnique(Edge& e) const {
        Node n;
        first(n);
        if( n.is_valid() ) {
            nextHelperEdge(n,e,false,true,true,false,false,true);
        } else {
            e.invalidate();
        }
    }
    void nextUnique(Edge& e) const {
        if( e.is_valid() ) {
            Node n = uf(e);
            if( n.is_valid() ) {
                nextHelperEdge(n, e, true,true,true,false,false,true);
            } else {
                e.invalidate();
            }
        }
    }
    
    void firstUniqueInc( Edge& e, bool& d, const Node& n ) const {
        if( n.is_valid() && degreeUnique(n) > 0 ) {
            Node nn(n);
            nextHelperEdge(nn, e, false, false, false, false, false, true);
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    void nextUniqueInc(Edge& e, bool& d) const {
        if( e.is_valid() ) {
            if( e.reverseDirection ) {
                Node n = vf(e);
                nextHelperEdge(n, e, true, false, false, false, false, true);
            } else {
                Node n = uf(e);
                nextHelperEdge(n, e, true, false, false, false, false, true);
            }
            d = !e.reverseDirection;
        } else {
            e.invalidate();
        }
    }
    
//    void firstNeighbor( Node &current, const Node& n ) const {
//        if( n.is_valid() && degreeUnique(n) > 0 ) {
//            Edge e;
//            bool d;
//            firstUniqueInc(e, d, n);
//            if( e.is_valid() ) {
//                current = oppositeNode(n, e);
//            } else {
//                current.invalidate();
//            }
//        } else {
//            current.invalidate();
//        }
//    }
//    void nextNeighbor( Node& current, const Node& n ) const {
//        if( current.is_valid() && n.is_valid() ) {
//            if( n.id == 30 || current.id == 30 ) {
//                std::cerr << n.id << " " << current.id << std::endl;
//            }
//            Edge e = getEdge(n, current);
//            e.place = adj_mat.at( n.id, current.id );
//            bool d;
//            nextUniqueInc(e,d);
//            if( e.is_valid() ) {
//                current = oppositeNode(n, e);
//            } else {
//                current.invalidate();
//            }
//        } else {
//            current.invalidate();
//        }
//    }


protected:
    // tries to find the next arc
    // can skip to the next node if necessary
    // if eandnrelated, then e cannot be invalid
    void nextHelperArc(Node& n, Arc& a, bool aandnrelated = false, bool allowIterate = true, bool forceDirection = false, bool forceOutward = false ) const {
        bool allowIterateOnce = !allowIterate;
        while( n.is_valid() && (allowIterate || allowIterateOnce) ) {
            allowIterateOnce = false;
            if( degree(n) > 0 ) {
                bool f_done = false;
                if( !aandnrelated || !INPLACE_HAS_BIT(a.id, 1) ) {
                    // if aandnrelated is false, then we start from the beginning with single-adjacent nodes
                    // if there is a relation, then e is a single edge and we continue with single-adjacent nodes
                    node_id_type start = (aandnrelated ? a.place+1 : 0);
                    for( node_id_type i = start; i < adj_lst.size(2*n.id); i++ ) {
                        node_id_type node_id = adj_lst.at(2*n.id, i);
                        if( n.id <= node_id || forceDirection ) {
                            // only report the arc if it is a forward arc, or we should force the direction
                            bool direction = (n.id <= node_id);
                            if( forceDirection && !forceOutward ) {
                                direction = (n.id > node_id);
                            } else {
                                direction = (n.id <= node_id);
                            }
                            a.id = nodes_to_arc(n.id, node_id, direction, false);
                            a.place = i;
                            f_done = true;
                            break;
                        }
                    }
                }
                //if( (!aandnrelated && !f_done) || (aandnrelated && INPLACE_HAS_BIT(a.id, 1)) || (aandnrelated && !INPLACE_HAS_BIT(a.id, 1) && !f_done) ) {
                if( !f_done ) {
                    // if aandnrelated is false and we were not done before, then we start from the beginning with double-adjacent nodes
                    // if aandnrelated is true and e is a double edge, then we find the next double edge
                    // if aandnrelated is true and e is not a double edge, but we are not done, then try to find the first double edge
                    node_id_type start = (aandnrelated && INPLACE_HAS_BIT(a.id, 1) ? a.place+1 : 0);
                    for( node_id_type i = start; i < adj_lst.size(2*n.id+1); i++ ) {
                        node_id_type node_id = adj_lst.at(2*n.id+1, i);
                        if( n.id <= node_id || forceDirection ) {
                            // only report the arc if it is a forward arc, or we should force the direction
                            bool direction = (n.id <= node_id);
                            if( forceDirection && !forceOutward ) {
                                direction = (n.id > node_id);
                            } else {
                                direction = (n.id <= node_id);
                            }
                            a.id = nodes_to_arc(n.id, node_id, direction, true);
                            a.place = i;
                            f_done = true;
                            break;
                        }
                    }
                }
                if( f_done ) {
                    break;
                }
            }
            if( allowIterate ) {
                next(n);
                aandnrelated = false;
            } else {
                a.invalidate();
                return;
            }
        }
        if( n.is_invalid() ) {
            a.invalidate();
        }
    }

public:
    
    // just a reminder
    // void nextHelperArc(Node Arc, aandnrelated = false, allowIterate = true, forceDirection = false, forceOutward = false )
    
    void first(Arc& a) const {
        Node n;
        first(n);
        if( n.is_valid() ) {
            // we are only interested in the naturally forward arcs
            nextHelperArc(n, a, false);
        } else {
            a.invalidate();
        }
    }
    void next(Arc& a) const {
        if( a.is_valid() ) {
            Node n = sourcef(a);
            if( n.is_valid() ) {
                // we are only interested in the naturally forward arcs
                nextHelperArc(n, a, true);
            } else {
                a.invalidate();
            }
        }
    }
  
    void firstOut(Arc& a, const Node &n) const {
        if( n.is_valid() ) {
            Node nn(n);
            nextHelperArc(nn, a, false, false, true, true);
        } else {
            a.invalidate();
        }
    }
    void nextOut(Arc& a) const {
        if( a.is_valid() ) {
            Node n = sourcef(a);
            // try to go to the next arc
            if( n.is_valid() ) {
                nextHelperArc(n, a, true, false, true, true);
            } else {
                a.invalidate();
            }
        }
    }
    
    void firstIn(Arc& a, const Node& n) const {
        if( n.is_valid() ) {
            Node nn(n);
            nextHelperArc(nn, a, false, false, true, false);
        } else {
            a.invalidate();
        }
    }
    void nextIn(Arc& a) const {
        if( a.is_valid() ) {
            Node n = sourcef(a);
            // try to go to the next arc
            if( n.is_valid() ) {
                nextHelperArc(n, a, true, false, true, false);
            } else {
                a.invalidate();
            }
        }
    }
    
    // The second parameter is dummy.
    Node fromId(int i, Node) const {
        return nodeFromId(i);
    }
    // The second parameter is dummy.
    Edge fromId(int i, Edge) const {
        return edgeFromId(i);
    }
    // The second parameter is dummy.
    Arc fromId(int i, Arc) const {
        return arcFromId(i);
    }
    
    // Dummy parameter.
    int maxId(Node) const {
        return maxNodeId();
    }
    // Dummy parameter.
    int maxId(Edge) const {
        return maxEdgeId();
    }
    // Dummy parameter.
    int maxId(Arc) const {
        return maxArcId();
    }
    
    Node addNode() const {
        return lemon::INVALID;
    }
    
    Edge addEdge(const Node& from, const Node& to) const {
        return lemon::INVALID;
    }
};

//template <typename Base>
class BaseGraphExtender : public BaseGraph {
    typedef BaseGraph Parent;
    
public:
    typedef BaseGraphExtender Graph;
    
    typedef lemon::True UndirectedTag;
    
    typedef typename Parent::Node Node;
    typedef typename Parent::Arc Arc;
    typedef typename Parent::Edge Edge;
    
    class SingleEdgeIt : public Parent::Edge {
        const Graph* _graph;
    public:
        
        SingleEdgeIt() { }
        
        SingleEdgeIt(lemon::Invalid i) : Edge(i) { }
        
        explicit SingleEdgeIt(const Graph& graph) : _graph(&graph) {
            _graph->firstSingle(static_cast<Edge&>(*this));
        }
        
        SingleEdgeIt(const Graph& graph, const Edge& edge) :
        Edge(edge), _graph(&graph) { }
        
        SingleEdgeIt& operator++() {
            _graph->nextSingle(*this);
            return *this;
        }
        
    };
    
    class IncSingleEdgeIt : public Parent::Edge {
        friend class BaseGraphExtender;
        const Graph* _graph;
        bool _direction;
    public:
        
        IncSingleEdgeIt() { }
        
        IncSingleEdgeIt(lemon::Invalid i) : Edge(i), _direction(false) { }
        
        IncSingleEdgeIt(const Graph& graph, const Node &node) : _graph(&graph) {
            _graph->firstSingleInc(*this, _direction, node);
        }
        
        IncSingleEdgeIt(const Graph& graph, const Edge &edge, const Node &node)
        : _graph(&graph), Edge(edge) {
            _direction = (_graph->source(edge) == node);
        }
        
        IncSingleEdgeIt& operator++() {
            _graph->nextSingleInc(*this, _direction);
            return *this;
        }
    };

    class DoubleEdgeIt : public Parent::Edge {
        const Graph* _graph;
    public:
        
        DoubleEdgeIt() { }
        
        DoubleEdgeIt(lemon::Invalid i) : Edge(i) { }
        
        explicit DoubleEdgeIt(const Graph& graph) : _graph(&graph) {
            _graph->firstDouble(static_cast<Edge&>(*this));
        }
        
        DoubleEdgeIt(const Graph& graph, const Edge& edge) :
        Edge(edge), _graph(&graph) { }
        
        DoubleEdgeIt& operator++() {
            _graph->nextDouble(*this);
            return *this;
        }
        
    };
    
    class IncDoubleEdgeIt : public Parent::Edge {
        friend class BaseGraphExtender;
        const Graph* _graph;
        bool _direction;
    public:
        
        IncDoubleEdgeIt() { }
        
        IncDoubleEdgeIt(lemon::Invalid i) : Edge(i), _direction(false) { }
        
        IncDoubleEdgeIt(const Graph& graph, const Node &node) : _graph(&graph) {
            _graph->firstDoubleInc(*this, _direction, node);
        }
        
        IncDoubleEdgeIt(const Graph& graph, const Edge &edge, const Node &node)
        : _graph(&graph), Edge(edge) {
            _direction = (_graph->source(edge) == node);
        }
        
        IncDoubleEdgeIt& operator++() {
            _graph->nextDoubleInc(*this, _direction);
            return *this;
        }
    };
    
    class UniqueEdgeIt : public Parent::Edge {
        const Graph* _graph;
    public:
        
        UniqueEdgeIt() { }
        
        UniqueEdgeIt(lemon::Invalid i) : Edge(i) { }
        
        explicit UniqueEdgeIt(const Graph& graph) : _graph(&graph) {
            _graph->firstUnique(static_cast<Edge&>(*this));
        }
        
        UniqueEdgeIt(const Graph& graph, const Edge& edge) :
        Edge(edge), _graph(&graph) { }
        
        UniqueEdgeIt& operator++() {
            _graph->nextUnique(*this);
            return *this;
        }
        
    };
    
    class IncUniqueEdgeIt : public Parent::Edge {
        friend class BaseGraphExtender;
        const Graph* _graph;
        bool _direction;
    public:
        
        IncUniqueEdgeIt() { }
        
        IncUniqueEdgeIt(lemon::Invalid i) : Edge(i), _direction(false) { }
        
        IncUniqueEdgeIt(const Graph& graph, const Node &node) : _graph(&graph) {
            _graph->firstUniqueInc(*this, _direction, node);
        }
        
        IncUniqueEdgeIt(const Graph& graph, const Edge &edge, const Node &node)
        : _graph(&graph), Edge(edge) {
            _direction = (_graph->source(edge) == node);
        }
        
        IncUniqueEdgeIt& operator++() {
            _graph->nextUniqueInc(*this, _direction);
            return *this;
        }
    };
    
//    class NeighborIt : public Node {
//        const Graph* _graph;
//        const Node itn;
//    public:
//        
//        NeighborIt() {}
//        
//        NeighborIt(lemon::Invalid i) : Node(i) { }
//        
//        explicit NeighborIt(const Graph& graph, const Node& n) : _graph(&graph), itn(n) {
//            _graph->firstNeighbor(static_cast<Node&>(*this),itn);
//        }
//        
//        NeighborIt(const Graph& graph, const Node& node, const Node& n)
//        : Node(node), _graph(&graph), itn(n) {}
//        
//        NeighborIt& operator++() {
//            _graph->nextNeighbor(*this,itn);
//            return *this;
//        }
//    };
    
    class NeighborIt : public Node {
        const Graph* _graph;
        IncUniqueEdgeIt it;
        const Node itn;
    public:
        
        NeighborIt() {}
        
        NeighborIt(lemon::Invalid i) : Node(i) { }
        
        explicit NeighborIt(const Graph& graph, const Node& n) : _graph(&graph), it(graph, n), itn(n) {
            this->id = _graph->oppositeNode(itn, it).id;
        }
        
        //NeighborIt(const Graph& graph, const Node& node, const Node& n)
        //: Node(node), _graph(&graph), itn(n) {}
        
        NeighborIt& operator++() {
            ++it;
            if( static_cast<Edge>(it).is_invalid() ) {
                this->invalidate();
            } else {
                this->id = _graph->oppositeNode(itn, it).id;
            }
            return *this;
        }
    };
    
    typedef lemon::AlterationNotifier<BaseGraphExtender, Node> NodeNotifier;
    typedef lemon::AlterationNotifier<BaseGraphExtender, Arc> ArcNotifier;
    typedef lemon::AlterationNotifier<BaseGraphExtender, Edge> EdgeNotifier;
    
    
protected:
    
    mutable NodeNotifier node_notifier;
    mutable ArcNotifier arc_notifier;
    mutable EdgeNotifier edge_notifier;
    
public:
    
    NodeNotifier& notifier(Node) const {
        return node_notifier;
    }
    
    ArcNotifier& notifier(Arc) const {
        return arc_notifier;
    }
    
    EdgeNotifier& notifier(Edge) const {
        return edge_notifier;
    }
    
    
    
    class NodeIt : public Node {
        const Graph* _graph;
    public:
        
        NodeIt() {}
        
        NodeIt(lemon::Invalid i) : Node(i) { }
        
        explicit NodeIt(const Graph& graph) : _graph(&graph) {
            _graph->first(static_cast<Node&>(*this));
        }
        
        NodeIt(const Graph& graph, const Node& node)
        : Node(node), _graph(&graph) {}
        
        NodeIt& operator++() {
            _graph->next(*this);
            return *this;
        }
    };
    
    
    class ArcIt : public Arc {
        const Graph* _graph;
    public:
        
        ArcIt() { }
        
        ArcIt(lemon::Invalid i) : Arc(i) { }
        
        explicit ArcIt(const Graph& graph) : _graph(&graph) {
            _graph->first(static_cast<Arc&>(*this));
        }
        
        ArcIt(const Graph& graph, const Arc& arc) :
        Arc(arc), _graph(&graph) { }
        
        ArcIt& operator++() {
            _graph->next(*this);
            return *this;
        }
        
    };
    
    
    class OutArcIt : public Arc {
        const Graph* _graph;
    public:
        
        OutArcIt() { }
        
        OutArcIt(lemon::Invalid i) : Arc(i) { }
        
        OutArcIt(const Graph& graph, const Node& node)
        : _graph(&graph) {
            _graph->firstOut(*this, node);
        }
        
        OutArcIt(const Graph& graph, const Arc& arc)
        : Arc(arc), _graph(&graph) {}
        
        OutArcIt& operator++() {
            _graph->nextOut(*this);
            return *this;
        }
        
    };
    
    
    class InArcIt : public Arc {
        const Graph* _graph;
    public:
        
        InArcIt() { }
        
        InArcIt(lemon::Invalid i) : Arc(i) { }
        
        InArcIt(const Graph& graph, const Node& node)
        : _graph(&graph) {
            _graph->firstIn(*this, node);
        }
        
        InArcIt(const Graph& graph, const Arc& arc) :
        Arc(arc), _graph(&graph) {}
        
        InArcIt& operator++() {
            _graph->nextIn(*this);
            return *this;
        }
        
    };
    
    
    class EdgeIt : public Parent::Edge {
        const Graph* _graph;
    public:
        
        EdgeIt() { }
        
        EdgeIt(lemon::Invalid i) : Edge(i) { }
        
        explicit EdgeIt(const Graph& graph) : _graph(&graph) {
            _graph->first(static_cast<Edge&>(*this));
        }
        
        EdgeIt(const Graph& graph, const Edge& edge) :
        Edge(edge), _graph(&graph) { }
        
        EdgeIt& operator++() {
            _graph->next(*this);
            return *this;
        }
        
    };
    
    class IncEdgeIt : public Parent::Edge {
        friend class BaseGraphExtender;
        const Graph* _graph;
        bool _direction;
    public:
        
        IncEdgeIt() { }
        
        IncEdgeIt(lemon::Invalid i) : Edge(i), _direction(false) { }
        
        IncEdgeIt(const Graph& graph, const Node &node) : _graph(&graph) {
            _graph->firstInc(*this, _direction, node);
        }
        
        IncEdgeIt(const Graph& graph, const Edge &edge, const Node &node)
        : _graph(&graph), Edge(edge) {
            _direction = (_graph->source(edge) == node);
        }
        
        IncEdgeIt& operator++() {
            _graph->nextInc(*this, _direction);
            return *this;
        }
    };
    
    // \brief Base node of the iterator
    //
    // Returns the base node (ie. the source in this case) of the iterator
    Node baseNode(const OutArcIt &arc) const {
        return Parent::source(static_cast<const Arc&>(arc));
    }
    // \brief Running node of the iterator
    //
    // Returns the running node (ie. the target in this case) of the
    // iterator
    Node runningNode(const OutArcIt &arc) const {
        return Parent::target(static_cast<const Arc&>(arc));
    }
    
    // \brief Base node of the iterator
    //
    // Returns the base node (ie. the target in this case) of the iterator
    Node baseNode(const InArcIt &arc) const {
        return Parent::target(static_cast<const Arc&>(arc));
    }
    // \brief Running node of the iterator
    //
    // Returns the running node (ie. the source in this case) of the
    // iterator
    Node runningNode(const InArcIt &arc) const {
        return Parent::source(static_cast<const Arc&>(arc));
    }
    
    // Base node of the iterator
    //
    // Returns the base node of the iterator
    Node baseNode(const IncEdgeIt &edge) const {
        return edge._direction ? this->u(edge) : this->v(edge);
    }
    // Running node of the iterator
    //
    // Returns the running node of the iterator
    Node runningNode(const IncEdgeIt &edge) const {
        return edge._direction ? this->v(edge) : this->u(edge);
    }
    
    // Mappable extension
    
    template <typename _Value>
    class NodeMap
    : public lemon::MapExtender<lemon::DefaultMap<Graph, Node, _Value> > {
        typedef lemon::MapExtender<lemon::DefaultMap<Graph, Node, _Value> > Parent;
        
    public:
        explicit NodeMap(const Graph& graph)
        : Parent(graph) {}
        NodeMap(const Graph& graph, const _Value& value)
        : Parent(graph, value) {}
        
    private:
        NodeMap& operator=(const NodeMap& cmap) {
            return operator=<NodeMap>(cmap);
        }
        
        template <typename CMap>
        NodeMap& operator=(const CMap& cmap) {
            Parent::operator=(cmap);
            return *this;
        }
        
    };
    
    template <typename _Value>
    class ArcMap
    : public lemon::MapExtender<lemon::DefaultMap<Graph, Arc, _Value> > {
        typedef lemon::MapExtender<lemon::DefaultMap<Graph, Arc, _Value> > Parent;
        
    public:
        explicit ArcMap(const Graph& graph)
        : Parent(graph) {}
        ArcMap(const Graph& graph, const _Value& value)
        : Parent(graph, value) {}
        
    private:
        ArcMap& operator=(const ArcMap& cmap) {
            return operator=<ArcMap>(cmap);
        }
        
        template <typename CMap>
        ArcMap& operator=(const CMap& cmap) {
            Parent::operator=(cmap);
            return *this;
        }
    };
    
    
    template <typename _Value>
    class EdgeMap
    : public lemon::MapExtender<lemon::DefaultMap<Graph, Edge, _Value> > {
        typedef lemon::MapExtender<lemon::DefaultMap<Graph, Edge, _Value> > Parent;
        
    public:
        explicit EdgeMap(const Graph& graph)
        : Parent(graph) {}
        
        EdgeMap(const Graph& graph, const _Value& value)
        : Parent(graph, value) {}
        
    private:
        EdgeMap& operator=(const EdgeMap& cmap) {
            return operator=<EdgeMap>(cmap);
        }
        
        template <typename CMap>
        EdgeMap& operator=(const CMap& cmap) {
            Parent::operator=(cmap);
            return *this;
        }
        
    };
    
    // Alteration extension
    
//    Node addNode() {
//        Node node = Parent::addNode();
//        notifier(Node()).add(node);
//        return node;
//    }
//    
//    Edge addEdge(const Node& from, const Node& to) {
//        Edge edge = Parent::addEdge(from, to);
//        notifier(Edge()).add(edge);
//        std::vector<Arc> ev;
//        ev.push_back(Parent::direct(edge, true));
//        ev.push_back(Parent::direct(edge, false));
//        notifier(Arc()).add(ev);
//        return edge;
//    }
    
//    void clear() {
//        notifier(Arc()).clear();
//        notifier(Edge()).clear();
//        notifier(Node()).clear();
//        Parent::clear();
//    }
    
    template <typename Graph, typename NodeRefMap, typename EdgeRefMap>
    void build(const Graph& graph, NodeRefMap& nodeRef,
               EdgeRefMap& edgeRef) {
        //Parent::build(graph, nodeRef, edgeRef);
        notifier(Node()).build();
        notifier(Edge()).build();
        notifier(Arc()).build();
    }
    
//    void erase(const Node& node) {
//        Arc arc;
//        Parent::firstOut(arc, node);
//        while (arc != INVALID ) {
//            erase(arc);
//            Parent::firstOut(arc, node);
//        }
//        
//        Parent::firstIn(arc, node);
//        while (arc != INVALID ) {
//            erase(arc);
//            Parent::firstIn(arc, node);
//        }
//        
//        notifier(Node()).erase(node);
//        Parent::erase(node);
//    }
//    
//    void erase(const Edge& edge) {
//        std::vector<Arc> av;
//        av.push_back(Parent::direct(edge, true));
//        av.push_back(Parent::direct(edge, false));
//        notifier(Arc()).erase(av);
//        notifier(Edge()).erase(edge);
//        Parent::erase(edge);
//    }
    
    BaseGraphExtender(int n) : Parent(n) {
        node_notifier.setContainer(*this);
        arc_notifier.setContainer(*this);
        edge_notifier.setContainer(*this);
    }
    
    ~BaseGraphExtender() {
        edge_notifier.clear();
        arc_notifier.clear();
        node_notifier.clear();
    }
};

//typedef BaseGraphExtender<BaseGraph<int,int>> ExtendedBaseGraph;
//typedef BaseGraphExtender<BaseGraph> ExtendedBaseGraph;
typedef BaseGraphExtender ExtendedBaseGraph;

class Graph : public ExtendedBaseGraph
{
    typedef ExtendedBaseGraph Parent;
    
    // Graphs require a size
    Graph() : ExtendedBaseGraph(0) {}
    /// Graphs are \e not copy constructible. Use GraphCopy instead.
    Graph(const Graph &) : ExtendedBaseGraph(0)  {}
    /// \brief Assignment of a graph to another one is \e not allowed.
    /// Use GraphCopy instead.
    void operator=(const Graph& ) {}
    
public:
    
    typedef lemon::True UndirectedTag;
    
    typedef typename Parent::Node Node;
    typedef typename Parent::Arc Arc;
    typedef typename Parent::Edge Edge;
    
    Graph( int n ) : Parent(n) {}
    
    template <typename _Graph>
    struct Constraints {
        void constraints() {
            lemon::checkConcept<lemon::concepts::BaseGraphComponent, _Graph>();
            lemon::checkConcept<lemon::concepts::IterableGraphComponent<>, _Graph>();
            lemon::checkConcept<lemon::concepts::IDableGraphComponent<>, _Graph>();
            lemon::checkConcept<lemon::concepts::MappableGraphComponent<>, _Graph>();
        }
    };
};

#endif /* newgraph_h */
