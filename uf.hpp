/*
 *  uf.hpp
 *  Created by Erik Jan van Leeuwen on 14/7/16.
 *
 * This file re-uses a large quantity of code that is a part of LEMON, a generic C++ optimization library.
 * Copyright (C) 2003-2013
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 */

#ifndef uf_h
#define uf_h

/// \ingroup auxdat
///
/// \brief A \e Union-Find data structure implementation
///
/// The class implements the \e Union-Find data structure.
/// The union operation uses rank heuristic, while
/// the find operation uses path compression.
/// This is a very simple but efficient implementation, providing
/// only four methods: join (union), find, insert and size.
/// For more features, see the \ref UnionFindEnum class.
///
/// It is primarily used in Kruskal algorithm for finding minimal
/// cost spanning tree in a graph.
/// \sa kruskal()
///
/// \pre You need to add all the elements by the \ref insert()
/// method.
template <class Item, class size_type=int>
class UnionFind {
   
private:
    // If the items vector stores negative value for an item then
    // that item is root item and it has -items[it] component size.
    // Else the items[it] contains the index of the parent.
    std::vector<size_type> items;
    size_type classCount;
    
    bool rep(size_type idx) const {
        return items[idx] < 0;
    }
    
    int repIndex(size_type idx) const {
        size_type k = idx;
        // find the item
        while (!rep(k)) {
            k = items[k] ;
        }
        // now do path-compression
        while (idx != k) {
            size_type next = items[idx];
            const_cast<size_type&>(items[idx]) = k;
            idx = next;
        }
        return k;
    }
    
public:
    
    /// \brief Constructor
    ///
    /// Constructor of the UnionFind class. You should give an item to
    /// integer map which will be used from the data structure. If you
    /// modify directly this map that may cause segmentation fault,
    /// invalid data structure, or infinite loop when you use again
    /// the union-find.
    UnionFind(size_type c) : items(c,-1), classCount(c) {}
    
    /// \brief Returns the index of the element's component.
    ///
    /// The method returns the index of the element's component.
    /// This is an integer between zero and the number of inserted elements.
    ///
    size_type find(const Item& a) {
        return repIndex(a.id);
    }
    
    /// \brief Clears the union-find data structure
    ///
    /// Erase each item from the data structure.
    void clear() {
        std::fill(items.begin(),items.end(),-1);
        classCount = items.size();
    }
    
    /// \brief Inserts a new element into the structure.
    ///
    /// This method inserts a new element into the data structure.
    ///
    /// The method returns the index of the new component.
//    void insert(const Item& a) {
//        items[a.id] = a.id;
//    }
    
    /// \brief Joining the components of element \e a and element \e b.
    ///
    /// This is the \e union operation of the Union-Find structure.
    /// Joins the component of element \e a and component of
    /// element \e b. If \e a and \e b are in the same component then
    /// it returns false otherwise it returns true.
    bool join(const Item& a, const Item& b) {
        size_type ka = repIndex(a.id);
        size_type kb = repIndex(b.id);
        
        if ( ka == kb )
            return false;
        
        if (items[ka] < items[kb]) {
            // ka is bigger
            items[ka] += items[kb];
            items[kb] = ka;
        } else {
            // kb is bigger
            items[kb] += items[ka];
            items[ka] = kb;
        }
        classCount--;
        return true;
    }
    
    size_type numberOfClasses() const {
        return classCount;
    }
    
    /// \brief Returns the size of the component of element \e a.
    ///
    /// Returns the size of the component of element \e a.
    size_type size(const Item& a) {
        size_type k = repIndex(a.id);
        return - items[k];
    }
    
};

#endif /* uf_h */
