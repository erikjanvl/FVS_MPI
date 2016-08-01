//
//  inplace.hpp
//  back
//
//  Created by Erik Jan van Leeuwen on 25/7/16.
//  Copyright © 2016 Erik Jan van Leeuwen. All rights reserved.
//

#ifndef inplace_h
#define inplace_h

#include <iostream>
#include <vector>
#include <queue>
#include <math.h>

#define INPLACE_HAS_BIT( number, bit ) ( ((number) & (1 << (bit))) == (1 << (bit)) )
//#define INPLACE_HAS_BIT( number, bit ) false

// provides a square matrix with fixed size length
// a protected version enables a layered version
template <class value_type>
class inplace_square_matrix {
    std::vector<value_type> mat;
    int side_length;
    
protected:
    inplace_square_matrix() : mat(0) {}
    inplace_square_matrix(const inplace_square_matrix&) : mat(0) {}
    
    typename std::vector<value_type>::reference realat( int row, int column ) {
        if( column >= side_length ) {
            throw std::out_of_range("At out of range");
        }
        //return const_cast<std::vector<value_type>::reference>(mat.at(row*side_length+column+layer));
        return mat.at(row*side_length+column);
    }
    
    typename std::vector<value_type>::const_reference realat( int row, int column ) const {
        if( row >= side_length || column >= side_length ) {
            throw std::out_of_range("At out of range");
        }
        //return const_cast<std::vector<value_type>::const_reference>(mat.at(row*side_length+column+layer));
        return mat.at(row*side_length+column);
    }
    
public:
    inplace_square_matrix( int n ) : mat(n * n), side_length(n) {}
    inplace_square_matrix( int n, const value_type &initial ) : mat(n * n, initial), side_length(n) {}
    
    value_type &at( int row, int column ) const {
        if( row >= side_length || column >= side_length ) {
            throw std::out_of_range("At out of range");
        }
        return const_cast<value_type&>(mat.at(row*side_length+column));
    }
};

template <class value_type>
class inplace_layered_square_matrix {
    std::vector<value_type> mat;
    int side_length, layers;
    
protected:
    inplace_layered_square_matrix() : mat(0) {}
    inplace_layered_square_matrix(const inplace_layered_square_matrix&) : mat(0) {}
    
    typename std::vector<value_type>::reference realat( int row, int column, int layer ) {
        if( row < 0 || column < 0 || layer < 0 || row >= side_length || column >= side_length || layer >= layers ) {
            throw std::out_of_range("At out of range");
        }
        //return const_cast<std::vector<value_type>::reference>(mat.at(row*side_length+column+layer));
        return mat.at(layers*(row*side_length+column)+layer);
    }
    
    /*typename std::vector<value_type>::const_reference realat( int row, int column, int layer ) const {
        if( row < 0 || column < 0 || layer < 0 || row >= side_length || column >= side_length || layer >= layers ) {
            throw std::out_of_range("At out of range");
        }
        //return const_cast<std::vector<value_type>::const_reference>(mat.at(row*side_length+column+layer));
        return mat.at(layers*(row*side_length+column)+layer);
    }*/
    
public:
    inplace_layered_square_matrix( int n, int _layers, const value_type& initial ) : mat( _layers * n * n, initial ), side_length(n), layers(_layers) {}
    
    const value_type &at( int row, int column, int layer ) const {
        if( row < 0 || column < 0 || layer < 0 || row >= side_length || column >= side_length || layer >= layers ) {
            throw std::out_of_range("At out of range");
        }
        return static_cast<const value_type&>(mat.at(layers*(row*side_length+column)+layer));
    }
    
    value_type &at( int row, int column, int layer ) {
        if( row < 0 || column < 0 || layer < 0 || row >= side_length || column >= side_length || layer >= layers ) {
            throw std::out_of_range("At out of range");
        }
        return const_cast<value_type&>(static_cast<const value_type&>(mat.at(layers*(row*side_length+column)+layer)));
    }
};

// a special square matrix that only allows values between 0 and 3
class inplace_square_matrix_4 {//: protected inplace_layered_square_matrix<bool> {
    //typedef inplace_layered_square_matrix<bool> Parent;
    std::vector<bool> mat;
    int side_length;

protected:
    inplace_square_matrix_4() {}
    inplace_square_matrix_4(const inplace_square_matrix_4&) {}
    
    int elem( int row, int column, int layer ) const {
        return 2*(row*side_length+column)+layer;
    }

public:
    enum value {
        ZERO = 0,
        ONE = 1,
        TWO = 2,
        THREE = 3
    };
    inplace_square_matrix_4( int n ) : side_length(n), mat(2*n*n, false) {}
    
    void set( int row, int column, value v ) {
        switch(v) {
            case ZERO:
                mat.at(elem(row,column,0)) = false;
                mat.at(elem(row,column,1)) = false;
                break;
            case ONE:
                mat.at(elem(row,column,0)) = true;
                mat.at(elem(row,column,1)) = false;
                break;
            case TWO:
                mat.at(elem(row,column,0)) = false;
                mat.at(elem(row,column,1)) = true;
                break;
            case THREE:
                mat.at(elem(row,column,0)) = true;
                mat.at(elem(row,column,1)) = true;
                break;
        }
    }
    
    /*value increase( int row, int column ) {
        if( !realat(row,column,0) || !realat(row,column,1) ) {
            if( realat(row,column,0) ) {
                realat(row,column,1) = !realat(row,column,1);
            }
            realat(row,column,0) = !realat(row,column,0);
        }
        return get(row,column);
    }
    
    value decrease( int row, int column ) {
        if( realat(row,column,0) || realat(row,column,1) ) {
            if( !realat(row,column,0) ) {
                realat(row,column,1) = !realat(row,column,1);
            }
            realat(row,column,0) = !realat(row,column,0);
        }
        return get(row,column);
    }*/
    
//    value at( int row, int column ) const {
//        return get(row, column);
//    }
    
    value get( int row, int column ) const {
        if( mat.at(elem(row,column,1)) ) {
            if( mat.at(elem(row,column,0)) ) {
                return value::THREE;
            } else {
                return value::TWO;
            }
        } else {
            if( mat.at(elem(row,column,0)) ) {
                return value::ONE;
            } else {
                return value::ZERO;
            }
        }
    }
};

template <class value_type, class size_type = int>
class inplace_layered_arraylist {
    std::vector<value_type> elems;
    std::vector<size_type> sizes;
    int layers;
    size_type max_size_per_layer;
    
    inplace_layered_arraylist() : inplace_layered_arraylist(0,0) {}
    inplace_layered_arraylist( const inplace_layered_arraylist& ) : inplace_layered_arraylist(0,0) {}
    
protected:
    int elem( int layer, size_type pos ) const {
        return layer * max_size(layer) + pos;
    }
    
    value_type* pointInternal( int layer, size_type pos ) const {
        return const_cast<value_type*>(&(elems.at(elem(layer,pos))));
    }
    
    value_type& atInternal( int layer, size_type pos ) const {
        return const_cast<value_type&>(elems.at(elem(layer,pos)));
    }
    
public:
    // _layers is the number of layers
    // _max_size is the maximum size of each layer
    inplace_layered_arraylist( int _layers, int _max_size_per_layer ) : elems(_layers * _max_size_per_layer), sizes(_layers,0), layers(_layers), max_size_per_layer(_max_size_per_layer) {}
    
    int no_of_layers() const {
        return layers;
    }
    
    size_type size( int layer ) const {
        return sizes[layer];
    }
    
    bool empty( int layer ) const {
        return size(layer) == 0;
    }
    
    size_type max_size( int layer ) const {
        return max_size_per_layer;
    }
    
    size_type capacity( int layer ) const {
        return max_size(layer);
    }
    
    value_type* point( int layer, size_type pos ) const {
        if( pos < 0 || pos >= size(layer) ) {
            throw std::out_of_range("Point to invalid object");
        }
        return const_cast<value_type*>(&(elems.at(elem(layer,pos))));
    }
    
    value_type& at( int layer, size_type pos ) const {
        if( pos >= size(layer) ) {
            throw std::out_of_range("At to invalid object");
        }
        return const_cast<value_type&>(elems.at(elem(layer,pos)));
    }
    
    void clear() {
        std::fill( sizes.begin(), sizes.end(), 0 );
    }
    
    void clear( int layer ) {
        sizes[layer] = 0;
    }
    
    void clear( int layer, size_t pos ) {
        if( pos < 0 || pos >= size(layer) ) {
            throw std::out_of_range("Clearing out of range object");
        } else if( pos != size(layer) - 1 ) {
            elems[elem(layer,pos)] = elems[elem(layer,size(layer)-1)];
        }
        sizes[layer]--;
    }
    
    void push_back( int layer, const value_type& val ) {
        if( size(layer) == max_size(layer) ) {
            throw std::overflow_error("Layered list overflowed");
        }
        elems[elem(layer,size(layer))] = val;
        sizes[layer]++;
    }
};

template <class value_type, class size_type = int>
class inplace_arraylist {
    inplace_layered_arraylist<value_type,size_type> array;
    inplace_arraylist() : array(0,0) {}
    inplace_arraylist( const inplace_arraylist& ) : array(0,0) {}
    
public:
    // _max_size is the maximum size of each layer
    inplace_arraylist( int _max_size_per_layer ) : array(1, _max_size_per_layer) {}
    
    size_type size() const {
        return array.size(0);
    }
    
    bool empty() const {
        return array.empty(0);
    }
    
    size_type max_size() const {
        return array.max_size(0);
    }
    
    size_type capacity() const {
        return array.capacity(0);
    }
    
    value_type* point( size_type pos ) const {
        return array.point(0,pos);
    }
    
    value_type& at( size_type pos ) const {
        return array.at(0,pos);
    }
    
    void clear() {
        array.clear();
    }
    
    void clear( size_t pos ) {
        array.clear(0,pos);
    }
    
    void push_back( const value_type& val ) {
        array.push_back(0,val);
    }
};

template <class value_type, class size_type = int>
class inplace_fifo {
    std::vector<value_type> array;
    size_type start, end, cap;
public:
    inplace_fifo( size_type c ) : array(c), start(0), end(0), cap(c) { }
    
    const value_type& front() const {
        if( start == -1 ) {
            throw std::out_of_range("Asking front of empty FIFO");
        }
        return array.at(start);
    }
    
    void push( const value_type& val ) {
        if( size() == cap ) {
            throw std::overflow_error("FIFO overflowed");
        }
        if( start == -1 ) {
            start = 0;
            end = 0;
        } else if( end == cap - 1 ) {
            end = 0;
        } else {
            end++;
        }
        array.at(end) = val;
    }
    
    void pop() {
        size_type s = size();
        if( s > 1 ) {
            if( start == cap - 1 ) {
                start = 0;
            } else {
                start++;
            }
        } else {
            start = -1;
        }
    }
    
    void empty() const {
        return start == -1;
    }
    
    size_type size() const {
        if( start == -1 ) {
            return 0;
        } else if( start <= end ) {
            return end-start+1;
        } else {
            return (end+1) + (cap-end+1);
        }
    }
    
    const value_type& at( size_type pos ) const {
        if( pos < 0 || pos >= size() ) {
            throw std::out_of_range("Asking FIFO non-existing element");
        }
        pos = pos+start;
        if( pos >= cap ) {
            pos -= cap;
        }
        return array.at(pos);
    }
    
    void clear() {
        start = -1;
    }
};


// an inplace stack
// value_type must support operator=
template <class value_type, class size_type = int>
class inplace_stack {
    inplace_arraylist<value_type,size_type> array;
    
public:
    inplace_stack( size_type c ) : array(c) {}
    
    void push( const value_type& val ) {
        if( array.size() == array.capacity() ) {
            throw std::out_of_range("inplace_stack overflowed");
        }
        array.push_back(val);
    }
    
    value_type& top() const {
        return array.at( array.size()-1 );
    }
    
    value_type& almost_top() const {
        return array.at( array.size()-2 );
    }
    
    void pop() {
        array.clear( array.size()-1 );
    }
    
    bool empty() const {
        return array.empty();
    }
    
    size_type size() const {
        return array.size();
    }
    
    size_type max_size() const {
        return array.max_size();
    }
    
    size_type capacity() const {
        return array.capacity();
    }
    
    void clear() {
        array.clear();
    }
    
    const value_type& at( size_type pos ) const {
        return array.at(pos);
    }
};



// this provides a linked-list in a fixed size array (in-place)
// it provides constant-time access to a specific element, as well as constant-time access to previous and next elements
template <class value_type, class size_type = int>
class inplace_list {
    
public:
    const static int Invalid = -1;
    
protected:
    // each element of the array will have a previous and next, and a value
    class ListElem {
        friend class inplace_list;
        value_type val;
        size_type prev = Invalid;
        size_type next = Invalid;
        
        void clear() {
            prev = next = Invalid;
        }
    };
    
    std::vector<ListElem> elems;
    size_type numberOfElems = 0;
    size_type maxNumberOfElems;
    size_type firstElem = Invalid;
    
    // internal helper function
    ListElem *pointInternal( size_type i ) const {
        return const_cast<ListElem*>(&(elems.at(i)));
    }
    
public:
    
    // a vectorlist can only be constructed with a given (fixed) size
    inplace_list(size_type c) : elems(c), maxNumberOfElems(c) {}
    
    size_type size() const {
        return numberOfElems;
    }
    
    bool empty() const {
        return size() == 0;
    }
    
    //    bool has_free() const {
    //        return size() < max_size();
    //    }
    
    size_type max_size() const {
        return maxNumberOfElems;
    }
    
    size_type capacity() const {
        return maxNumberOfElems;
    }
    
    value_type* point( size_type i ) const {
        if( !is_valid(i) ) {
            throw std::out_of_range("Point to invalid object");
        }
        return const_cast<value_type*>(&(elems.at(i).val));
    }
    
    value_type& at( size_type i ) const {
        if( !is_valid(i) ) {
            throw std::out_of_range("At to invalid object");
        }
        return const_cast<value_type&>(elems.at(i).val);
    }
    
    //    const value_type& operator[]( size_type i ) const {
    //        return elems[i].val;
    //    }
    
    value_type& front() const {
        return at(firstElem);
    }
    
    void pop_front() {
        disable( first() );
    }
    
    // enable the element at position i
    // returns false only if i is not a valid index
    // returns true if succesfull or position i is already active
    // WARNING: worst-case it runs in time linear in the size of the list
    bool enable( size_type i ) {
        // first check whether i is a sane number
        // then check whether i is already
        if( !is_valid(i,false) ) {
            return false;
        }
        if( is_valid(i) ) {
            return true;
        }
        
        if( size() == 0 ) {
            // first element, so no neighbors in the list
            return enable( i, Invalid, Invalid );
        }
        if( i <= first() ) {
            // i comes before the first element, so prev is Invalid
            return enable( i, Invalid, first() );
        }
        // now we know that we have at least one element, and first() < i
        // look for the two elements that enclose i
        size_type prev_i = first();
        size_type next_i = first();
        while( next_i != Invalid && next_i <= i ) {
            prev_i = next_i;
            next_i = next(i);
        }
        return enable( i, prev_i, next_i );
    }
    
    // enable the element at position i
    // prev_i indicates the previous element
    // next_i indicates the next element
    // returns true iff the element was successfully activated
    bool enable( size_type i, size_type prev_i, size_type next_i ) {
        // first check that i is a sane number
        // then check whether i is already active
        if( !is_valid(i,false) ) {
            return false;
        }
        if( is_valid(i) ) {
            return true;
        }
        ListElem* elem = pointInternal(i);
        elem->clear();
        
        if( prev_i == Invalid ) {
            // if prev_i is Invalid, then we try to insert a new first element
            if( next_i == Invalid ) {
                // if next_i is also Invalid, then the list should be empty now
                if( size() != 0 ) {
                    return false;
                }
                firstElem = i;
                numberOfElems++;
                return true;
            }
            // next_i is valid, and since prev_i is Invalid, it should be the first element
            // also add a sanity check that next_i comes after i
            if( next_i != first() || !(i < next_i) ) {
                return false;
            }
            // new first element
            ListElem* next = pointInternal(next_i);
            elem->next = next_i;
            next->prev = i;
            firstElem = i;
            numberOfElems++;
            return true;
        }
        
        // now prev_i should be valid
        if( !is_valid(prev_i) || !(prev_i < i) ) {
            return false;
        }
        ListElem* prev = pointInternal(prev_i);
        
        if( next_i == Invalid ) {
            // if next_i is Invalid, then we try to insert a new last element
            // make sure that prev_i indeed has no next
            if( next(prev_i) != Invalid ) {
                return false;
            }
            prev->next = i;
            elem->prev = prev_i;
            numberOfElems++;
            return true;
        }
        
        // now next_i should be valid
        if( !is_valid(next_i) ) {
            return false;
        }
        ListElem* next = pointInternal(next_i);
        
        // now make sure that prev and next are indeed consecutive
        if( prev->next != next_i || next->prev != prev_i ) {
            return false;
        }
        elem->prev = prev_i;
        elem->next = next_i;
        prev->next = i;
        next->prev = i;
        numberOfElems++;
        return true;
    }
    
    // disable the element at position i
    // same as removing it
    void disable( size_type i ) {
        remove_at(i);
    }
    
    // remove the element at position i
    void remove_at( size_type i ) {
        if( is_valid(i) ) {
            // i is active. Now fix the linked list
            ListElem* elem = pointInternal(i);
            if( elem->next != Invalid ) {
                ListElem* next = pointInternal(elem->next);
                if( elem->prev != Invalid ) {
                    next->prev = elem->prev;
                } else {
                    //no previous element. elem must be the first element, so now elem->next is the new first element!
                    firstElem = elem->next;
                }
            } else if( elem->prev == Invalid ) {
                // both elem->next and elem->prev are invalid
                // then elem is the first element
                // so clear first elem
                firstElem = Invalid;
            }
            if( elem->prev != Invalid ) {
                ListElem* prev = pointInternal(elem->prev);
                prev->next = elem->next;
            }
            elem->clear();
            numberOfElems--;
        }
    }
    
    size_type first() const {
        return firstElem;
    }
    
    size_type next( size_type i ) const {
        if( !is_valid(i) ) {
            return Invalid;
        }
        ListElem* n = pointInternal(i);
        return n->next;
    }
    
    size_type prev( size_type i ) const {
        if( !is_valid(i) ) {
            return Invalid;
        }
        ListElem* n = pointInternal(i);
        return n->prev;
    }
    
    void clear() {
        size_type i = first();
        while( i != Invalid ) {
            size_type j = next(i);
            pointInternal(i)->clear();
            i = j;
        }
    }
    
    inline bool is_valid(size_type i, bool fail_on_cleared = true) const {
        if( i < 0 || i >= max_size() ) {
            return false;
        } else if( i == firstElem ) {
            return true;
        } else if( fail_on_cleared && elems.at(i).next == Invalid && elems.at(i).prev == Invalid ) {
            return false;
        } else {
            return true;
        }
    }
};

/*// an inplace stack
// value_type must support operator=
template <class value_type, class size_type = int>
class inplace_stack {
    inplace_list<value_type,size_type> list;
    
public:
    inplace_stack( size_type c ) : list(c) {}
    
    void push( const value_type& val ) {
        if( list.size() == list.capacity() ) {
            throw std::out_of_range("inplace_stack overflowed");
        }
        size_type pos = list.size();
        if( pos == 0 ) {
            list.enable( pos, list.Invalid, list.Invalid );
        } else {
            list.enable( pos, pos-1, list.Invalid );
        }
        list.at(pos) = val;
    }
    
    value_type& top() const {
        return list.at( list.size()-1 );
    }
    
    void pop() {
        list.disable( list.size()-1 );
    }
    
    size_type size() const {
        return list.size();
    }
    
    size_type max_size() const {
        return list.max_size();
    }
    
    size_type capacity() const {
        return list.capacity();
    }
};*/

/*
// an inplace linked list that remembers its open spaces
// gives constant time insertion and deletion
template <class value_type>
class inplace_list_with_memory : inplace_list<value_type> {

    
public:
    //<#member functions#>
};


template <class value_type>
class inplace_matrix {
    
public:
    const static int Invalid = -1;
    
    
protected:
    class ListElem {
        friend class matrixlist;
        value_type val;
        int prev = Invalid;
        int next = Invalid;
        
        void clear() {
            prev = next = Invalid;
        }
    };
    
    class RowInfo {
        friend class matrixlist;
        int firstElem = Invalid;
    };
    
    std::vector<ListElem> elems;
    int numberOfElems = 0;
    int maxNumberOfElems;
    inplace_list<RowInfo> rows;
    int numberOfRows = 0;
    int maxNumberOfRows;
    int firstRow = Invalid;
    
    
    ListElem *pointInternal( int i ) {
        return const_cast<ListElem*>(&(elems.at(i)));
    }
    
public:
    
    inplace_matrix(int c) : elems(c*c), maxNumberOfElems(c*c), rows(c), maxNumberOfRows(c) {}
    
    int size() const {
        return numberOfElems;
    }
    
    int size_row() const {
        return numberOfRows;
    }
    
    bool empty() const {
        return size() == 0;
    }
    
    int max_size() const {
        return maxNumberOfElems;
    }
    
    int max_size_row() const {
        return maxNumberOfRows;
    }
    
    int capacity() const {
        return maxNumberOfElems;
    }
    
    int capacity_row() const {
        return maxNumberOfRows;
    }
    
    bool has_free() const {
        return size() < max_size();
    }
    
    bool has_free_row() const {
        return size_row() < max_size_row();
    }
    
    value_type* point( int i ) const {
        if( !is_valid(i) ) {
            throw std::out_of_range("Pointing to invalid object");
        }
        return const_cast<value_type*>(&(elems.at(i).val));
    }
    
    const value_type& at( int i ) const {
        if( !is_valid(i) ) {
            throw std::out_of_range("At invalid object");
        }
        return elems.at(i).val;
    }
    
    void clear() {
        int i = first();
        while( i != Invalid ) {
            int j = next(i);
            pointInternal(i)->clear();
            i = j;
        }
    }
    
    void remove_row_at( int row ) {
        if( is_valid_row(row) ) {
            int i = first_in_row(row);
            while( i != Invalid ) {
                int j = next_in_row(i);
                pointInternal(i)->clear();
                i = j;
            }
            rows.remove_at(row);
        }
    }
    
    void remove_at( int row, int column ) {
        remove_at(row * maxNumberOfRows + column);
    }
    
    void remove_at( int i ) {
        if( is_valid(i) ) {
            ListElem* elem = pointInternal(i);
            if( elem->next != Invalid ) {
                ListElem* next = pointInternal(elem->next);
                if( elem->prev != Invalid ) {
                    next->prev = elem->prev;
                } else {
                    //no previous element. next is the new first element!
                    int row = floor( i / maxNumberOfRows );
                    rows.at(row).firstElem = elem->next;
                }
            }
            if( elem->prev != Invalid ) {
                ListElem* prev = pointInternal(elem->prev);
                prev->next = elem->next;
            }
            numberOfElems--;
        }
    }
    
    int first_row() const {
        return rows.first();
    }
    
    int first_in_row( int row ) const {
        if( !rows.is_valid(row) ) {
            return Invalid;
        } else {
            return rows.at(row).firstElem;
        }
    }
    
    int first() const {
        return first_in_row(first_row());
    }
    
    int next_row( int i ) const {
        if( !is_valid_row() ) {
            return false;
        }
        return rows.next(i);
    }
    
    int next_in_row( int row, int column ) const {
        return next_in_row( row * maxNumberOfRows + column );
    }
    
    int next_in_row( int i ) const {
        if( !is_valid(i) ) {
            return Invalid;
        }
        ListElem* n = pointInternal(i);
        return n->next;
    }
    
    int next( int i ) const {
        if( !is_valid(i) ) {
            return Invalid;
        }
        
        int j = next_in_row(i);
        int row = floor(i / maxNumberOfRows);
        while( j == Invalid ) {
            row = next_row(row);
            if( row == Invalid ) {
                return Invalid;
            } else {
                j = first_in_row(row);
            }
        }
        return j;
    }
    
    bool is_valid(int i) const {
        if( i < 0 || i >= max_size() ) {
            return false;
        }
        int row = floor(i / maxNumberOfRows);
        int column = i % maxNumberOfRows;
        if( !is_valid_row(row) ) {
            return false;
        } else if( column == rows.at(row).firstElem ) {
            return true;
        } else if( elems.at(i).next == Invalid && elems.at(i).prev == Invalid ) {
            return false;
        } else {
            return true;
        }
    }
    
    bool is_valid_row(int i) const {
        return rows.isValid(i);
    }
    
};*/


#endif /* inplace_h */
