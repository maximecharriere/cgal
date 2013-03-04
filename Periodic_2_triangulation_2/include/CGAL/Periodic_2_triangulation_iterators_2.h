// Copyright (c) 1999,2008-2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_ITERATORS_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_ITERATORS_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/iterator.h>

namespace CGAL {

template < class T >
class Periodic_2_triangulation_triangle_iterator_2 {
  // Iterates over the primitives in a periodic triangulation.
  // Options:
  // - STORED: output each primitive from the Tds exactly once
  // - UNIQUE: output exactly one periodic copy of each primitive, no matter
  //     whether the current tds stores a n-sheeted covering for n!=1.
  // - STORED_COVER_DOMAIN: output each primitive whose intersection with the
  //     actually used periodic domain is non-zero.
  // - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
  //     with the original domain that the user has given is non-zero
  //
  // Comments:
  // When computing in 1-sheeted covering, there will be no difference in the
  // result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
  // UNIQUE_COVER_DOMAIN.
  
public:
  
  typedef typename T::Periodic_triangle                   value_type;
  typedef const typename T::Periodic_triangle *           pointer;
  typedef const typename T::Periodic_triangle &           reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;
  
  typedef typename T::Periodic_triangle                   Periodic_triangle;
  typedef Periodic_2_triangulation_triangle_iterator_2<T> Periodic_triangle_iterator;
  typedef typename T::Face                                Face;
  typedef typename T::Face_iterator                       Face_iterator;
  
  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;
  
  Periodic_2_triangulation_triangle_iterator_2(Iterator_type it = T::STORED)
  : _t(NULL), _it(it), _off(0) {}
  
  Periodic_2_triangulation_triangle_iterator_2(const T * t,
                                               Iterator_type it = T::STORED)
  : _t(t), pos(_t->faces_begin()), _it(it), _off(0) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->faces_end() && !is_canonical() )
        ++pos;
    }
  }
  
  // used to initialize the past-the-end iterator
  Periodic_2_triangulation_triangle_iterator_2(const T* t, int,
                                                  Iterator_type it = T::STORED)
  : _t(t), pos(_t->faces_end()), _it(it), _off(0) {}
  
  Periodic_triangle_iterator& operator++() {
    switch (_it) {
      case T::STORED:
        ++pos;
        break;
      case T::UNIQUE:
        do { ++pos; } while (pos != _t->faces_end() && !is_canonical());
        break;
      case T::STORED_COVER_DOMAIN:
      case T::UNIQUE_COVER_DOMAIN:
        increment_domain();
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return *this;
  }
  
  Periodic_triangle_iterator& operator--() {
    switch (_it) {
      case T::STORED:
        --pos;
        break;
      case T::UNIQUE:
        do { --pos; } while (pos != _t->cells_begin() && !is_canonical());
        break;
      case T::STORED_COVER_DOMAIN:
      case T::UNIQUE_COVER_DOMAIN:
        decrement_domain();
    };
    return *this;
  }
  
  Periodic_triangle_iterator operator++(int)
  {
    Periodic_triangle_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
  
  Periodic_triangle_iterator operator--(int)
  {
    Periodic_triangle_iterator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Periodic_triangle_iterator& ti) const
  {
    // We are only allowed to compare iterators of the same type.
    CGAL_triangulation_assertion(_it == ti._it);
    return _t == ti._t && pos == ti.pos && _off == ti._off;
  }
  
  bool operator!=(const Periodic_triangle_iterator& ti) const
  {
    return !(*this == ti);
  }
  
  reference operator*() const
  {
    periodic_triangle = construct_periodic_triangle();
    return periodic_triangle;
  }
  
  pointer operator->() const
  {
    periodic_triangle = construct_periodic_triangle();
    return &periodic_triangle;
  }
  
  Face_iterator get_face() const
  {
    return pos;
  }
  
private:
  const T*  _t;
  Face_iterator pos; // current cell.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_triangle periodic_triangle; // current triangle.
  
private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1, off2;
    get_edge_offsets(off0, off1, off2);
    
    if (_t->number_of_sheets() != make_array(1,1)) {
      // If there is one offset with entries larger than 1 then we are
      // talking about a vertex that is too far away from the original
      // domain to belong to a canonical triangle.
      if (off0.x() > 1) return false;
      if (off0.y() > 1) return false;
      if (off1.x() > 1) return false;
      if (off1.y() > 1) return false;
      if (off2.x() > 1) return false;
      if (off2.y() > 1) return false;
    }
    
    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x() & off2.x();
    int offy = off0.y() & off1.y() & off2.y();
    
    return (offx == 0 && offy == 0);
  }
  
  // Artificial incrementation function that takes periodic
  // copies into account.
  void increment_domain() {
    int off = get_drawing_offsets();
    CGAL_triangulation_assertion(_off <= off);
    if (_off == off) {
      _off = 0;
      do { ++pos; } while (_it == T::UNIQUE_COVER_DOMAIN
                           && pos != _t->faces_end() && !is_canonical());
    } else {
      do {
        ++_off;
      } while ((((~_off)|off)&3)!=3); // Increment until a valid
                                      // offset has been found
    }
  }
  
  // Artificial decrementation function that takes periodic
  // copies into account.
  void decrement_domain() {
    if (_off == 0) {
      if (pos == _t->cells_begin()) return;
      do { --pos; } while (_it == T::UNIQUE_COVER_DOMAIN && !is_canonical());
      _off = get_drawing_offsets();
    } else {
      int off = get_drawing_offsets();
      do {
        --_off;
      } while ((((~_off)|off)&3)!=3); // Decrement until a valid
                                      // offset has been found
    }
  }
  
  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1,
                        Offset &off2) const {
    Offset face_off0 = _t->int_to_off(pos->offset(0));
    Offset face_off1 = _t->int_to_off(pos->offset(1));
    Offset face_off2 = _t->int_to_off(pos->offset(2));
    Offset diff_off((face_off0.x() == 1 
                     && face_off1.x() == 1 
                     && face_off2.x() == 1)?-1:0,
                    (face_off0.y() == 1 
                     && face_off1.y() == 1
                     && face_off2.y() == 1)?-1:0);
    off0 = _t->combine_offsets(_t->get_offset(pos,0), diff_off);
    off1 = _t->combine_offsets(_t->get_offset(pos,1), diff_off);
    off2 = _t->combine_offsets(_t->get_offset(pos,2), diff_off);
  }
  
  // return an integer that encodes the translations which have to be
  // applied to the edge *pos
  int get_drawing_offsets() {
    Offset off0, off1, off2;
    // Choose edges that are to be duplicated. These are edges that
    // intersect the boundary of the periodic domain. In UNIQUE mode
    // this means that the offset with respect to drawing should
    // differ in some entries. Otherwise we consider the offsets
    // internally stored inside the cell telling us that this cell
    // wraps around the domain.
    if (_it == T::UNIQUE_COVER_DOMAIN)
      get_edge_offsets(off0,off1,off2);
    else {
      CGAL_triangulation_assertion(_it == T::STORED_COVER_DOMAIN);
      off0 = _t->int_to_off(pos->offset(0));
      off1 = _t->int_to_off(pos->offset(1));
      off2 = _t->int_to_off(pos->offset(2));
    }
    
    CGAL_triangulation_assertion(off0.x() == 0 || off0.x() == 1);
    CGAL_triangulation_assertion(off0.y() == 0 || off0.y() == 1);
    CGAL_triangulation_assertion(off1.x() == 0 || off1.x() == 1);
    CGAL_triangulation_assertion(off1.y() == 0 || off1.y() == 1);
    CGAL_triangulation_assertion(off2.x() == 0 || off2.x() == 1);
    CGAL_triangulation_assertion(off2.y() == 0 || off2.y() == 1);
    
    int offx = ( ((off0.x() == 0 && off1.x() == 0 
                   && off2.x() == 0)
                  || (off0.x() == 1 && off1.x() == 1 
                      && off2.x() == 1)) ? 0 : 1);
    int offy = ( ((off0.y() == 0 && off1.y() == 0 
                   && off2.y() == 0)
                  || (off0.y() == 1 && off1.y() == 1 
                      && off2.y() == 1)) ? 0 : 1);
    
    return( 2*offx + offy );
  }
  
  Periodic_triangle construct_periodic_triangle() const {
    CGAL_triangulation_assertion(pos != typename T::Face_handle());
    Offset off0, off1, off2;
    get_edge_offsets(off0, off1, off2);
    Offset transl_off = Offset((((_off>>2)&1)==1 ? -1:0),
                               (((_off>>1)&1)==1 ? -1:0));
    if (_it == T::STORED_COVER_DOMAIN) {
      off0 = _t->combine_offsets(off0,transl_off);
      off1 = _t->combine_offsets(off1,transl_off);
      off2 = _t->combine_offsets(off2,transl_off);
    }
    if (_it == T::UNIQUE_COVER_DOMAIN) {
      off0 += transl_off;
      off1 += transl_off;
      off2 += transl_off;
    }
    return make_array(std::make_pair(pos->vertex(0)->point(),off0),
                      std::make_pair(pos->vertex(1)->point(),off1),
                      std::make_pair(pos->vertex(2)->point(),off2));
  }
};

template < class T >
class Periodic_2_triangulation_segment_iterator_2 {
  // Iterates over the primitives in a periodic triangulation.
  // Options:
  // - STORED: output each primitive from the Tds exactly once
  // - UNIQUE: output exactly one periodic copy of each primitive, no matter
  //     whether the current tds stores a n-sheeted covering for n!=1.
  // - STORED_COVER_DOMAIN: output each primitive whose intersection with the
  //     actually used periodic domain is non-zero.
  // - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
  //     with the original domain that the user has given is non-zero
  //
  // Comments:
  // When computing in 1-sheeted covering, there will be no difference in the
  // result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
  // UNIQUE_COVER_DOMAIN.
  
public:
  
  typedef typename T::Periodic_segment                    value_type;
  typedef const typename T::Periodic_segment *            pointer;
  typedef const typename T::Periodic_segment &            reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;
  
  typedef typename T::Periodic_segment                    Periodic_segment;
  typedef Periodic_2_triangulation_segment_iterator_2<T>
  Periodic_segment_iterator;
  typedef typename T::Edge                                Edge;
  typedef typename T::Edge_iterator                       Edge_iterator;
  
  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;
  
  Periodic_2_triangulation_segment_iterator_2(Iterator_type it = T::STORED)
  : _t(NULL), _it(it), _off(0) {}
  
  Periodic_2_triangulation_segment_iterator_2(const T * t,
                                              Iterator_type it = T::STORED)
  : _t(t), pos(_t->edges_begin()), _it(it), _off(0) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->edges_end() && !is_canonical() )
        ++pos;
    }
  }
  
  // used to initialize the past-the-end iterator
  Periodic_2_triangulation_segment_iterator_2(const T* t, int,
                                              Iterator_type it = T::STORED)
  : _t(t), pos(_t->edges_end()), _it(it), _off(0) {}
  
  Periodic_segment_iterator& operator++() {
    switch (_it) {
      case T::STORED:
        ++pos;
        break;
      case T::UNIQUE:
        do { ++pos; } while (pos != _t->edges_end() && !is_canonical());
        break;
      case T::STORED_COVER_DOMAIN:
      case T::UNIQUE_COVER_DOMAIN:
        increment_domain();
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return *this;
  }
  
  Periodic_segment_iterator& operator--() {
    switch (_it) {
      case T::STORED:
        --pos;
        break;
      case T::UNIQUE:
        do { --pos; } while (pos != _t->edges_begin() && !is_canonical());
        break;
      case T::STORED_COVER_DOMAIN:
      case T::UNIQUE_COVER_DOMAIN:
        decrement_domain();
    };
    return *this;
  }
  
  Periodic_segment_iterator operator++(int)
  {
    Periodic_segment_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
  
  Periodic_segment_iterator operator--(int)
  {
    Periodic_segment_iterator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Periodic_segment_iterator& ti) const
  {
    // We are only allowed to compare iterators of the same type.
    CGAL_triangulation_assertion(_it == ti._it);
    return _t == ti._t && pos == ti.pos && _off == ti._off;
  }
  
  bool operator!=(const Periodic_segment_iterator& ti) const
  {
    return !(*this == ti);
  }
  
  reference operator*() const
  {
    periodic_segment = construct_periodic_segment();
    return periodic_segment;
  }
  
  pointer operator->() const
  {
    periodic_segment = construct_periodic_segment();
    return &periodic_segment;
  }
  
  Edge_iterator get_edge() const
  {
    return pos;
  }
private:
  const T*  _t;
  Edge_iterator pos; // current edge.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_segment periodic_segment; // current segment.
  
private:
  // check whether pos points onto a unique edge or not.
  // If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    // fetch all offsets
    Offset off0, off1;
    get_edge_offsets(off0, off1);
    
    if (_t->number_of_sheets() != make_array(1,1)) {
      // If there is one offset with entries larger than 1 then we are
      // talking about a vertex that is too far away from the original
      // domain to belong to a canonical triangle.
      if (off0.x() > 1) return false;
      if (off0.y() > 1) return false;
      if (off1.x() > 1) return false;
      if (off1.y() > 1) return false;
    }
    
    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x();
    int offy = off0.y() & off1.y();
    
    return (offx == 0 && offy == 0);
  }
  
  // Artificial incrementation function that takes periodic
  // copies into account.
  void increment_domain() {
    int off = get_drawing_offsets();
    CGAL_triangulation_assertion(_off <= off);
    if (_off == off) {
      _off = 0;
      do { ++pos; } while (_it == T::UNIQUE_COVER_DOMAIN
                           && pos != _t->edges_end() && !is_canonical());
    } else {
      do {
        ++_off;
      } while ((((~_off)|off)&7)!=7); // Increment until a valid
                                      // offset has been found
    }
  }
  
  // Artificial decrementation function that takes periodic
  // copies into account.
  void decrement_domain() {
    if (_off == 0) {
      if (pos == _t->edges_begin()) return;
      do { --pos; } while (_it == T::UNIQUE_COVER_DOMAIN && !is_canonical());
      _off = get_drawing_offsets();
    } else {
      int off = get_drawing_offsets();
      do {
        --_off;
      } while ((((~_off)|off)&7)!=7); // Decrement until a valid
                                      // offset has been found
    }
  }
  
  // Get the canonicalized offsets of an edge.
  // This works in any cover that is encoded in _t->combine_offsets
  void get_edge_offsets(Offset &off0, Offset &off1) const {
    Offset cell_off0 = _t->int_to_off(pos->first->offset(_t->cw(pos->second)));
    Offset cell_off1 = _t->int_to_off(pos->first->offset(_t->ccw(pos->second)));
    Offset diff_off((cell_off0.x()==1 && cell_off1.x()==1)?-1:0,
                    (cell_off0.y()==1 && cell_off1.y()==1)?-1:0);
    off0 = _t->combine_offsets(_t->get_offset(pos->first,_t->cw(pos->second)),
                               diff_off);
    off1 = _t->combine_offsets(_t->get_offset(pos->first,_t->ccw(pos->second)),
                               diff_off);
  }
  
  // return an integer that encodes the translations which have to be
  // applied to the edge *pos
  int get_drawing_offsets() {
    Offset off0, off1;
    // Choose edges that are to be duplicated. These are edges that
    // intersect the boundary of the periodic domain. In UNIQUE mode
    // this means that the offset with respect to drawing should
    // differ in some entries. Otherwise we consider the offsets
    // internally stored inside the cell telling us that this cell
    // wraps around the domain.
    if (_it == T::UNIQUE_COVER_DOMAIN)
      get_edge_offsets(off0,off1);
    else {
      CGAL_triangulation_assertion(_it == T::STORED_COVER_DOMAIN);
      off0 = _t->int_to_off(pos->first->offset(_t->cw(pos->second)));
      off1 = _t->int_to_off(pos->first->offset(_t->ccw(pos->second)));
    }
    Offset diff_off = off0 - off1;
    
    CGAL_triangulation_assertion(diff_off.x() >= -1 || diff_off.x() <= 1);
    CGAL_triangulation_assertion(diff_off.y() >= -1 || diff_off.y() <= 1);
    
    return( 2*(diff_off.x() == 0 ? 0:1)
           + (diff_off.y() == 0 ? 0:1));
  }
  
  Periodic_segment construct_periodic_segment() const {
    CGAL_triangulation_assertion(pos->first != typename T::Face_handle());
    Offset off0, off1;
    get_edge_offsets(off0, off1);
    Offset transl_off = Offset((((_off>>1)&1)==1 ? -1:0),
                               (( _off    &1)==1 ? -1:0));
    if (_it == T::STORED_COVER_DOMAIN) {
      off0 = _t->combine_offsets(off0,transl_off);
      off1 = _t->combine_offsets(off1,transl_off);
    }
    if (_it == T::UNIQUE_COVER_DOMAIN) {
      off0 += transl_off;
      off1 += transl_off;
    }
    return make_array(
                      std::make_pair(pos->first->vertex(_t->cw(pos->second))->point(),off0),
                      std::make_pair(pos->first->vertex(_t->ccw(pos->second))->point(),off1));
  }
};

template < class T >
class Periodic_2_triangulation_point_iterator_2 {
  // Iterates over the primitives in a periodic triangulation.
  // Options:
  // - STORED: output each primitive from the Tds exactly once
  // - UNIQUE: output exactly one periodic copy of each primitive, no matter
  //     whether the current tds stores a n-sheeted covering for n!=1.
  // - STORED_COVER_DOMAIN: output each primitive whose intersection with the
  //     actually used periodic domain is non-zero.
  // - UNIQUE_COVER_DOMAIN: output each primitive whose intersection
  //     with the original domain that the user has given is non-zero
  //
  // Comments:
  // When computing in 1-sheeted covering, there will be no difference in the
  // result of STORED and UNIQUE as well as STORED_COVER_DOMAIN and
  // UNIQUE_COVER_DOMAIN.
  
public:
  typedef typename T::Periodic_point                      value_type;
  typedef const typename T::Periodic_point *              pointer;
  typedef const typename T::Periodic_point &              reference;
  typedef std::size_t                                     size_type;
  typedef std::ptrdiff_t                                  difference_type;
  typedef std::bidirectional_iterator_tag                 iterator_category;
  
  typedef typename T::Periodic_point                      Periodic_point;
  typedef Periodic_2_triangulation_point_iterator_2<T>  Periodic_point_iterator;
  
  typedef typename T::Vertex                              Vertex;
  typedef typename T::Vertex_iterator                     Vertex_iterator;
  
  typedef typename T::Offset                              Offset;
  typedef typename T::Iterator_type                       Iterator_type;
  
  Periodic_2_triangulation_point_iterator_2(Iterator_type it = T::STORED)
  : _t(NULL), _it(it) {}
  
  Periodic_2_triangulation_point_iterator_2(const T * t,
                                            Iterator_type it = T::STORED)
  : _t(t), pos(_t->vertices_begin()), _it(it) {
    if (_it == T::UNIQUE || _it == T::UNIQUE_COVER_DOMAIN) {
      while (pos != _t->vertices_end() && !is_canonical() )
        ++pos;
    }
  }
  
  // used to initialize the past-the-end iterator
  Periodic_2_triangulation_point_iterator_2(const T* t, int,
                                            Iterator_type it = T::STORED)
  : _t(t), pos(_t->vertices_end()), _it(it) {}
  
  Periodic_point_iterator& operator++() {
    switch (_it) {
      case T::STORED:
      case T::STORED_COVER_DOMAIN:
        ++pos;
        break;
      case T::UNIQUE:
      case T::UNIQUE_COVER_DOMAIN:
        do { ++pos; } while (pos != _t->vertices_end() && !is_canonical());
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return *this;
  }
  
  Periodic_point_iterator& operator--() {
    switch (_it) {
      case T::STORED:
      case T::STORED_COVER_DOMAIN:
        --pos;
        break;
      case T::UNIQUE:
      case T::UNIQUE_COVER_DOMAIN:
        do { --pos; } while (pos != _t->vertices_begin() && !is_canonical());
        break;
      default:
        CGAL_triangulation_assertion(false);
    };
    return *this;
  }
  
  Periodic_point_iterator operator++(int)
  {
    Periodic_point_iterator tmp(*this);
    ++(*this);
    return tmp;
  }
  
  Periodic_point_iterator operator--(int)
  {
    Periodic_point_iterator tmp(*this);
    --(*this);
    return tmp;
  }
  
  bool operator==(const Periodic_point_iterator& pi) const
  {
    // We are only allowed to compare iterators of the same type.
    CGAL_triangulation_assertion(_it == pi._it);
    return _t == pi._t && pos == pi.pos;
  }
  
  bool operator!=(const Periodic_point_iterator& pi) const
  {
    return !(*this == pi);
  }
  
  reference operator*() const
  {
    periodic_point = construct_periodic_point();
    return periodic_point;
  }
  
  pointer operator->() const
  {
    periodic_point = construct_periodic_point();
    return &periodic_point;
  }
  
  Vertex_iterator get_vertex() const
  {
    return pos;
  }
private:
  const T*  _t;
  Vertex_iterator pos; // current vertex.
  Iterator_type _it;
  int _off; // current offset
  mutable Periodic_point periodic_point; // current point.
  
private:
  // check whether pos points onto a vertex inside the original
  // domain. If we are computing in 1-sheeted covering this should
  // always be true.
  bool is_canonical() {
    return (_t->get_offset(pos).is_null());
  }
  
  Periodic_point construct_periodic_point() const {
    CGAL_triangulation_assertion(pos != typename T::Vertex_handle());
    Offset off = _t->get_offset(pos);
    return std::make_pair(pos->point(),off);
  }
};

template <class T>
class Domain_tester {  
  const T *t;
  
public:
  Domain_tester() {}
  Domain_tester(const T *tr) : t(tr) {}
  
  bool operator()(const typename T::Vertex_iterator & v) const {
    return (t->get_offset(v) != typename T::Offset(0,0));
  }
};

// Iterates over the vertices in a periodic triangulation that are
// located inside the original cube.
// Derives from Filter_iterator in order to add a conversion to handle
//
// Comments:
// When computing in 1-sheeted covering, there will be no difference
// between a normal Vertex_iterator and this iterator
template <class T>
class Periodic_2_triangulation_unique_vertex_iterator_2
: public Filter_iterator<typename T::Vertex_iterator, Domain_tester<T> > {
  
  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Vertex_iterator Vertex_iterator;
  
  typedef Filter_iterator<Vertex_iterator, Domain_tester<T> > Base;
  typedef Periodic_2_triangulation_unique_vertex_iterator_2 Self;
public:
  
  Periodic_2_triangulation_unique_vertex_iterator_2() : Base() {}
  Periodic_2_triangulation_unique_vertex_iterator_2(const Base &b) : Base(b) {}
  
  Self & operator++() { Base::operator++(); return *this; }
  Self & operator--() { Base::operator--(); return *this; }
  Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
  Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
  
  operator Vertex_handle() const { return Base::base(); }
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_ITERATORS_2_H
