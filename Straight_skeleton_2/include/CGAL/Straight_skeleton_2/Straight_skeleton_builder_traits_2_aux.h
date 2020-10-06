// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/certified_quotient_predicates.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Handle.h>
#include <CGAL/Quotient.h>
#include <CGAL/tags.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Unfiltered_predicate_adaptor.h>

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/optional/optional.hpp>
#include <boost/none.hpp>
#include <boost/mpl/has_xxx.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace CGAL {

namespace CGAL_SS_i {

template<class T>
T const& validate ( boost::optional<T> const& o )
{
  if ( !o )
    throw std::overflow_error("Arithmetic overflow");
  return *o ;
}

template<class NT>
NT const& validate( NT const& n )
{
  if ( !CGAL_NTS is_finite(n) )
    throw std::overflow_error("Arithmetic overflow");
  return n ;
}

// boost::make_optional is provided in Boost >= 1.34, but not before, so we define our own versions here.
template<class T>
boost::optional<T> cgal_make_optional( T const& v )
{
  return boost::optional<T>(v) ;
}

template<class T>
boost::optional<T> cgal_make_optional( bool cond, T const& v )
{
  return cond ? boost::optional<T>(v) : boost::optional<T>() ;
}

template<class K>
struct Is_filtering_kernel
{
  typedef Tag_false type ;
} ;

template<>
struct Is_filtering_kernel< Exact_predicates_inexact_constructions_kernel >
{
  typedef Tag_true type ;
} ;

//
// This is the same as Filtered_construction but uses optional<result> instead of exceptions.
//
template <class AC
         ,class EC
         ,class FC
         ,class C2E
         ,class C2F
         ,class E2C
         ,class F2C
         ,bool Protection = true
>
class Exceptionless_filtered_construction
{
private:
  EC Exact_construction;
  FC Filter_construction;
  C2E To_Exact;
  C2F To_Filtered;
  E2C From_Exact;
  F2C From_Filtered;

  typedef typename AC::result_type  AC_result_type;
  typedef typename FC::result_type  FC_result_type;
  typedef typename EC::result_type  EC_result_type;
  typedef typename C2F::Target_kernel FK;

  bool has_enough_precision(const typename FK::Point_2& point, double precision) const
  {
    return has_smaller_relative_precision(point.x(), precision) &&
           has_smaller_relative_precision(point.y(), precision);
  }

  bool has_enough_precision(const boost::tuple<typename FK::FT, typename FK::Point_2>& time_and_point, double precision) const
  {
    return has_smaller_relative_precision(boost::get<0>(time_and_point), precision) &&
           has_enough_precision(boost::get<1>(time_and_point), precision);
  }

public:

  Exceptionless_filtered_construction() {}

  Exceptionless_filtered_construction(const EC& Exact_construction, const FC& Filter_construction)
    : Exact_construction(Exact_construction)
    , Filter_construction(Filter_construction)
  {}

  typedef AC_result_type           result_type;

  template <class ... A>
  result_type
  operator()(A&& ... a) const
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      FC_result_type fr = Filter_construction(To_Filtered(std::forward<A>(a))...);

      const double precision =
        Lazy_exact_nt<double>::get_relative_precision_of_to_double();

      if ( fr && has_enough_precision(*fr, precision) )
        return From_Filtered(fr);
    }
    catch (Uncertain_conversion_exception&) {}

    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    EC_result_type er = Exact_construction(To_Exact(std::forward<A>(a))...) ;
    return From_Exact(er);
  }
};


//
// This number type is provided because unlike Quotient<> it allows you to create it
// with a zero denominator. Of course you can't evaluate it in that case, but is convenient because it allows client code
// to handle the "error" itself, which in this context is useful.
//
template<class NT>
class Rational
{
  public:

    Rational( NT aN, NT aD ) : mN(aN), mD(aD) {}

    NT n() const { return mN ; }
    NT d() const { return mD ; }

    CGAL::Quotient<NT> to_quotient() const { return CGAL::Quotient<NT>(mN,mD) ; }

    NT to_nt() const { return mN / mD ; }

    friend std::ostream& operator << ( std::ostream& os, Rational<NT> const& rat )
    {
      if ( ! CGAL_NTS is_zero(rat.d()) )
           return os << n2str(rat.n()/rat.d());
      else return os << "INF_RATIONAL" ;
    }

  private:

    NT mN, mD ;
} ;

template <class K>
struct Segment_2_with_ID
  : public Segment_2<K>
{
  typedef Segment_2<K> Base;
  typedef typename K::Point_2 Point_2;

public:
  Segment_2_with_ID() : Base(), mID(-1) { }
  Segment_2_with_ID(Base const& aS) : Base(aS), mID(-1) { }
  Segment_2_with_ID(Base const& aS, const std::size_t aID) : Base(aS), mID(aID) { }
  Segment_2_with_ID(Point_2 const& aP, Point_2 const& aQ, const std::size_t aID) : Base(aP, aQ), mID(aID) { }

public:
  std::size_t mID;
};

// @fixme this has to somehow be documented for the traits concept to make sense
//
// A straight skeleton event is the simultaneous collision of 3 offseted oriented straight line segments
// e0*,e1*,e2* [e* denotes an _offseted_ edge].
//
// This record stores the segments corresponding to the INPUT edges (e0,e1,e2) whose offsets intersect
// at the event along with their collinearity.
//
// If the event is an edge-event, then e0*->e1*->e2* must be consecutive right before the event so that
// after the event e0* and e2* become consecutive. Thus, there are _offset_ vertices (e0*,e1*) and (e1*,e2*)
// in the offset polygon which not necessarily exist in the original polygon.
//
// If the event is a split-event, e0*->e1* must be consecutive right before the event so that after the event
// e0*->right(e2*) and left(e2*)->e1* become consecutive. Thus, there is an offset vertex (e0*,e1*) in the
// offset polygon which does not necessarily exist in the original polygon.
//
// The offset vertices (e0*,e1*) and (e1*,e2*) are called the left and right seeds for the event.
// A seed is a contour node if the vertex is already present in the input polygon, otherwise is a skeleton node.
// If a seed is a skeleton node it is produced by a previous event so it is itself defined as a trisegment, thus,
// a trisegment is actually a node in a binary tree.
// Since trisegments are tree nodes they must always be handled via the nested smart pointer type: Self_ptr.
//
template<class K>
class Trisegment_2 : public Ref_counted_base
{
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2_with_ID;
  typedef Trisegment_2<K> Self;

public:
  typedef boost::intrusive_ptr<Trisegment_2> Self_ptr ;

  Trisegment_2 ( Segment_2_with_ID const&        aE0
               , Segment_2_with_ID const&        aE1
               , Segment_2_with_ID const&        aE2
               , Trisegment_collinearity aCollinearity
               , std::size_t aID
               )
    : mID(aID)
  {
    mCollinearity = aCollinearity ;

    mE[0] = aE0 ;
    mE[1] = aE1 ;
    mE[2] = aE2 ;

    switch ( mCollinearity )
    {
      case TRISEGMENT_COLLINEARITY_01:
        mCSIdx=0; mNCSIdx=2; break ;

      case TRISEGMENT_COLLINEARITY_12:
        mCSIdx=1; mNCSIdx=0; break ;

      case TRISEGMENT_COLLINEARITY_02:
        mCSIdx=0; mNCSIdx=1; break ;

      case TRISEGMENT_COLLINEARITY_ALL:
        mCSIdx = mNCSIdx = (std::numeric_limits<unsigned>::max)(); break ;

      case TRISEGMENT_COLLINEARITY_NONE:
        mCSIdx = mNCSIdx = (std::numeric_limits<unsigned>::max)(); break ;
    }
  }

  std::size_t& id() { return mID; }
  const std::size_t& id() const { return mID; }

  static Trisegment_2 null() { return Self_ptr() ; }

  Trisegment_collinearity collinearity() const { return mCollinearity ; }

  Segment_2_with_ID const& e( unsigned idx ) const { CGAL_precondition(idx<3) ; return mE[idx] ; }

  Segment_2_with_ID const& e0() const { return e(0) ; }
  Segment_2_with_ID const& e1() const { return e(1) ; }
  Segment_2_with_ID const& e2() const { return e(2) ; }

  // If 2 out of the 3 edges are collinear they can be reclassified as 1 collinear edge (any of the 2) and 1 non-collinear.
  // These methods returns the edges according to that classification.
  // PRECONDITION: Exactly 2 out of 3 edges are collinear
  Segment_2_with_ID const& collinear_edge    () const { return e(mCSIdx) ; }
  Segment_2_with_ID const& non_collinear_edge() const { return e(mNCSIdx) ; }
  Segment_2_with_ID const& other_collinear_edge() const
  {
    switch ( mCollinearity )
    {
      case TRISEGMENT_COLLINEARITY_01:
        return e(1);
      case TRISEGMENT_COLLINEARITY_12:
        return e(2);
      case TRISEGMENT_COLLINEARITY_02:
        return e(2);
    }
  }

  Self_ptr child_l() const { return mChildL ; }
  Self_ptr child_r() const { return mChildR ; }

  void set_child_l( Self_ptr const& aChild ) { mChildL = aChild ; }
  void set_child_r( Self_ptr const& aChild ) { mChildR = aChild ; }

  enum SEED_ID { LEFT, RIGHT, UNKNOWN } ;

  // Indicates which of the seeds is collinear for a normal collinearity case.
  // PRECONDITION: The collinearity is normal.
  SEED_ID degenerate_seed_id() const
  {
    Trisegment_collinearity c = collinearity();

    return c == TRISEGMENT_COLLINEARITY_01 ? LEFT : c == TRISEGMENT_COLLINEARITY_12 ? RIGHT : UNKNOWN  ;
  }

  friend std::ostream& operator << ( std::ostream& os, Self const& aTrisegment )
  {
    return os << "[" << s2str(aTrisegment.e0())
              << " " << s2str(aTrisegment.e1())
              << " " << s2str(aTrisegment.e2())
              << " " << trisegment_collinearity_to_string(aTrisegment.collinearity())
              << "]";
  }

  friend std::ostream& operator << ( std::ostream& os, Self_ptr const& aPtr )
  {
    recursive_print(os,aPtr,0);
    return os ;
  }

  static void recursive_print ( std::ostream& os, Self_ptr const& aTriPtr, int aDepth )
  {
    os << "\n" ;

    for ( int i = 0 ; i < aDepth ; ++ i )
      os << "  " ;

    if ( aTriPtr )
    {
      os << *aTriPtr ;

      if ( aTriPtr->child_l() )
        recursive_print(os,aTriPtr->child_l(),aDepth+1);

      if ( aTriPtr->child_r() )
        recursive_print(os,aTriPtr->child_r(),aDepth+1);
    }
    else
    {
      os << "{null}" ;
    }
  }

private :
  std::size_t             mID;
  Segment_2_with_ID       mE[3];
  Trisegment_collinearity mCollinearity ;
  unsigned                mCSIdx, mNCSIdx ;

  Self_ptr mChildL ;
  Self_ptr mChildR ;
} ;

template <class Info>
struct No_cache
{
  bool IsCached ( std::size_t ) { return false; }

  Info const& Get ( std::size_t )
  {
    CGAL_error();
    return Info();
  }

  void Set ( std::size_t, Info const& ) { }
  void Reset ( std::size_t ) { }
};

//TODO: call reserve, but how? #input vertices + n*m estimation?
template <class Info>
struct Info_cache
{
  std::vector<Info> mValues ;
  std::vector<bool> mAlreadyComputed ;

  bool IsCached ( std::size_t i )
  {
    return ( (mAlreadyComputed.size() > i) && mAlreadyComputed[i] ) ;
  }

  Info const& Get(std::size_t i)
  {
    CGAL_precondition ( IsCached(i) ) ;
    return mValues[i] ;
  }

  void Set ( std::size_t i, Info const& aValue)
  {
    if (mValues.size() <= i )
    {
      mValues.resize(i+1) ;
      mAlreadyComputed.resize(i+1, false) ;
    }

    mAlreadyComputed[i] = true ;
    mValues[i] = aValue ;
  }

  void Reset ( std::size_t i )
  {
    if ( IsCached(i) ) // needed if approx info is set but not exact info
      mAlreadyComputed[i] = false ;
  }
};

template <typename K>
using Time_cache = Info_cache< boost::optional< CGAL_SS_i::Rational< typename K::FT > > > ;

template <typename K>
using Coeff_cache = Info_cache< boost::optional< Line_2<K> > > ;

template<class K>
struct Functor_base_2
{
  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  typedef typename K::Vector_2 Vector_2 ;
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2_with_ID ;

  typedef CGAL_SS_i::Trisegment_2<K> Trisegment_2 ;

  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;
};

template<class Converter>
struct SS_converter : Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Source_kernel::FT Source_FT ;
  typedef typename Target_kernel::FT Target_FT ;

  typedef typename Source_kernel::Point_2 Source_point_2 ;
  typedef typename Target_kernel::Point_2 Target_point_2 ;

  typedef typename Source_kernel::Vector_2 Source_vector_2 ;
  typedef typename Target_kernel::Vector_2 Target_vector_2 ;

  typedef typename Source_kernel::Segment_2 Source_segment_2 ;
  typedef typename Target_kernel::Segment_2 Target_segment_2 ;

  typedef Segment_2_with_ID<Source_kernel> Source_segment_2_with_ID ;
  typedef Segment_2_with_ID<Target_kernel> Target_segment_2_with_ID ;

  typedef Trisegment_2<Source_kernel> Source_trisegment_2 ;
  typedef Trisegment_2<Target_kernel> Target_trisegment_2 ;

  typedef boost::tuple<Source_FT,Source_point_2> Source_time_and_point_2 ;
  typedef boost::tuple<Target_FT,Target_point_2> Target_time_and_point_2 ;

  typedef boost::optional<Source_FT> Source_opt_FT ;
  typedef boost::optional<Target_FT> Target_opt_FT ;

  typedef boost::optional<Source_point_2> Source_opt_point_2 ;
  typedef boost::optional<Target_point_2> Target_opt_point_2 ;

  typedef boost::optional<Source_time_and_point_2> Source_opt_time_and_point_2 ;
  typedef boost::optional<Target_time_and_point_2> Target_opt_time_and_point_2 ;

  typedef boost::optional<Source_segment_2> Source_opt_segment_2 ;
  typedef boost::optional<Target_segment_2> Target_opt_segment_2 ;

  typedef typename Source_trisegment_2::Self_ptr Source_trisegment_2_ptr ;
  typedef typename Target_trisegment_2::Self_ptr Target_trisegment_2_ptr ;


  Target_FT cvt_n(Source_FT const& n) const  { return this->Converter::operator()(n); }

  Target_opt_FT cvt_n(Source_opt_FT const& n) const
  {
    Target_opt_FT r ;
    if ( n )
      r = cvt_n(*n);
    return r ;
  }

  Target_point_2   cvt_p(Source_point_2 const& p) const  { return this->Converter::operator()(p); }

  Target_vector_2 cvt_v( Source_vector_2 const& v) const {
    return Target_vector_2(cvt_p(Source_point_2(CGAL::ORIGIN)), cvt_p(Source_point_2(CGAL::ORIGIN) + v) ) ;
  }

  Target_segment_2 cvt_s( Source_segment_2 const& e) const {
    return Target_segment_2(cvt_p(e.source()), cvt_p(e.target())) ;
  }

  Target_segment_2_with_ID cvt_s( Source_segment_2_with_ID const& e) const {
    return Target_segment_2_with_ID(cvt_p(e.source()), cvt_p(e.target()), e.mID) ;
  }

  Target_time_and_point_2 cvt_t_p( Source_time_and_point_2 const& v ) const
  {
    Source_FT      t ;
    Source_point_2 p ;
    boost::tie(t,p) = v ;
    return Target_time_and_point_2(cvt_n(t),cvt_p(p));
  }

  Target_trisegment_2_ptr cvt_single_trisegment( Source_trisegment_2_ptr const& tri ) const
  {
    CGAL_precondition( tri!= Source_trisegment_2_ptr() ) ;

    return Target_trisegment_2_ptr ( new Target_trisegment_2(cvt_s(tri->e0())
                                                            ,cvt_s(tri->e1())
                                                            ,cvt_s(tri->e2())
                                                            ,tri->collinearity()
                                                            ,tri->id()
                                                            )
                                   ) ;
  }

  Target_trisegment_2_ptr cvt_trisegment( Source_trisegment_2_ptr const& tri ) const
  {
    Target_trisegment_2_ptr res ;

    if ( tri )
    {
      res = cvt_single_trisegment(tri) ;

      if ( tri->child_l() )
        res->set_child_l( cvt_trisegment(tri->child_l()) ) ;

      if ( tri->child_r() )
        res->set_child_r( cvt_trisegment(tri->child_r() ) ) ;
    }

    return res ;
  }

  bool operator()( bool v ) const { return v ; }

  Trisegment_collinearity  operator()(Trisegment_collinearity c) const { return c ; }

  Oriented_side operator()(Oriented_side s) const { return s ; }

  Target_FT        operator()(Source_FT const& n) const { return cvt_n(n) ; }

  Target_opt_FT    operator()(Source_opt_FT const& n) const { return cvt_n(n) ; }

  Target_point_2   operator()( Source_point_2 const& p) const { return cvt_p(p) ; }

  Target_vector_2 operator()( Source_vector_2 const& v) const { return cvt_v(v); }

  Target_segment_2 operator()( Source_segment_2 const& s) const { return cvt_s(s); }

  Target_segment_2_with_ID operator()( Source_segment_2_with_ID const& s) const { return cvt_s(s); }

  Target_trisegment_2_ptr operator()( Source_trisegment_2_ptr const& tri ) const
  {
    return cvt_trisegment(tri);
  }

  Target_time_and_point_2 operator() ( Source_time_and_point_2 const& v ) const
  {
    return cvt_t_p(v);
  }

  Target_opt_point_2 operator()( Source_opt_point_2 const& p) const
  {
    if ( p )
         return Target_opt_point_2(cvt_p(*p));
    else return Target_opt_point_2();
  }

  Target_opt_segment_2 operator()( Source_opt_segment_2 const& s) const
  {
    if ( s )
         return Target_opt_segment_2(cvt_s(*s));
    else return Target_opt_segment_2();
  }

  Target_opt_time_and_point_2 operator()( Source_opt_time_and_point_2 const& v) const
  {
    if ( v )
         return Target_opt_time_and_point_2(cvt_t_p(*v));
    else return Target_opt_time_and_point_2();
  }

};

BOOST_MPL_HAS_XXX_TRAIT_DEF(Filters_split_events_tag)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Segment_2_with_ID)

} // namespace CGAL_SS_i


//
// This macro defines a global functor adapter which allows users to use it in the followig ways:
//
// Given a 'Functor' provided by a given 'Traits' (or Kernel):
//
//   typedef typename CGAL::Functor<Traits>::type Functor ;
//   result r = CGAL::Functor<Traits>(traits)(a,b,c);
//
#define CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(functor) \
        template<class K> \
        typename K :: functor functor ( K const& aK ) \
        { \
          return aK.get((typename K :: functor const*)0);  \
        }


} // end namespace CGAL


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H //

// EOF //
