//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef __VC_MATH_VECN_HH__
#define __VC_MATH_VECN_HH__

#include "array.hh"

#include <cmath>
#include <cassert>
#include <cstdio>

#include <type_traits>
#include <limits>
#include <typeinfo>
#include <exception>
#include <iostream>

//#include "../system/platform.hh" // isnan, isinf, isfinite

# if defined(_MSC_VER)
#  pragma warning( disable : 4290 ) // throw(/*type*/) is parsed but handled as throw(/*...*/)
#  define isnan(x) _isnan(x)
#  define isinf(x) (!_finite(x))
#endif

namespace VC {
namespace math {

/** \file VecN.hh
    math vectors, [\ref vc_math_lav]
*/

/** @name Utilities for initialization by VecN::operator=()
    \ingroup vc_math_lav
    @{
*/

/// \enum VC::math::LIST_T VecN::operator=()
enum LIST_T {
  LIST //!< begin intialization by comma seperated list (VecN::operator=())
};
/// \enum VC::math::FILL_T VecN::operator=()
enum FILL_T {
  FILL //!< fill until end  with the last data value provided (see VecN::operator=())
};
/// \enum VC::math::END_T VecN::operator=()
enum END_T {
  ENDL //!< end list initialization and return array type (see VecN::operator=())
};

# ifndef DOXYGEN_SKIP

template <typename T,unsigned L> class _ScalarListRep {
public:
  _ScalarListRep(const T& _value) : value(_value) {}
  T value;
};

/// LIST initialization: repeat `_value` `L` times, see VecN::operator=()
template <typename T,unsigned L> _ScalarListRep<T,L>
REP(const T& _value) { return _ScalarListRep<T,L>(_value); }
template <unsigned L> _ScalarListRep<double,L>
REP(const double& _value) { return _ScalarListRep<double,L>(_value); }

/// LIST initialization used by VecN::operator=()
template <typename Ary,typename T,unsigned K,int L>
class ScalarList {
public:
  ScalarList(T* _p) : m_p(_p) {}

  ScalarList<Ary,T,K,L-1>& operator,(const T& _t) {
    static_assert(L>0,"too many elements in list initialization");
    //printf("m_p[%d]=%f\n",K-L,_t);
    m_p[K-L]=_t;
    return *((ScalarList<Ary,T,K,L-1>*) this);
  }
  ScalarList<Ary,T,K,0>& operator,(FILL_T) {
    static_assert(K-L>0,"fill requires preceding element in list initialization");
    return (*this).operator,(_ScalarListRep<T,L>(m_p[K-L-1]));
  }
  template <typename S,unsigned N>
  ScalarList<Ary,T,K,L-N>& operator,(const _ScalarListRep<S,N>& _r) {
    static_assert(L>0,"too many elements in list initialization");
    //printf("m_p[%d]=%f\n",K-L,T(_r.value));
    m_p[K-L]=T(_r.value);
    return *((ScalarList<Ary,T,K,L-1>*) this),_ScalarListRep<T,N-1>(_r.value);
  }
  template <typename S>
  ScalarList<Ary,T,K,L>& operator,(const _ScalarListRep<S,0>&) {
    return *this;
  }
  Ary& operator,(END_T) {
    static_assert(K-L>0,"fill requires preceding element in list initialization");
    return *((Ary*) m_p);
  }

private:
  T* m_p;
};

# else

/// LIST initialization: repeat `_value` `L` times, see VecN::operator=()
template <unsigned N> some_type REP(const T& _value);

# endif // DOXYGEN_SKIP

/// @}


# ifndef DOXYGEN_SKIP
template <typename T,unsigned K>
struct VecN_align {};

# if defined(__GNUC__) && defined(__SSE__) && !defined(__CUDACC__)

//
// How about aliasing?
//

template<>
struct VecN_align<float,2> {
  float __attribute__ ((__vector_size__ (8))) m_v2f;
};
template<>
struct VecN_align<float,4> {
  float __attribute__ ((__vector_size__ (16))) m_v4f;
};
template<>
struct VecN_align<double,2> {
  double __attribute__ ((__vector_size__ (16))) m_v2d;
};

# endif

# if defined(__GNUC__) && defined(__AVX__) && !defined(__CUDACC__)

template<>
struct VecN_align<float,8> {
  float __attribute__ ((__vector_size__ (32))) m_v8f;
};
template<>
struct VecN_align<double,4> {
  double __attribute__ ((__vector_size__ (32))) m_v4d;
};

# endif

//
// Can we define a "kernel" to have several operations work on these data types?
// -> element wise +,-,*,/ (Array.hh)
//



# endif // DOXYGEN_SKIP

/** \def VC_VECN_CHECK_BOUNDS
    \ingroup vc_math_lav
    Turn on bounds checking (using \c assert) in [\ref vc_math_lav].
*/

# ifdef DOXYGEN_SKIP
#  define VC_VECN_CHECK_BOUNDS "this is an external switch"
#  error "doxygen only"
# endif


/** \class VC::math::VecN VecN.hh
    \brief N-Vectors.

    \ingroup vc_math_lav

    **A vector is simply an array:** The VecN<T,K> vector data
    structure is guaranteed to have the same layout as a
    K-dimensional C-array of type T.

    VecN defines operations on these arrays, many of them conforming
    with the algebraic operations defined on a K-dimensional
    Eucledian vector space. Additional operations include
    construction, permutation, sorting, and other higher-level
    operations. (Some operations are defined externally as \c friend
    functions.)

    Operations are based on respective functions on arrays defined
    in Array.hh (see [\ref vc_math_lav], FixedAryOps).

    **Summary on operators**

    \code
    VecN<double,2> a,b,c;
    VecN<float,2>  f;

    // -- elementwise operations --

    // addition/subtraction
    a+=b;    a-=b;   // vector addition
    a+=1.0;  a-=1.0; // add/subtract scalar to all components

    c=a+b;   a=a-b;  // vector addition
    c=a+1.0; c=a-1;  // add scalar to all components

    a=-a;            // unary minus, negation

    // multiplication and division
    a*=b;    a/=b;   // This is not a scalar product or matrix inverse!
    a*=2.0;  a/=2.0; // scale vector

    c=a*b;   c=a/b;
    c=a*2.0; c=a/2.0;

    a.invert();      // elementwise 1.0/a[i]
    a.invert(2.0);   // elementwise 2.0/a[i]

    // -- vector/vector operations

    double d=(a|b);  // scalar product
    c=(a%b);         // cross product, defined only for N=3

    \endcode

    Remarks

    \arg Operators +,-,*,/ all work \a elementwise! Especially,
    VecN*=VecN and VecN/=Vecn (operator*=(),operator/=()) are \a
    elementwise multiplication and division, respectively.
    \code
    VecN<double,2> a(1,2),b(3,4),c;
    a+=b;  // yields (4,6)
    a*=b;  // yields (3,8)
    c=a/b; // yields (1/3,2/4)
    \endcode

    \arg The scalar product is defined as dot() or operator|().  For
    the latter, note that C++ operator preceedence does \a not conform
    with <|>: we recommend to write "(a|b)" rather than "a|b".
    \code
    VecN<double,2> a(1,2),b(3,4);
    double d=(a|b); // yields 1*3+2*4
    \endcode

    \arg Addition/substraction of \a scalars VecN+T adds/subtracts
    the scalar to/from every vector \a element (see, e.g.,
    operator+=()).
    \code
    VecN<double,2> a(1,2);
    a+=1;  // yields (2,3)
    a*=2;  // yields (2,4)
    \endcode

    \arg Binary operators are defined externally and declared
    `friend`, e.g., `operator+()`, `operator*()`.

    \arg For binary operations, the scalar is always the right (second)
    operand.

    \arg The use of assignment operators, e.g., `a+=b`
    (`operator+=())`, is generally more \a efficient than using binary
    operators `a+b`.

    \arg Consider the use of lincomb() for computing linear
    combinations.
    \code
    VecN<double,2> a(1,2),b(3,4), c;
    c.lincomb(0.25,a,0.75,b); // yields c=a*0.25+b*0.75
    \endcode

    \arg `operator=()` enables initialization by a comma separated list
    of values.
    \code
    VecN<double,4> x;
    x=LIST,1,2,3,4;    // see operator=() below for details
    \endcode

    \arg to_v() and to_cv() cast any arrays to VecN vectors.
    (External version are VC::math::to_v() and
    VC::math::to_cv()).
    \code VecN<double,2> a(1,2), b(3,4);
    double b[]={3,4};

    a+=Vector<double,2>::to_cv(b); // yields a=(4,6)
    Vector<double,2>::to_v(b)=c;   // yields c[0]=3, c[1]=4

    a+=to_cv<2>(b);                // same externally
    to_v<2>(b)=c;                  // same externally
    \endcode

    \arg In many binary operations Result is used to determine the
    resulting scalar type, i.e., `VecN<float,N>+VecN<double,N>`
    yields a `VecN<double,N>`. -- Mind automatic type conversion!
    \code
    VecN<float,2>  a;
    VecN<double,2> b;
    a+b; // yields Vector<double,2> as result
    \endcode

    \arg Use map() to apply a function `f(a)` (or `f(a,b)`) element-wise
    to a (and b, respectively).
    \code
    VecN<double,3>  a(1,2,3),x;
    x.map(a,sqrt);              // load elementwise sqrt(a[i]) to x
    x*=x;                       // yields a again, ||x-a||=eps
    \endcode

    \arg Several functions are defined only for floating point data
    types (e.g., normalize()) or certain dimensions (e.g.,
    the cross product `operator%=()`).

    \arg Method swap() and VC::math::swap() exchange vector contents
    element-wise. The latter is an alternative to std::swap().

    \arg If VC_VECN_CHECK_BOUNDS is defined then operator[]()
    and function nth() include an assertion on index within bounds.

    \sa [\ref vc_math_lav_vops], [\ref vc_math_lav]
*/

template <typename T,unsigned K>
class VecN {
public:
  typedef T value_type;       //!< scalar type
  typedef FixedAryOps<K> Ops;

  typedef VecN<T,K> Self;     //!< synonym for `VecN<T,K>`

  enum Dimensions {
    ROWS=K,  //!< number of rows (=K=n_rows())
    COLS=1,  //!< number of columns (=1=n_cols())
    N=K,     //!< N=K, synonym for ROWS
    SIZE=K   //!< number of elements (=K=N=size())
  };

  /** @name initialization
      Specializations for `Vec<T,2>`, `VecN<T,3>`, `VecN<T,4>` take 2,3,4 scalar
      arguments for initialization.
      @{
  */

  VecN() {} //!< uninitialized vector

  VecN(const Self& _v) { memcpy(m_v,_v.m_v,sizeof(*this)); } //!< copy constructor

  /// elementwise copy with type cast
  template <typename S>
  VecN(const VecN<S,K>& _v) {
    for (unsigned i=0;i<K;++i) m_v[i]=T(_v.data()[i]);
  }
  /// set *all* elements to `_s`
  explicit VecN(const T& _s) {
    for (unsigned i=0;i<K;++i) m_v[i]=_s;
  }
  /// construct from array (elementwise copy)
  explicit VecN(const T* _v) {
    for (unsigned i=0;i<K;++i) m_v[i]=_v[i];
  }
  /// construct from array (elementwise copy)
  template <typename S>
  explicit VecN(const S* _v) {
    for (unsigned i=0;i<K;++i) m_v[i]=T(_v[i]);
  }

  /// construct 2-vector
  VecN(const T& _v0,const T& _v1) {
    static_assert(K==2,"dimension mismatch");
    m_v[0]=_v0; m_v[1]=_v1;
  }
  /// construct 3-vector
  VecN(const T& _v0,const T& _v1,const T& _v2) {
    static_assert(K==3,"dimension mismatch");
    m_v[0]=_v0; m_v[1]=_v1; m_v[2]=_v2;
  }
  /// construct 4-vector
  VecN(const T& _v0,const T& _v1,const T& _v2,const T& _v3) {
    static_assert(K==4,"dimension mismatch");
    m_v[0]=_v0; m_v[1]=_v1; m_v[2]=_v2; m_v[3]=_v3;
  }
  /// construct 5-vector
  VecN(const T& _v0,const T& _v1,const T& _v2,const T& _v3,const T& _v4) {
    static_assert(K==5,"dimension mismatch");
    m_v[0]=_v0; m_v[1]=_v1; m_v[2]=_v2; m_v[3]=_v3; m_v[4]=_v4;
  }
  /// construct 6-vector
  VecN(const T& _v0,const T& _v1,const T& _v2,const T& _v3,const T& _v4,const T& _v5) {
    static_assert(K==6,"dimension mismatch");
    m_v[0]=_v0; m_v[1]=_v1; m_v[2]=_v2; m_v[3]=_v3; m_v[4]=_v4; m_v[5]=_v5;
  }

  ~VecN() {} //!< (empty)


  /// `(*this)=_v`
  Self& operator=(const Self& _v) {
    memcpy(m_v,_v.m_v,sizeof(*this)); return *this;
  }
  /// `(*this)=_v`
  template <typename S>
  Self& operator=(const VecN<S,K>& _v) {
    for (unsigned i=0;i<K;++i) m_v[i]=_v.data()[i];
    return *this;
  }

  /** Initialize with comma separated list of scalars.
      \arg Use LIST, REP, and LIST to repeat values or to fill up to end() of
      vector. Bounds checking is done at compile time.

      \code
      VecN<double,3> x;

      // initialize by comma separated list
      // equivalent to x=VecN<double,3>(1.0,2.0,3.0);
      x=LIST,1.0,2.0,3.0;

      // use Rep<K> to repeat scalar K times in initialization
      // equivalent to  x=VecN<double,3>(2.0,2.0,1.0)
      x=LIST,REP<2>(2.0),1.0;

      // use Fill to fill remainder of vector with scalar
      // equivalent to x=VecN<double,3>(1.0,1.0,1.0)
      x=LIST,1.0,FILL;

      // using will return the modified array
      (x=LIST,1.0,2.0,3.0,ENDL)+=1;  // yields 2.0,3.0,4.0
      // Mind operator preceedence!

      \endcode
      LIST, REP, FILL, ENDL (ScalarList used internally)

      Caveat
      \arg Mind operator precedence!
      \code
      VecN<double,3> x;
      x=(LIST,1.0,2.0,3.0);     // fails at runtime (!)
      //  because list is collapsed before assignment
      \endcode
      \arg List initialization is available only for the assignment operator,
      \a not for the constructor!
      \code
      VecN<double,3> x=LIST,1.0,2.0,3.0; // fails at compile time
      // because constructor is called!
      \endcode
      \arg The purpose of LIST assignment is use within a single
      assignment statement, e.g., for initialization.
  */

  ScalarList<Self,T,K,K> operator=(LIST_T) { return ScalarList<Self,T,K,K>(m_v); }

  /// concatenate vectors "vertically", returns `[_v0;_v1]`
  template <typename T0,unsigned K0,typename T1,unsigned K1>
  friend VecN<T0,K0+K1> v_cat(const VecN<T0,K0>& _v0,const VecN<T1,K1>& _v1);

  /// concatenate vectors "vertically", returns `[_v0;_v1;_v2]`
  template <typename T0,unsigned K0,typename T1,unsigned K1,typename T2,unsigned K2>
  friend VecN<T0,K0+K1+K2> v_cat(const VecN<T0,K0>& _v0,
                                 const VecN<T1,K1>& _v1,
                                 const VecN<T2,K2>& _v2);

  /// concatenate vectors "vertically", assigns `*this=[_v0;_v1]`
  template <typename T0,unsigned K0,typename T1,unsigned K1>
  Self& v_cat(const VecN<T0,K0>& _v0,const VecN<T1,K1>& _v1) {
    static_assert(K0+K1==K,"dimension mismatch");
    unsigned i;
    for (i=0;i<K0   ;++i)  m_v[i]=T(_v0.data()[i]);
    for (   ;i<K0+K1;++i)  m_v[i]=T(_v1.data()[i-K0]);
    return *this;
  }
  /// concatenate vectors "vertically", assigns `*this=[_v0;_v1;_v2]`
  template <typename T0,unsigned K0,typename T1,unsigned K1,typename T2,unsigned K2>
  Self& v_cat(const VecN<T0,K0>& _v0,const VecN<T1,K1>& _v1,const VecN<T2,K2>& _v2) {
    static_assert(K0+K1+K2==K,"dimension mismatch");
    unsigned i;
    for (i=0;i<K0      ;++i) m_v[i]=T(_v0.data()[i]);
    for (   ;i<K0+K1   ;++i) m_v[i]=T(_v1.data()[i-K0]);
    for (   ;i<K0+K1+K2;++i) m_v[i]=T(_v2.data()[i-K0-K1]);
    return *this;
  }

  /// @}

  /** @name dimensions
      @{
  */

  unsigned n_rows() const { return K; } //!< number of rows (`=K=N`)
  unsigned n_cols() const { return 1; } //!< number of columns (`=1`)
  unsigned size() const { return K; }   //!< number of elements (`=K=N`)

  /// @}

  /** @name access and assignment
      If \c VC_VECN_CHECK_BOUNDS is defined then
      `VecN::operator[]()` and function nth() include an assertion on index
      within bounds.
      @{
  */

  T* data() { return m_v; }              //!< get raw data
  const T* data() const { return m_v; }  //!< get raw data

  T* begin() { return m_v; }              //!< synonym for data()
  const T* begin() const { return m_v; }  //!< synonym for data()

  T* end() { return m_v+K; }              //!< synonym for data()+`N`
  const T* end() const { return m_v+K; }  //!< synonym for data()+`N`


# ifdef VC_VECN_CHECK_BOUNDS
  /// `distance(_p,begin())` assumed within `[0,K-1]` (assertion)
  int nth(const T* _p) {
    int i=_p-&m_v[0];
    assert(0<=i); assert(i<int(K));
    return i;
  }
  /// access element (assertion)
  T& operator[](unsigned _i) {
    assert(_i<K); return m_v[_i];
  }
  /// access element (assertion)
  const T& operator[](unsigned _i) const {
    assert(_i<K); return m_v[_i];
  }
  /// access element (assertion)
  T& operator[](int _i) {
    assert(_i>=0); assert(_i<int(K)); return m_v[_i];
  }
  /// access element (assertion)
  const T& operator[](int _i) const {
    assert(_i>=0); assert(_i<int(K)); return m_v[_i];
  }
# else
  /// distance(_p,begin()) assumed within [0,K-1]
  int nth(const T* _p) { return _p-&m_v[0]; }
  T& operator[](unsigned _i) { return m_v[_i]; } //!< access element
  const T& operator[](unsigned _i) const { return m_v[_i]; } //!< access element
#endif

  /// cast C-array (also defined externally as function)
  static Self&
  to_v(T* _p) { return (Self&) *_p; }

  /// cast const C-array (also defined externally as function)
  static const Self&
  to_cv(const T* _p) { return *((const Self*) _p); }


  /// swap contents of (*this) and _v
  template <typename S>
  inline Self& swap(VecN<S,K>& _v)  {
    Ops::swap(m_v, _v.m_v);
    return *this;
  }

  /// alternative to std::swap(_a,_b)
  template <typename A,unsigned N>
  friend void
  swap(VecN<A,N>& _a,VecN<A,N>& _b);

  /// external definition of `VecN<T,K>::to_v()`
  template <unsigned J,typename U>
  VecN<U,J>& to_v(U* _v); // for documentation only
  /// external definition of `VecN<T,K>::to_cv()`
  template <unsigned J,typename U>
  const VecN<U,J>& to_cv(const U* _v); // for documentation only

  /// @}


  /** @name arithmetic: elementwise operations
      See also \ref vc_math_lav_vops
      @{
  */

  //
  // elementwise OP=SCALAR
  //

  /// \a element-wise addition of scalar `(*this)+=_s (adds())`
  template <typename S>
  inline Self& operator+=(const S& _s)  {
    Ops::adds(m_v, m_v,_s);
    return *this;
  }
  /// \a element-wise subtraction of scalar `(*this)-=_s (subs())`
  template <typename S>
  inline Self& operator-=(const S& _s)  {
    Ops::subs(m_v, m_v,_s);
    return *this;
  }
  /// multiplication with scalar `(*this)*=_s (muls())`
  template <typename S>
  inline Self& operator*=(const S& _s)  {
    Ops::muls(m_v, m_v,_s);
    return *this;
  }
  /// division by scalar `(*this)/=_s (divs())`
  template <typename S>
  inline Self& operator/=(const S& _s)  {
    Ops::divs(m_v, m_v,_s);
    return *this;
  }

  /// \a element-wise `(*this)= -(*this) (subs2())`
  inline Self& negate() {
    Ops::muls(m_v, m_v,T(-1));
    return *this;
  }
  /// \a element-wise division `(*this)=_s/(*this) (divs2())`
  template <typename S>
  inline Self& invert(const S& _s)  {
    Ops::divs2(m_v, _s,m_v);
    return *this;
  }
  /// \a element-wise `(*this)=1/(*this) (divs2())`
  inline Self& invert() {
    Ops::divs2(m_v, T(1),m_v);
    return *this;
  }

  ///

  //
  // element=wise OP=VecN
  //

  /// `(*this)+=_v (addv())`
  template <typename S>
  inline Self& operator+=(const VecN<S,K>& _v)  {
    Ops::addv(m_v, m_v,_v.m_v);
    return *this;
  }
  /// `(*this)-=_v (subv())`
  template <typename S>
  inline Self& operator-=(const VecN<S,K>& _v)  {
    Ops::subv(m_v, m_v,_v.m_v);
    return *this;
  }
  /// \a element-wise multiplication `(*this)*=_v (mulv())`
  template <typename S>
  inline Self& operator*=(const VecN<S,K>& _v)  {
    Ops::mulv(m_v, m_v,_v.m_v);
    return *this;
  }
  /// \a element-wise division `(*this)-=_v (divv())`
  template <typename S>
  inline Self& operator/=(const VecN<S,K>& _v)  {
    Ops::divv(m_v, m_v,_v.m_v);
    return *this;
  }

  /// @}

  /** @name arithmetic: operations on this VecN
      @{
  */

  /// sum of elements
  T sum() const { return Ops::sum(m_v); }

  /// arithmetic mean of elements (only `float`, `double`)
  T mean() const { return Ops::mean(m_v); }

  /// square `<(*this)|(*this)>`
  T sqr() const { return Ops::sqr(m_v); }

  /// Euclidean norm `sqrt(sqr())`
  T norm() const {
    static_assert(std::is_floating_point<T>::value,"require floating point type");
    return Ops::norm(m_v);
  }

  /// get vector of absolute values
  Self abs() const {
    Self v(*this);
    for (int i=0;i<K;++i)
      v[i]=fabs(v[i]);
    return v;
  }

  /// 1-norm `sum(abs())`
  T norm1() const { return Ops::norm1(m_v); }

  /// 2-norm, synonym for norm()
  T norm2() const { return Ops::norm(m_v); }

  /// maximum norm
  T norminf() const { return Ops::norminf(m_v); }

  /// Euclidean norm `sqrt(sqr(_v))`
  template <typename S,unsigned N>
  friend S norm(const VecN<S,N>& _v);

  /// 1-norm `sum(abs(_v))`
  template <typename S,unsigned N>
  friend S norm1(const VecN<S,N>& _v);

  /// 2-norm, synonym for norm()
  template <typename S,unsigned N>
  friend S norm2(const VecN<S,N>& _v);

  /// maximum norm VecN::norminf()
  template <typename S,unsigned N>
  friend S norminf(const VecN<S,N>& _v);

  /// synomym for `isfinite(sum())`
  bool is_finite() const {
    static_assert(std::is_floating_point<T>::value,"require floating point type");
    return std::isfinite(Ops::sum(m_v));
  }

  /// unary minus: return -_v (negate())
  template <typename S,unsigned N>
  friend  VecN<S,N>
  operator-(const VecN<S,N>& _v);

  /** Normalize `this` to unit vector.
      Normalizing a null vector leads to `!is_finite(this))`.
  */
  Self& normalize() {
    static_assert(std::is_floating_point<T>::value,"require floating point type");
    T s=T(1)/Ops::norm(m_v);
    Ops::muls(m_v,m_v,s);
    return *this;
  }

  /// @}

  /** @name arithmetic: vector-vector operations
      See also \ref vc_math_lav_vops
      @{
  */

  /// dot product `<(*this)|_v> (dot())`
  template <typename S>
  typename Result<T,S>::type operator|(const VecN<S,K>& _v) const {
    return Ops::dot(m_v,_v.m_v);
  }

  /// load linear combination `(*this)=_a*_sa+_b*_sb (FixedAryOps::lincomb())`
  template <typename A,typename B>
  Self& lincomb(const A& _sa,const VecN<A,N>& _a,
                const B& _sb,const VecN<B,N>& _b) {
    FixedAryOps<N>::lincomb(m_v, _sa,_a.m_v, _sb,_b.m_v);
    return *this;
  }
  /// load linear combination `(*this)=_v1*_s1+_v2*_s2+_v3*_s3 (FixedAryOps::lincomb())`
  template <typename T1,typename T2,typename T3>
  Self& lincomb(const T1& _s1,const VecN<T1,N>& _v1,
                const T2& _s2,const VecN<T2,N>& _v2,
                const T3& _s3,const VecN<T3,N>& _v3) {
    FixedAryOps<N>::lincomb(m_v, _s1,_v1.m_v, _s2,_v2.m_v, _s3,_v3.m_v);
    return *this;
  }
  /// load linear combination `(*this)=_v1*_s1+_v2*_s2+_v3*_s3+_v4*_s4 (FixedAryOps::lincomb())`
  template <typename T1,typename T2,typename T3,typename T4>
  Self& lincomb(const T1& _s1,const VecN<T1,N>& _v1,
                const T2& _s2,const VecN<T2,N>& _v2,
                const T3& _s3,const VecN<T3,N>& _v3,
                const T4& _s4,const VecN<T4,N>& _v4) {
    FixedAryOps<N>::lincomb(m_v,
                            _s1,_v1.m_v, _s2,_v2.m_v, _s3,_v3.m_v, _s4,_v4.m_v);
    return *this;
  }

  /// load linear combination `(*this)=_v0*_w[0]+_v1*_w[1] (FixedAryOps::lincomb())`
  template <typename W,typename T0,typename T1>
  Self& lincomb(const VecN<W,2>& _w,const VecN<T0,N>& _v0,const VecN<T1,N>& _v1) {
    FixedAryOps<N>::lincomb(m_v, _w.data()[0],_v0.m_v, _w.data()[1],_v1.m_v);
    return *this;
  }
  /// load linear combination `(*this)=_v0*_w[0]+_v1*_w[1]+_v2*_w[2]` (FixedAryOps::lincomb())
  template <typename W,typename T0,typename T1,typename T2>
  Self& lincomb(const VecN<W,3>& _w,
                const VecN<T0,N>& _v0,const VecN<T1,N>& _v1,const VecN<T2,N>& _v2) {
    FixedAryOps<N>::lincomb(m_v,
                            _w.data()[0],_v0.m_v,_w.data()[1],_v1.m_v,
                            _w.data()[2],_v2.m_v);
    return *this;
  }
  /// load linear combination `(*this)= { sum _i=0,...,3 _vi*_w[i] }` (FixedAryOps::lincomb())
  template <typename W,typename T0,typename T1,typename T2,typename T3>
  Self& lincomb(const VecN<W,4>& _w,
                const VecN<T0,N>& _v0,const VecN<T1,N>& _v1,
                const VecN<T2,N>& _v2,const VecN<T3,N>& _v3) {
    FixedAryOps<N>::lincomb(m_v,
                            _w.data()[0],_v0.m_v,_w.data()[1],_v1.m_v,
                            _w.data()[2],_v2.m_v,_w.data()[3],_v3.m_v);
    return *this;
  }


  /// load cross product `*this = *this x _v` (`N=3`) \ingroup vc_math_lav_vops
  template <typename S>
  Self& operator%=(const VecN<S,3>& _v) {
    static_assert(K==3,"dimension mismatch");
    T x=m_v[1]*_v.m_v[2]-m_v[2]*_v.m_v[1];
    T y=m_v[2]*_v.m_v[0]-m_v[0]*_v.m_v[2];
    T z=m_v[0]*_v.m_v[1]-m_v[1]*_v.m_v[0];
    m_v[0]=x; m_v[1]=y; m_v[2]=z;
    return *this;
  }

  ///  vector addition `_a+_b` (addv()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type,N>
  operator+(const VecN<A,N>& _a,const VecN<B,N>& _b);

  /// vector subtraction `_a+_b` (addv()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type ,N>
  operator-(const VecN<A,N>& _a,const VecN<B,N>& _b);

  /// \a element-wise multiplication (mulv()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type ,N>
  operator*(const VecN<A,N>& _a,const VecN<B,N>& _b);

  /// \a element=wise division (mulv()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type ,N>
  operator/(const VecN<A,N>& _a,const VecN<B,N>& _b);

  /// add scalar `_b` to all components of `_a` (adds()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type,N>
  operator+(const VecN<A,N>& _a,const B& _b);

  /// subtract scalar `_b` from all components of `_a` (subs()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type,N>
  operator-(const VecN<A,N>& _a,const B& _b);

  /// multiply scalar `_b` to all components of `_a` (muls()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type,N>
  operator*(const VecN<A,N>& _a,const B& _b);

  /// divide components of `_a` by `_b` (divs()) \ingroup vc_math_lav_vops
  template <typename A,typename B,unsigned N>
  friend  VecN<typename Result<A,B>::type,N>
  operator/(const VecN<A,N>& _a,const B& _b);

  /// cross productof 3-vectors \ingroup vc_math_lav_vops
  template <typename A,typename B>
  friend VecN<typename Result<A,B>::type,3>
  operator%(const VecN<A,3>& _a,const VecN<B,3>& _b);

  /// output
  template <typename S,unsigned N>
  friend std::ostream& operator<<(std::ostream&,const VecN<S,N>&);
  /// input
  template <typename S,unsigned N>
  friend std::istream& operator>>(std::istream&,VecN<S,N>&);

  /// @}

  /** @name arithmetic: evaluate functions elementwise
      @{
  */

  /// replace each element `v[i]` by `_op(v[i])`
  template <typename Op>
  Self& map(Op& _op) {
    for (unsigned i=0;i<K;++i)
      m_v[i]=_op(m_v[i]);
    return *this;
  }
  /// replace each element `v[i]` by `_op(_v[i])`
  template <typename Op,typename S>
  Self& map(Op& _op,const VecN<S,K>& _v) {
    for (unsigned i=0;i<K;++i)
      m_v[i]=_op(_v.m_v[i]);
    return *this;
  }
  /// replace each element `v[i]` by `_op(_a[i],_b[i])`
  template <typename Op,typename A,typename B>
  Self& map(Op& _op,const VecN<A,K>& _a,const VecN<B,K>& _b) {
    for (unsigned i=0;i<K;++i)
      m_v[i]=_op(_a.m_v[i],_b.m_v[i]);
    return *this;
  }

  /// load random values (FixedAryOps::rand())
  Self& rand() {
    ::VC::math::rand(m_v,K);
    return *this;
  }

  /// @}


  /** @name Extrema, sorting and permutation.
      @{
  */

  /// index of minimum element (FixedAryOps::min())
  unsigned imin() const { return Ops::min(m_v); }
  /// minium element  (FixedAryOps::min())
  const T& min() const { return m_v[Ops::min(m_v)]; }
  /// minimum element, set index `_imin`  (FixedAryOps::min())
  const T& min(unsigned& _imin) const { return m_v[_imin=Ops::min(m_v)]; }

  /// load minimum elements of `_a` and `_b` to `*this`
  template <typename S1,typename S2>
  Self& min(const VecN<S1,K>& _a,const VecN<S2,K>& _b) {
    Ops::min(m_v,_a.m_v,_b.m_v);
    return *this;
  }

  /// index of maximum element (FixedAryOps::max())
  unsigned imax() const { return Ops::max(m_v); }
  /// maximum element (FixedAryOps::max())
  const T& max() const { return m_v[Ops::max(m_v)]; }
  /// maximum element, set index `_imax`  (FixedAryOps::max())
  const T& max(unsigned& _imax) const { return m_v[_imax=Ops::max(m_v)]; }

  /// load maximum elements of _a and _b to *this
  template <typename S1,typename S2>
  Self& max(const VecN<S1,K>& _a,const VecN<S2,K>& _b) {
    Ops::max(m_v,_a.m_v,_b.m_v);
    return *this;
  }

  /// indices of minimum and maximum elements (FixedAryOps::minmax())
  void minmax(unsigned* _minmax_idx) { Ops::minmax(_minmax_idx,m_v); }
  /// minimum and maximum elements (FixedAryOps::minmax())
  template <typename S>
  void minmax(S* _minmax) {
    unsigned idx[2];
    minmax(_minmax,idx);
  }
  /// minimum and maximum elements and indices (FixedAryOps::minmax())
  template <typename S>
  void minmax(S* _minmax,unsigned* _minmax_idx) {
    Ops::minmax(_minmax_idx,m_v);
    _minmax[0]=m_v[_minmax_idx[0]];
    _minmax[1]=m_v[_minmax_idx[1]];
  }

  /** Return selection/permutation `_idx` (select(),FixedAryOps::select()).
      Check `_idx` bounds if `VC_VECN_CHECK_BOUNDS` is defined.
  */
  template <typename Idx,unsigned L>
  VecN<T,L> operator()(const VecN<Idx,L>& _idx) {
    VecN<T,L> v;
# ifdef VC_VECN_CHECK_BOUNDS
    for (unsigned i=0;i<L;++i)
      assert(0<=int(_idx[i]) && int(_idx[i])<int(K));
# endif
    FixedAryOps<L>::select(v.data(),m_v,_idx.data());
    return v;
  }
  /** Load selection/permutation `_idx` of `_v` (FixedAryOps::select()).
      Requires `_v != *this` and integer type `_idx`.
      Check `_idx` bounds if `VC_VECN_CHECK_BOUNDS` is defined.
  */
  template <typename S,typename Idx,unsigned L>
  Self& select(const VecN<S,L>& _v,const VecN<Idx,K>& _idx) {
    assert(_v.data()!=m_v);
# ifdef VC_VECN_CHECK_BOUNDS
    for (unsigned i=0;i<K;++i)
      assert(0<=int(_idx[i]) && int(_idx[i])<int(L));
# endif
    Ops::select(m_v,_v.data(),_idx.data());
    return *this;
  }
  /** Get selection/apply permutation `_idx` to `*this` (FixedAryOps::select()).
      Check `_idx` bounds if `VC_VECN_CHECK_BOUNDS` is defined.
  */
  template <typename Idx>
  Self& select(const VecN<Idx,K>& _idx) {
    Self v(*this);
# ifdef VC_VECN_CHECK_BOUNDS
    for (unsigned i=0;i<K;++i)
      assert(0<=int(_idx[i]) && int(_idx[i])<int(K));
# endif
    Ops::select(m_v,v.m_v,_idx.data());
    return *this;
  }
  /** Load elements in selection `_idx` as linear elements from `_v` (FixedAryOps::assign()).
      Requires `_v != *this` and integer type `_idx`.
      Check `_idx` bounds if `VC_VECN_CHECK_BOUNDS` is defined.
  */
  template <typename S,typename Idx,unsigned L>
  Self& assign(const VecN<S,L>& _v,const VecN<Idx,L>& _idx) {
    assert(_v.data()!=m_v);
# ifdef VC_VECN_CHECK_BOUNDS
    for (unsigned i=0;i<K;++i)
      assert(0<=int(_idx[i]) && int(_idx[i])<int(K));
# endif
    FixedAryOps<L>::assign(m_v,_v.data(),_idx.data());
    return *this;
  }

  /// get permutation sorting elements in ascending order (FixedAryOps::isort())
  VecN<unsigned,K> isort() const {
    VecN<unsigned,K> idx;
    Ops::isort(idx.data(),m_v);
    return idx;
  }

  /// get permutation _idx sorting elements in ascending order (FixedAryOps::isort())
  template <typename Idx>
  Self& isort(VecN<Idx,K>& _idx) {
    Ops::isort(_idx.data(),m_v);
    return *this;
  }

  /// sort elements in ascending order (FixedAryOps::isort())
  Self& sort() {
    Self v(*this);
    VecN<unsigned,K> idx;
    Ops::isort(idx.data(),m_v);
    Ops::select(m_v,v.m_v,idx.data());
    return *this;
  }

  /// sort elements in ascending order and yield permutation `_idx` (FixedAryOps::isort())
  template <typename Idx>
  Self& sort(VecN<Idx,K>& _idx) {
    Self v(*this);
    Ops::isort(_idx.data(),m_v);
    Ops::select(m_v,v.m_v,_idx.data());
    return *this;
  }

  /// @}


private:
  union {
    T m_v[K];                //!< data fields
    VecN_align<T,K> m_align; //!< force alignment
  };
};


/** \defgroup vc_math_lav_vops Linear algebra: externally defined operations on VecN
    \ingroup vc_math_lav
    \sa VecN
*/

/** @name Externally defined operations on VecN
    @{
*/

/// vector addition `_a+_b` (addv()) \ingroup  vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator+(const VecN<A,N>& _a,const VecN<B,N>& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::addv(c.m_v, _a.m_v,_b.m_v);
  return c;
}
/// vector subtraction `_a-_b` (subv()) \ingroup  vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator-(const VecN<A,N>& _a,const VecN<B,N>& _b) {
  VecN<typename Result<A,B>::type,N> c;
  FixedAryOps<N>::subv(c.m_v, _a.m_v,_b.m_v);
  return c;
}
/// \a element-wise multiplication (mulv()) \ingroup  vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator*(const VecN<A,N>& _a,const VecN<B,N>& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::mulv(c.m_v, _a.m_v,_b.m_v);
  return c;
}
/// \a element-wise division (divv()) \ingroup  vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator/(const VecN<A,N>& _a,const VecN<B,N>& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::divv(c.m_v, _a.m_v,_b.m_v);
  return c;
}
/// unary minus: return `-_v` (VecN<N,t>::negate()) \ingroup vc_math_lav_vops
template <typename T,unsigned N>
VecN<T,N>
operator-(const VecN<T,N>& _v) {
  VecN<T,N> v=_v;
  FixedAryOps<N>::muls(v.m_v, _v.m_v,T(-1));
  return v;
}

/// add scalar `_b` to all components of `_a` (adds()) \ingroup vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator+(const VecN<A,N>& _a,const B& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::adds(c.m_v, _a.m_v,_b);
  return c;
}
/// add scalar `_a` to all components of `_b` (adds()) \ingroup vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator+(const A& _a,const VecN<A,N>& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::adds(c.m_v, _b.m_v,_a);
  return c;
}
/// subtract scalar `_b` from all components of `_a` (subs()) \ingroup vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator-(const VecN<A,N>& _a,const B& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::subs(c.m_v, _a.m_v,_b);
  return c;
}
/// multiply scalar `_b` to all components of `_a` (muls()) \ingroup vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator*(const VecN<A,N>& _a,const B& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::muls(c.m_v, _a.m_v,_b);
  return c;
}
/// divide components of `_a` by `_b` (divs()) \ingroup vc_math_lav_vops
template <typename A,typename B,unsigned N>
VecN<typename Result<A,B>::type,N>
operator/(const VecN<A,N>& _a,const B& _b) {
  VecN<typename  Result<A,B>::type,N> c;
  FixedAryOps<N>::divs(c.m_v, _a.m_v,_b);
  return c;
}

/// cross product of 3-vectors `_a x _b`
template <typename A,typename B>
VecN<typename Result<A,B>::type,3>
operator%(const VecN<A,3>& _a,const VecN<B,3>& _b) {
  return VecN<typename  Result<A,B>::type,3>
    (_a.m_v[1]*_b.m_v[2]-_a.m_v[2]*_b.m_v[1],
     _a.m_v[2]*_b.m_v[0]-_a.m_v[0]*_b.m_v[2],
     _a.m_v[0]*_b.m_v[1]-_a.m_v[1]*_b.m_v[0]);
}

// // begin: using rvalues

// ambiguous overloads -- fix is nontrivial (Result<>)

// template <typename T,unsigned N>
// VecN<T,N> operator+(VecN<T,N>&& _a,const VecN<T,N>& _b) {
//   VecN<T,N> c=_a; c+=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator+(const VecN<T,N>& _b,VecN<T,N>&& _a) {
//   VecN<T,N> c=_a; c+=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator-(VecN<T,N>&& _a,const VecN<T,N>& _b) {
//   VecN<T,N> c=_a; c+=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator-(const VecN<T,N>& _b,VecN<T,N>&& _a) {
//   VecN<T,N> c=_a; c+=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator*(VecN<T,N>&& _a,const VecN<T,N>& _b) {
//   VecN<T,N> c=_a; c*=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator*(const VecN<T,N>& _b,VecN<T,N>&& _a) {
//   VecN<T,N> c=_a; c/=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator/(VecN<T,N>&& _a,const VecN<T,N>& _b) {
//   VecN<T,N> c=_a; c*=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator/(const VecN<T,N>& _b,VecN<T,N>&& _a) {
//   VecN<T,N> c=_a; c/=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator%(VecN<T,N>&& _a,const VecN<T,N>& _b) {
//   VecN<T,N> c=_a; c*=_b; return std::move(c);
// }
// template <typename T,unsigned N>
// VecN<T,N> operator%(const VecN<T,N>& _b,VecN<T,N>&& _a) {
//   VecN<T,N> c=_a; c/=_b; return std::move(c);
// }

// // end: using rvalues


/// Euclidean norm `sqrt(sqr(_v))`
template <typename T,unsigned N>
T norm(const VecN<T,N>& _v) { return FixedAryOps<N>::norm(_v.m_v); }

/// 1-norm `sum(abs(_v))`
template <typename T,unsigned N>
T norm1(const VecN<T,N>& _v) { return FixedAryOps<N>::norm1(_v.m_v); }

/// 2-norm, synonym for norm()
template <typename T,unsigned N>
T norm2(const VecN<T,N>& _v) { return FixedAryOps<N>::norm(_v.m_v); }

/// maximum norm VecN::norminf()
template <typename T,unsigned N>
T norminf(const VecN<T,N>& _v) { return _v.m_v[FixedAryOps<N>::max(_v.m_v)]; }


/// swap contents of `_a` and `_b` \ingroup vc_math_lav_vops
template <typename T,unsigned N>
void swap(VecN<T,N>& _a,VecN<T,N>& _b) {
  FixedAryOps<N>::swap(_a.m_v,_b.m_v);
}

/// external definition of `VecN::to_v()` \ingroup vc_math_lav_vops
template <unsigned K,typename T>
VecN<T,K>& to_v(T* _v) {
  return (VecN<T,K>&) *_v;
}
/// external definition of `VecN::to_cv()` \ingroup vc_math_lav_vops
template <unsigned K,typename T>
const VecN<T,K>& to_cv(const T* _v) {
  return *((const VecN<T,K>*) _v);
}

/** Concatenate/stack vectors `_v0` and `_v1` vertically. \ingroup vc_math_lav_vops
    Creates a new `(K0+K1)`-vector `[_v0;_v1]` with `_v0`'s scalar type.
*/
template <typename T0,unsigned K0,typename T1,unsigned K1>
VecN<T0,K0+K1> v_cat(const VecN<T0,K0>& _v0,const VecN<T1,K1>& _v1) {
  VecN<T0,K0+K1> v;
  unsigned i;
  for (i=0;i<K0   ;++i)  v.m_v[i]=   _v0.m_v[i];
  for (   ;i<K0+K1;++i)  v.m_v[i]=T0(_v1.m_v[i-K0]);
  return v;
}

/** Concatenate/stack vectors `_v0`, `_v1` and `_v2` vertically. \ingroup vc_math_lav_vops
    Creates a new `(K0+K1_K2)`-vector `[_v0;_v1;_v2]` with `_v0`'s scalar type.
*/
template <typename T0,unsigned K0,typename T1,unsigned K1,typename T2,unsigned K2>
VecN<T0,K0+K1+K2> v_cat(const VecN<T0,K0>& _v0,const VecN<T1,K1>& _v1,
                        const VecN<T2,K2>& _v2) {
  VecN<T0,K0+K1+K2> v;
  unsigned i;
  for (i=0;i<K0      ;++i) v.m_v[i]=   _v0.m_v[i];
  for (   ;i<K0+K1   ;++i) v.m_v[i]=T0(_v1.m_v[i-K0]);
  for (   ;i<K0+K1+K2;++i) v.m_v[i]=T0(_v2.m_v[i-K0-K1]);
  return v;
}

/// get angle between vectors `_a` and `_b` (calls angle_u()) \ingroup vc_math_lav_vops
template <typename A,typename B,unsigned N>
typename Result<A,B>::type
angle(const VecN<A,N>& _a,const VecN<B,N>& _b) {
  return angle_u(_a/_a.norm(),_b/_b.norm());
}
/** Get angle between _unit vectors `_a` and `_b`.
    \ingroup vc_math_lav_vops
    Clamps if `0<|acos(<_a|_b>)-1|<eps`
*/
template <typename A,typename B,unsigned N>
typename Result<A,B>::type
angle_u(const VecN<A,N>& _a,const VecN<B,N>& _b) {
  double d(_a|_b); // force double
  assert(fabs(d)<1.0+std::numeric_limits<double>::epsilon()*1000.0);
  if      (d<-1.0) d=-1.0;
  else if (d>+1.0) d=+1.0;
  return acos(d);
}

/// get determinant of 2x2 matrix `[_a,_b]` \ingroup vc_math_lav_vops
template <typename T>
T det(const VecN<T,2>& _a,const VecN<T,2>& _b) {
  return _a[0]*_b[1]-_a[1]*_b[0];
}
/// get determinant of 3x3 matrix `[_a,_b,_c]` \ingroup vc_math_lav_vops
template <typename T>
T det(const VecN<T,3>& _a,const VecN<T,3>& _b,const VecN<T,3>& _c) {
  return (_a|(_b%_c));
}

// h_cat creates matrices

// perp, external dot (+dotc taking diagonal matrix), rot

// utility: setup local coordinate system (N+1 pts)

// IO + format (pmatlab)

/*
  compare (lexicographic)
  hash (?!)

  inverse permutation

  left/right shift, rotate

*/

/** Stream output of VC::math::VecN.
    \ingroup vc_math_lav
*/
template <typename T,unsigned K>
std::ostream& operator<<(std::ostream& _s,const ::VC::math::VecN<T,K>& _v) {
  for (unsigned i=0;i<K-1;++i)
    _s << _v[i] << ' ';
  _s << _v[K-1];
  return _s;
}
/** Stream output of VC::math::VecN.
    \ingroup vc_math_lav
*/
template <typename T,unsigned K>
std::istream& operator>>(std::istream& _s,::VC::math::VecN<T,K>& _v) {
  for (unsigned i=0;i<K;++i) _s >> _v[i];
  return _s;
}

/// @}

//-----------------------------------------------------------------------------

} // namespace math
} // namespace VC


#endif // __VC_MATH_VECN_HH__
