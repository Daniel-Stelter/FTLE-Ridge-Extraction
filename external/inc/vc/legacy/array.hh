//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: array.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================

#ifndef __VC_MATH_ARRAY_HH
#define __VC_MATH_ARRAY_HH


#include <cmath>
#include <cstdlib>
#include <cstring>

// #include "../system/platform.hh" // drand48

#include <algorithm>

namespace VC {
namespace math {

/** \file array.hh
      math operations on arrays, [\ref vc_math_lav].
  */

/** \defgroup vc_math_lav Linear algebra: fixed-size vectors
    \ingroup vc_math
    Data types, arithmetic and utility functions for vector-vector
    and vector-scalar operations.
    \arg VecN <T,N> defines a N-vector data type.
    \arg blas::vector_reference_t defines an interface to using BLAS
     operations on vectors.
    \arg blas::Vec_N <T> defines a N-vector using the blas::vector_reference_t interface.
    \arg blas::Vector <T> defines variable size vector based on blas::vector_reference_t.
    \arg Core functions addv(), etc. work on general n-arrays.
    \arg Result < TypeA,TypeB > determines type of result.
    \arg FixedAryOps < N > define operations for fixed N, possibly
    specialized for small N.
    \sa VecN, addv(), Result, FixedAryOps, [\ref vc_blas]

    \todo specialization using SIMD operations for dimensions 2 and 4
    \todo use __restrict
    \todo multiply_add (how to integrate this in VecN?)
    \todo left-/right-rotate by n (low priority: we got select)
    \todo matrix multiplication, dyadic product
    \todo ddot accumulate double, also for VecN
    \todo LD versions

*/

#  define _LOOP0(n,expr) for (unsigned i=0;i<(n);++i) { expr; }
#  define _LOOP1(n,expr) for (unsigned i=1;i<(n);++i) { expr; }
#  define _LOOP1COND(n,cond,expr) \
  for (unsigned i=1;i<(n);++i) { if (cond) { expr; } }
#  define _LOOP1COND2(n,cond1,cond2,expr1,expr2) \
  for (unsigned i=1;i<(n);++i) { if (cond1) { expr1; } else if (cond2) { expr2; } }

  /** @name core functions to manipulate arrays
      @{
      \arg Functions which return vector values take the \a destination as
      \a first argument.
      \arg Functions which need to determine the type of the return value use
      ::VC::math::Result<Type1,Type2>.
      \arg The functions are also \c static memebers of
      ::VC::math::FixedAryOps<N>, which may be specialized for small N.
      \ingroup vc_math_lav
  */

  /// determine type of result (e.g. float+double=double by specialization)
  template <typename A,typename B> struct Result { typedef A type; };
# ifndef DOXYGEN_SKIP
  template <typename B> struct Result<double,B> { typedef double type; };
  template <typename B> struct Result<float,B> { typedef float type; };
  template <> struct Result<double,float> { typedef double type; };
# endif

  //
  // Note: gcc won't unroll loops nicely when using pointers!
  //

  /// elementwise c=a+b \ingroup vc_math_lav
  template <typename A,typename B,typename C>

  inline void addv(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]+_b[i]));
  }
  /// elementwise c=a-b \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void subv(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]-_b[i]));
  }
  /// elementwise c=a*b \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void mulv(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]*_b[i]));
  }
  /// elementwise c=a/b \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void divv(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]/_b[i]));
  }

  /// elementwise c=a+b (scalar b) \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void adds(C* _c,const A* _a,const B& _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]+_b));
  }
  /// elementwise c=a-b (scalar b) \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void subs(C* _c,const A* _a,const B& _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]-_b));
  }
  /// elementwise c=a-b (scalar a) \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void subs2(C* _c,const A& _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a-_b[i]));
  }
  /// elementwise c=a*b (scalar b) \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void muls(C* _c,const A* _a,const B& _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]*_b));
  }
  /// elementwise c=a/b (scalar b) \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void divs(C* _c,const A* _a,const B& _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]/_b));
  }
  /// elementwise c=a/b (scalar a) \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void divs2(C* _c,const A& _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a/_b[i]));
  }

  /// sum \ingroup vc_math_lav
  template <typename A>
  A sum(const A* _a,unsigned _n) {
    A x=_a[0];
    _LOOP1(_n,(x+=_a[i]));
    return x;
  }

  /// dot product <a|b> \ingroup vc_math_lav
  template <typename A,typename B>
  typename Result<A,B>::type dot(const A* _a,const B* _b,unsigned _n) {
    typename Result<A,B>::type x=_a[0]*_b[0];
    _LOOP1(_n,(x+=_a[i]*_b[i]));
    return x;
  }
  /// <a|a> = ||a||_2^2
  template <typename A>
  A sqr(const A* _a,unsigned _n) {
    A x=_a[0]*_a[0];
    _LOOP1(_n,(x+=_a[i]*_a[i]));
    return x;
  }
  /// Euclidean norm ||a||_2 \ingroup vc_math_lav
  template <typename A>
  A norm(const A* _a,unsigned _n) { return sqrt(sqr(_a,_n)); }
  /// city block norm ||a||_1=sum(fabs(a)) \ingroup vc_math_lav
  template <typename A>
  A norm1(const A* _a,unsigned _n) {
    A x=fabs(_a[0]);
    _LOOP1(_n,(x+=fabs(_a[i])));
    return x;
  }
  /// infinity norm ||a||_inf=max{fabs(a[i])} \ingroup vc_math_lav
  template <typename A>
  A norminf(const A* _a,unsigned _n) {
    /*register*/ A x=fabs(_a[0]),y;
    _LOOP1COND(_n,((y=fabs(_a[i]))>x),(x=y));
    return x;
  }
  /// mean (floating point only) \ingroup math_
  inline double mean(const double* _a,unsigned _n) { return sum(_a,_n)/double(_n); }
  inline float mean(const float* _a,unsigned _n) { return sum(_a,_n)/float(_n); }

  /// elementwise swap(_a,_b) \ingroup vc_math_lavr
  template <typename T>
  inline void swap(T* _a,T* _b,unsigned _n) {
    /*register*/ T c;
    _LOOP0(_n,(c=_a[i], _a[i]=_b[i], _b[i]=c));
  }

  /// map operator Op to all elements _a = _op(_b) \ingroup vc_math_lav
  template <typename Op,typename A,typename B>
  inline void map(Op _op,A* _a, const B* _b,unsigned _n) {
    _LOOP0(_n,(_a[i]=_op(_b[i])));
  }
  /// map operator Op to all elements _c = _op(_a,_b) \ingroup vc_math_lav
  template <typename Op,typename A,typename B,typename C>
  inline void map(Op _op,C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_op(_a[i],_b[i])));
  }

  /// linear combination c=a*sa+b*sb \ingroup vc_math_lav
  template <typename A,typename B,typename C>
  inline void lincomb(C* _c,
                      const A& _sa,const A* _a,
                      const B& _sb,const B* _b, unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[i]*_sa+_b[i]*_sb));
  }
  /// linear combination c=v1*s1+v2*s2+v3*s3 \ingroup vc_math_lav
  template <typename T1,typename T2,typename T3,typename C>
  inline void lincomb(C* _c,
                      const T1& _s1,const T1* _v1,
                      const T2& _s2,const T2* _v2,
                      const T3& _s3,const T3* _v3, unsigned _n) {
    _LOOP0(_n,(_c[i]=_v1[i]*_s1+_v2[i]*_s2+_v3[i]*_s3));
  }
  /// linear combination c=v1*s1+v2*s2+v3*s3+v4*s4 \ingroup vc_math_lav
  template <typename T1,typename T2,typename T3,typename T4,typename C>
  inline void lincomb(C* _c,
                      const T1& _s1,const T1* _v1,
                      const T2& _s2,const T2* _v2,
                      const T3& _s3,const T3* _v3,
                      const T4& _s4,const T4* _v4, unsigned _n) {
    _LOOP0(_n,(_c[i]=_v1[i]*_s1+_v2[i]*_s2+_v3[i]*_s3+_v4[i]*_s4));
  }

  /// get index of minimal element (first occurence) \ingroup vc_math_lav
  template <typename T>
  unsigned min(const T* _v,unsigned _n) {
    T        t=_v[0];
    unsigned imin=0;
    _LOOP1COND(_n,(_v[i]<t),(t=_v[i],imin=i));
    return imin;
  }
  /// get index of maximal element (first occurence) \ingroup vc_math_lav
  template <typename T>
  unsigned max(const T* _v,unsigned _n) {
    T        t=_v[0];
    unsigned imax=0;
    _LOOP1COND(_n,(_v[i]>t),(t=_v[i],imax=i));
    return imax;
  }
  /** Get index of minimal and maximal elements (first occurences).
      \param _minmax_idx array of size 2 storing indices to min and max
      \param _v vector
      \param _n length of _v
      \ingroup vc_math_lav
   */
  template <typename T>
  void minmax(unsigned* _minmax_idx,const T* _v,unsigned _n) {
    T        tmin=_v[0], tmax=_v[0];
    unsigned imax=0, imin=0;

    _LOOP1COND2(_n,(_v[i]<tmin),(_v[i]>tmax),(tmin=_v[imin=i]),(tmax=_v[imax=i]));

    _minmax_idx[0]=imin;
    _minmax_idx[1]=imax;
  }
  /// element-wise copy of minimum values
  template <typename A,typename B,typename C>
  void min(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=(_a[i]<=_b[i]) ? _a[i] : _b[i]));
  }
  /// element-wise copy maximum values
  template <typename A,typename B,typename C>
  void max(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=(_a[i]>=_b[i]) ? _a[i] : _b[i]));
  }

  /// copy elements selected by index _b: _c[i]=_a[_b[i]] (B is an integer type)
  template <typename A,typename B,typename C>
  void select(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[i]=_a[_b[i]]));
  }
  /// assign elements selected by index _b: _c[_b[i]]=_a[i] (B is an integer type)
  template <typename A,typename B,typename C>
  void assign(C* _c,const A* _a,const B* _b,unsigned _n) {
    _LOOP0(_n,(_c[_b[i]]=_a[i]));
  }

# ifndef DOXYGEN_SKIP
  /// (used internally by isort()) \internal
  template <typename A,typename B>
  struct _isort_compare {
    const A* base;
    _isort_compare(const A* _base) : base(_base) {}
    bool operator()(const B& _a,const B& _b) const {
      return base[unsigned(_a)]<base[unsigned(_b)];
    }
  };
# endif // DOXYGEN_SKIP

  /** Find permutation _b of _a such that _a[_b[i]] is sorted.
      \arg The sequence _a[_b[i]] is sorted in \a ascending order.
      \arg Sorting is \a not neccessarily stable.
   */
  template <typename A,typename B>
  void isort(B* _b,const A* _a,unsigned _n) {
    for (unsigned i=0;i<_n;++i) _b[i]=B(i);
    _isort_compare<A,B> cmp(_a);
    std::sort(_b,_b+_n,cmp);
  }

  /// load _n random values in [0,1] to _a (use drand48())
  inline void rand(double* _a,unsigned _n) {
    _LOOP0(_n,(_a[i]=drand48()));
  }
  /// load _n random values in [0,1] to _a (use drand48())
  inline void rand(float* _a,unsigned _n) {
    _LOOP0(_n,(_a[i]=float(drand48())));
  }

  // BLAS


# define FixedAryOps_abc(x) \
  template <typename A,typename B,typename C>\
  static void x(C* _c,const A* _a,const B* _b) { ::VC::math::x(_c,_a,_b,N); }
# define FixedAryOps_abcs(x) \
  template <typename A,typename B,typename C>\
  static void x(C* _c,const A* _a,const B& _b) { ::VC::math::x(_c,_a,_b,N); }

  /** \class FixedAryOps_base array.hh
      \brief (internal use as default implementation of FixedAryOps)
      \ingroup vc_math_lav
  */
  template <unsigned N>
  struct FixedAryOps_base {
    FixedAryOps_abc(addv)
    FixedAryOps_abc(subv)
    FixedAryOps_abc(mulv)
    FixedAryOps_abc(divv)

    FixedAryOps_abc(select)
    FixedAryOps_abc(assign)

    FixedAryOps_abcs(adds)
    FixedAryOps_abcs(subs)
    FixedAryOps_abcs(muls)
    FixedAryOps_abcs(divs)

    template <typename A,typename B,typename C>
    static void subs2(C* _c,const A& _a,const B* _b) { VC::math::subs2(_c,_a,_b,N); }
    template <typename A,typename B,typename C>
    static void divs2(C* _c,const A& _a,const B* _b) { VC::math::divs2(_c,_a,_b,N); }

    template <typename A,typename B,typename C>
    static void lincomb(C* _c,
                        const A& _sa,const A* _a,
                        const B& _sb,const B* _b) {
      ::VC::math::lincomb(_c,_sa,_a,_sb,_b,N);
    }
    template <typename T1,typename T2,typename T3,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3) {
      ::VC::math::lincomb(_c,_s1,_v1,_s2,_v2,_s3,_v3,N);
    }
    template <typename T1,typename T2,typename T3,typename T4,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3,
                        const T4& _s4,const T4* _v4) {
      ::VC::math::lincomb(_c,_s1,_v1,_s2,_v2,_s3,_v3,_s4,_v4,N);
    }
    template <typename A>
    static A sum(const A* _a) { return ::VC::math::sum(_a,N); }

    template <typename A,typename B>
    static typename Result<A,B>::type dot(const A* _a,const B* _b) {
      return ::VC::math::dot(_a,_b,N);
    }
    template <typename A>
    static A sqr(const A* _a) { return ::VC::math::sqr(_a,N); }

    template <typename A>
    static A norm(const A* _a) { return ::VC::math::norm(_a,N); }

    template <typename A>
    static A norm1(const A* _a) { return ::VC::math::norm1(_a,N); }

    template <typename A>
    static A norminf(const A* _a) { return ::VC::math::norminf(_a,N); }

    static double mean(const double* _a) { return ::VC::math::mean(_a,N); }
    static float mean(const float* _a) { return ::VC::math::mean(_a,N); }

    template <typename T>
    static void swap(T* _a,T* _b)  { ::VC::math::swap(_a,_b,N); }

    template <typename Op,typename A,typename B>
    static void map(Op _op,A* _a, const B* _b) { ::VC::math::map(_op,_a,_b,N); }
    template <typename Op,typename A,typename B,typename C>
    static void map(Op _op,C* _c,const A* _a,const B* _b) {
      ::VC::math::map(_op,_c,_a,_b,N);
    }

    template <typename T>
    static unsigned min(const T* _v) { return ::VC::math::min(_v,N); }
    template <typename T>
    static unsigned max(const T* _v) { return ::VC::math::max(_v,N); }
    template <typename T>
    static void minmax(unsigned* _minmax_idx,const T* _v) {
      ::VC::math::minmax(_minmax_idx,_v,N);
    }
    template <typename A,typename B,typename C>
    static void min(C* _c,const A* _a,const B* _b) {
      ::VC::math::min(_c,_a,_b,N);
    }
    template <typename A,typename B,typename C>
    static void max(C* _c,const A* _a,const B* _b) {
      ::VC::math::max(_c,_a,_b,N);
    }
    template <typename A,typename B>
    static void isort(B* _b,const A* _a) {
      ::VC::math::isort(_b,_a,N);
    }
  };
# undef FixedAryOps_abc
# undef FixedAryOps_abcs

/** \class FixedAryOps array.hh
    \brief operations on fixed sized arrays
    \ingroup vc_math_lav


    This class is specialized for compiler optimization for small
    dimensions N.
*/
  template <unsigned N>
  struct FixedAryOps : public FixedAryOps_base<N> {
  };

# ifndef DOXYGEN_SKIP

  // N=1
# define FixedAryOps_abcv(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B* _b) { \
  _c[0]=_a[0] op _b[0]; }
# define FixedAryOps_abcs(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B& _b) { \
  _c[0]=_a[0] op _b;  }

  template <>
  struct FixedAryOps<1> {
    FixedAryOps_abcv(addv,+)
    FixedAryOps_abcv(subv,-)
    FixedAryOps_abcv(mulv,*)
    FixedAryOps_abcv(divv,/)

    FixedAryOps_abcs(adds,+)
    FixedAryOps_abcs(subs,-)
    FixedAryOps_abcs(muls,*)
    FixedAryOps_abcs(divs,/)

    template <typename A,typename B,typename C>
    static void subs2(C* _c,const A& _a,const B* _b) { _c[0]=_a-_b[0]; }
    template <typename A,typename B,typename C>
    static void divs2(C* _c,const A& _a,const B* _b) { _c[0]=_a/_b[0]; }

    template <typename A,typename B,typename C>
    static void lincomb(C* _c,
                               const A& _sa,const A* _a,
                               const B& _sb,const B* _b) {
      _c[0]=_a[0]*_sa+_b[0]*_sb;
    }
    template <typename T1,typename T2,typename T3,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3;
    }
    template <typename T1,typename T2,typename T3,typename T4,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3,
                        const T4& _s4,const T4* _v4) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3+_v4[0]*_s4;
    }
    template <typename A>
    static A sum(const A* _a) { return _a[0]; }

    template <typename A,typename B>
    static typename Result<A,B>::type dot(const A* _a,const B* _b) {
      typename Result<A,B>::type x=_a[0]*_b[0];
      return x;
    }
    template <typename A>
    static A sqr(const A* _a) { return _a[0]*_a[0]; }

    template <typename A>
    static A norm(const A* _a) { return fabs(_a[0]); }

    template <typename A>
    static A norm1(const A* _a) { return fabs(_a[0]); }

    template <typename A>
    static A norminf(const A* _a) { return fabs(_a[0]); }

    static double mean(const double* _a) { return _a[0]; }
    static float mean(const float* _a) { return _a[0]; }

    template <typename T>
    static void swap(T* _a,T* _b) {
      /*register*/ T c;
      c=_a[0]; _a[0]=_b[0]; _b[0]=c;
    }
    template <typename Op,typename A,typename B>
    static void map(Op _op,A* _a, const B* _b) {
      _a[0]=_op(_b[0]);
    }
    template <typename Op,typename A,typename B,typename C>
    static void map(Op _op,C* _c,const A* _a,const B* _b) {
      _c[0]=_op(_a[0],_b[0]);
    }

    template <typename T>
    static unsigned min(const T*) { return 0; }
    template <typename T>
    static unsigned max(const T*) { return 0; }
    template <typename T>
    static void minmax(unsigned* _minmax_idx,const T*) {
      _minmax_idx[0]=_minmax_idx[1]=0;
    }
    template <typename A,typename B,typename C>
    static void min(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]<=_b[0]) ? _a[0] : _b[0];
    }
    template <typename A,typename B,typename C>
    static void max(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]>=_b[0]) ? _a[0] : _b[0];
    }

    template <typename A,typename B,typename C>
    static void select(C* _c,const A* _a,const B* _b) {
      _c[0]=_a[_b[0]];
    }
    template <typename A,typename B,typename C>
    static void assign(C* _c,const A* _a,const B* _b) {
      _c[_b[0]]=_a[0];
    }
    template <typename A,typename B>
    static void isort(B* _b,const A*) {
      _b[0]=B(0);
    }
  };
# undef FixedAryOps_abcv
# undef FixedAryOps_abcs

  // N=2
# define FixedAryOps_abcv(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B* _b) { \
  _c[0]=_a[0] op _b[0]; _c[1]=_a[1] op _b[1]; }
# define FixedAryOps_abcs(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B& _b) { \
  _c[0]=_a[0] op _b; _c[1]=_a[1] op _b; }

  template <>
  struct FixedAryOps<2> {
    FixedAryOps_abcv(addv,+)
    FixedAryOps_abcv(subv,-)
    FixedAryOps_abcv(mulv,*)
    FixedAryOps_abcv(divv,/)

    FixedAryOps_abcs(adds,+)
    FixedAryOps_abcs(subs,-)
    FixedAryOps_abcs(muls,*)
    FixedAryOps_abcs(divs,/)

    template <typename A,typename B,typename C>
    static void subs2(C* _c,const A& _a,const B* _b) {
      _c[0]=_a-_b[0]; _c[1]=_a-_b[1];
    }
    template <typename A,typename B,typename C>
    static void divs2(C* _c,const A& _a,const B* _b) {
      _c[0]=_a/_b[0]; _c[1]=_a/_b[1];
    }
    template <typename A,typename B,typename C>
    static void lincomb(C* _c,
                               const A& _sa,const A* _a,
                               const B& _sb,const B* _b) {
      _c[0]=_a[0]*_sa+_b[0]*_sb;
      _c[1]=_a[1]*_sa+_b[1]*_sb;
    }
    template <typename T1,typename T2,typename T3,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3;
      _c[1]=_v1[1]*_s1+_v2[1]*_s2+_v3[1]*_s3;
    }
    template <typename T1,typename T2,typename T3,typename T4,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3,
                        const T4& _s4,const T4* _v4) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3+_v4[0]*_s4;
      _c[1]=_v1[1]*_s1+_v2[1]*_s2+_v3[1]*_s3+_v4[1]*_s4;
    }
    template <typename A>
    static A sum(const A* _a) { return _a[0]+_a[1]; }

    template <typename A,typename B>
    static typename Result<A,B>::type dot(const A* _a,const B* _b) {
      typename Result<A,B>::type x=_a[0]*_b[0]+_a[1]*_b[1];
      return x;
    }
    template <typename A>
    static A sqr(const A* _a) { return _a[0]*_a[0]+_a[1]*_a[1]; }

    template <typename A>
    static A norm(const A* _a) { return sqrt(_a[0]*_a[0]+_a[1]*_a[1]);  }

    template <typename A>
    static A norm1(const A* _a) { return fabs(_a[0])+fabs(_a[1]); }

    template <typename A>
    static A norminf(const A* _a) { return std::max(fabs(_a[0]),fabs(_a[1])); }

    static double mean(const double* _a) { return (_a[0]+_a[1])*0.5f; }
    static float mean(const float* _a) { return (_a[0]+_a[1])*0.5f; }

    template <typename T>
    static void swap(T* _a,T* _b) {
      /*register*/ T c;
      c=_a[0]; _a[0]=_b[0]; _b[0]=c;
      c=_a[1]; _a[1]=_b[1]; _b[1]=c;
    }
    template <typename Op,typename A,typename B>
    static void map(Op _op,A* _a, const B* _b) {
      _a[0]=_op(_b[0]); _a[1]=_op(_b[1]);
    }
    template <typename Op,typename A,typename B,typename C>
    static void map(Op _op,C* _c,const A* _a,const B* _b) {
      _c[0]=_op(_a[0],_b[0]); _c[1]=_op(_a[1],_b[1]);
    }
    template <typename T>
    static unsigned min(const T* _v) { return (_v[0]<=_v[1]) ? 0 : 1; }
    template <typename T>
    static unsigned max(const T* _v) { return (_v[0]>=_v[1]) ? 0 : 1; }
    template <typename T>
    static void minmax(unsigned* _minmax_idx,const T* _v) {
      if (_v[0]<=_v[1]) { _minmax_idx[0]=0; _minmax_idx[1]=1; }
      else              { _minmax_idx[0]=1; _minmax_idx[1]=0; }
    }
    template <typename A,typename B,typename C>
    static void min(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]<=_b[0]) ? _a[0] : _b[0];
      _c[1]=(_a[1]<=_b[1]) ? _a[1] : _b[1];
    }
    template <typename A,typename B,typename C>
    static void max(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]>=_b[0]) ? _a[0] : _b[0];
      _c[1]=(_a[1]>=_b[1]) ? _a[1] : _b[1];
    }

    template <typename A,typename B,typename C>
    static void select(C* _c,const A* _a,const B* _b) {
      _c[0]=_a[_b[0]]; _c[1]=_a[_b[1]];
    }
    template <typename A,typename B,typename C>
    static void assign(C* _c,const A* _a,const B* _b) {
      _c[_b[0]]=_a[0]; _c[_b[1]]=_a[1];
    }

    template <typename A,typename B>
    static void isort(B* _b,const A* _a) {
      if (_a[0]<=_a[1]) {
        _b[0]=B(0); _b[1]=B(1);
      }
      else {
        _b[0]=B(1); _b[1]=B(0);
      }
    }
  };
# undef FixedAryOps_abcv
# undef FixedAryOps_abcs

  // N=3
# define FixedAryOps_abcv(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B* _b) { \
  _c[0]=_a[0] op _b[0]; _c[1]=_a[1] op _b[1]; _c[2]=_a[2] op _b[2]; }
# define FixedAryOps_abcs(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B& _b) { \
  _c[0]=_a[0] op _b; _c[1]=_a[1] op _b; _c[2]=_a[2] op _b; }

  template <>
  struct FixedAryOps<3> {
    FixedAryOps_abcv(addv,+)
    FixedAryOps_abcv(subv,-)
    FixedAryOps_abcv(mulv,*)
    FixedAryOps_abcv(divv,/)

    FixedAryOps_abcs(adds,+)
    FixedAryOps_abcs(subs,-)
    FixedAryOps_abcs(muls,*)
    FixedAryOps_abcs(divs,/)

    template <typename A,typename B,typename C>
    static void subs2(C* _c,const A& _a,const B* _b) {
      _c[0]=_a-_b[0]; _c[1]=_a-_b[1]; _c[2]=_a-_b[2];
    }
    template <typename A,typename B,typename C>
    static void divs2(C* _c,const A& _a,const B* _b) {
      _c[0]=_a/_b[0]; _c[1]=_a/_b[1]; _c[2]=_a/_b[2];
    }
    template <typename A,typename B,typename C>
    static void lincomb(C* _c,
                               const A& _sa,const A* _a,
                               const B& _sb,const B* _b) {
      _c[0]=_a[0]*_sa+_b[0]*_sb;
      _c[1]=_a[1]*_sa+_b[1]*_sb;
      _c[2]=_a[2]*_sa+_b[2]*_sb;
    }
    template <typename T1,typename T2,typename T3,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3;
      _c[1]=_v1[1]*_s1+_v2[1]*_s2+_v3[1]*_s3;
      _c[2]=_v1[2]*_s1+_v2[2]*_s2+_v3[2]*_s3;
    }
    template <typename T1,typename T2,typename T3,typename T4,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3,
                        const T4& _s4,const T4* _v4) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3+_v4[0]*_s4;
      _c[1]=_v1[1]*_s1+_v2[1]*_s2+_v3[1]*_s3+_v4[1]*_s4;
      _c[2]=_v1[2]*_s1+_v2[2]*_s2+_v3[2]*_s3+_v4[2]*_s4;
    }
    template <typename A>
    static A sum(const A* _a) { return _a[0]+_a[1]+_a[2]; }

    template <typename A,typename B>
    static typename Result<A,B>::type dot(const A* _a,const B* _b) {
      typename Result<A,B>::type x=_a[0]*_b[0]+_a[1]*_b[1]+_a[2]*_b[2];
      return x;
    }
    template <typename A>
    static A sqr(const A* _a) { return _a[0]*_a[0]+_a[1]*_a[1]+_a[2]*_a[2]; }

    template <typename A>
    static A norm(const A* _a) {
      return sqrt(_a[0]*_a[0]+_a[1]*_a[1]+_a[2]*_a[2]);
    }
    template <typename A>
    static A norm1(const A* _a) {
      return fabs(_a[0])+fabs(_a[1])+fabs(_a[2]);
    }
    template <typename A>
    static A norminf(const A* _a) {
      return std::max(std::max(fabs(_a[0]),fabs(_a[1])),fabs(_a[2]));
    }
    static double mean(const double* _a) { return (_a[0]+_a[1]+_a[2])/3.0; }
    static float mean(const float* _a) { return (_a[0]+_a[1]+_a[2])/3.0f; }

    template <typename T>
    static void swap(T* _a,T* _b) {
      /*register*/ T c;
      c=_a[0]; _a[0]=_b[0]; _b[0]=c;
      c=_a[1]; _a[1]=_b[1]; _b[1]=c;
      c=_a[2]; _a[2]=_b[2]; _b[2]=c;
    }
    template <typename Op,typename A,typename B>
    static void map(Op _op,A* _a, const B* _b) {
      _a[0]=_op(_b[0]); _a[1]=_op(_b[1]); _a[2]=_op(_b[2]);
    }
    template <typename Op,typename A,typename B,typename C>
    static void map(Op _op,C* _c,const A* _a,const B* _b) {
      _c[0]=_op(_a[0],_b[0]);
      _c[1]=_op(_a[1],_b[1]);
      _c[2]=_op(_a[2],_b[2]);
    }
    template <typename T>
    static unsigned min(const T* _v) {
      if (_v[0]<=_v[1])
        return _v[0]<=_v[2] ? 0 : 2;
      else
        return _v[1]<=_v[2] ? 1 : 2;
    }
    template <typename T>
    static unsigned max(const T* _v) {
      if (_v[0]>=_v[1])
        return _v[0]>=_v[2] ? 0 : 2;
      else
        return _v[1]>=_v[2] ? 1 : 2;
    }
    template <typename T>
    static void minmax(unsigned* _minmax_idx,const T* _v) {
      if (_v[0]<=_v[1]) {
        if (_v[0]<=_v[2]) { _minmax_idx[0]=0; _minmax_idx[1]=(_v[1]>=_v[2] ? 1 : 2); }
        else              { _minmax_idx[0]=2; _minmax_idx[1]=(_v[0]>=_v[1] ? 0 : 1); }
      }
      else {
        if (_v[1]<=_v[2]) { _minmax_idx[0]=1; _minmax_idx[1]=(_v[0]>=_v[2] ? 0 : 2); }
        else              { _minmax_idx[0]=2; _minmax_idx[1]=(_v[0]>=_v[1] ? 0 : 1); }
      }
    }
    template <typename A,typename B,typename C>
    static void min(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]<=_b[0]) ? _a[0] : _b[0];
      _c[1]=(_a[1]<=_b[1]) ? _a[1] : _b[1];
      _c[2]=(_a[2]<=_b[2]) ? _a[2] : _b[2];
    }
    template <typename A,typename B,typename C>
    static void max(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]>=_b[0]) ? _a[0] : _b[0];
      _c[1]=(_a[1]>=_b[1]) ? _a[1] : _b[1];
      _c[2]=(_a[2]>=_b[2]) ? _a[2] : _b[2];
    }

    template <typename A,typename B,typename C>
    static void select(C* _c,const A* _a,const B* _b) {
      _c[0]=_a[_b[0]]; _c[1]=_a[_b[1]]; _c[2]=_a[_b[2]];
    }
    template <typename A,typename B,typename C>
    static void assign(C* _c,const A* _a,const B* _b) {
      _c[_b[0]]=_a[0]; _c[_b[1]]=_a[1]; _c[_b[2]]=_a[2];
    }

    template <typename A,typename B>
    static void isort(B* _b,const A* _a) {
      _b[0]=B(0); _b[1]=B(1); _b[2]=B(2);
      if (_a[_b[0]]>_a[_b[1]]) std::swap(_b[0],_b[1]);
      if (_a[_b[0]]>_a[_b[2]]) std::swap(_b[0],_b[2]);
      if (_a[_b[1]]>_a[_b[2]]) std::swap(_b[1],_b[2]);
    }
  };
# undef FixedAryOps_abcv
# undef FixedAryOps_abcs

  // N=4
# define FixedAryOps_abcv(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B* _b) { \
  _c[0]=_a[0] op _b[0]; _c[1]=_a[1] op _b[1]; \
  _c[2]=_a[2] op _b[2]; _c[3]=_a[3] op _b[3]; }
# define FixedAryOps_abcs(x,op) \
template <typename A,typename B,typename C> \
  static void x(C* _c,const A* _a,const B& _b) { \
  _c[0]=_a[0] op _b; _c[1]=_a[1] op _b; \
  _c[2]=_a[2] op _b; _c[3]=_a[3] op _b; }

  template <>
  struct FixedAryOps<4> {
    FixedAryOps_abcv(addv,+)
    FixedAryOps_abcv(subv,-)
    FixedAryOps_abcv(mulv,*)
    FixedAryOps_abcv(divv,/)

    FixedAryOps_abcs(adds,+)
    FixedAryOps_abcs(subs,-)
    FixedAryOps_abcs(muls,*)
    FixedAryOps_abcs(divs,/)

    template <typename A,typename B,typename C>
    static void subs2(C* _c,const A& _a,const B* _b) {
      _c[0]=_a-_b[0]; _c[1]=_a-_b[1]; _c[2]=_a-_b[2];  _c[3]=_a-_b[3];
    }
    template <typename A,typename B,typename C>
    static void divs2(C* _c,const A& _a,const B* _b) {
      _c[0]=_a/_b[0]; _c[1]=_a/_b[1]; _c[2]=_a/_b[2]; _c[3]=_a/_b[3];
    }
    template <typename A,typename B,typename C>
    static void lincomb(C* _c,
                               const A& _sa,const A* _a,
                               const B& _sb,const B* _b) {
      _c[0]=_a[0]*_sa+_b[0]*_sb;
      _c[1]=_a[1]*_sa+_b[1]*_sb;
      _c[2]=_a[2]*_sa+_b[2]*_sb;
      _c[3]=_a[3]*_sa+_b[3]*_sb;
    }
    template <typename T1,typename T2,typename T3,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3;
      _c[1]=_v1[1]*_s1+_v2[1]*_s2+_v3[1]*_s3;
      _c[2]=_v1[2]*_s1+_v2[2]*_s2+_v3[2]*_s3;
      _c[3]=_v1[3]*_s1+_v2[3]*_s2+_v3[3]*_s3;
    }
    template <typename T1,typename T2,typename T3,typename T4,typename C>
    static void lincomb(C* _c,
                        const T1& _s1,const T1* _v1,
                        const T2& _s2,const T2* _v2,
                        const T3& _s3,const T3* _v3,
                        const T4& _s4,const T4* _v4) {
      _c[0]=_v1[0]*_s1+_v2[0]*_s2+_v3[0]*_s3+_v4[0]*_s4;
      _c[1]=_v1[1]*_s1+_v2[1]*_s2+_v3[1]*_s3+_v4[1]*_s4;
      _c[2]=_v1[2]*_s1+_v2[2]*_s2+_v3[2]*_s3+_v4[2]*_s4;
      _c[3]=_v1[3]*_s1+_v2[3]*_s2+_v3[3]*_s3+_v4[3]*_s4;
    }
    template <typename A>
    static A sum(const A* _a) { return _a[0]+_a[1]+_a[2]+_a[3]; }

    template <typename A,typename B>
    static typename Result<A,B>::type dot(const A* _a,const B* _b) {
      typename Result<A,B>::type x=
        _a[0]*_b[0]+_a[1]*_b[1]+_a[2]*_b[2]+_a[3]*_b[3];
      return x;
    }

    template <typename A>
    static A sqr(const A* _a) {
      return _a[0]*_a[0]+_a[1]*_a[1]+_a[2]*_a[2]+_a[3]*_a[3];
    }
    template <typename A>
    static A norm(const A* _a) {
      return sqrt(_a[0]*_a[0]+_a[1]*_a[1]+_a[2]*_a[2]+_a[3]*_a[3]);
    }
    template <typename A>
    static A norm1(const A* _a) {
      return fabs(_a[0])+fabs(_a[1])+fabs(_a[2])+fabs(_a[3]);
    }
    template <typename A>
    static A norminf(const A* _a) {
      return std::max(std::max(fabs(_a[0]),fabs(_a[1])),std::max(fabs(_a[2]),fabs(_a[3])));
    }
    static double mean(const double* _a) { return (_a[0]+_a[1]+_a[2]+_a[3])*0.25; }
    static float mean(const float* _a) { return (_a[0]+_a[1]+_a[2]+_a[3])*0.25f; }

    template <typename T>
    static void swap(T* _a,T* _b) {
      /*register*/ T c;
      c=_a[0]; _a[0]=_b[0]; _b[0]=c;
      c=_a[1]; _a[1]=_b[1]; _b[1]=c;
      c=_a[2]; _a[2]=_b[2]; _b[2]=c;
      c=_a[3]; _a[3]=_b[3]; _b[3]=c;
    }
    template <typename Op,typename A,typename B>
    static void map(Op _op,A* _a, const B* _b) {
      _a[0]=_op(_b[0]); _a[1]=_op(_b[1]); _a[2]=_op(_b[2]); _a[3]=_op(_b[3]);
    }
    template <typename Op,typename A,typename B,typename C>
    static void map(Op _op,C* _c,const A* _a,const B* _b) {
      _c[0]=_op(_a[0],_b[0]);
      _c[1]=_op(_a[1],_b[1]);
      _c[2]=_op(_a[2],_b[2]);
      _c[3]=_op(_a[3],_b[3]);
    }
    template <typename T>
    static unsigned min(const T* _v) {
      unsigned i=FixedAryOps<2>::min(_v), j=FixedAryOps<2>::min(_v+2);
      return _v[i]<=_v[j+2] ? i : j+2;
    }
    template <typename T>
    static unsigned max(const T* _v) {
      unsigned i=FixedAryOps<2>::max(_v), j=FixedAryOps<2>::max(_v+2);
      return _v[i]>=_v[j+2] ? i : j+2;
    }
    template <typename T>
    static void minmax(unsigned* _minmax_idx,const T* _v) {
      ::VC::math::minmax(_minmax_idx,_v,4); // don't specialize
    }
    template <typename A,typename B,typename C>
    static void min(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]<=_b[0]) ? _a[0] : _b[0];
      _c[1]=(_a[1]<=_b[1]) ? _a[1] : _b[1];
      _c[2]=(_a[2]<=_b[2]) ? _a[2] : _b[2];
      _c[3]=(_a[3]<=_b[3]) ? _a[3] : _b[3];
    }
    template <typename A,typename B,typename C>
    static void max(C* _c,const A* _a,const B* _b) {
      _c[0]=(_a[0]>=_b[0]) ? _a[0] : _b[0];
      _c[1]=(_a[1]>=_b[1]) ? _a[1] : _b[1];
      _c[2]=(_a[2]>=_b[2]) ? _a[2] : _b[2];
      _c[3]=(_a[3]>=_b[3]) ? _a[3] : _b[3];
    }

    template <typename A,typename B,typename C>
    static void select(C* _c,const A* _a,const B* _b) {
      _c[0]=_a[_b[0]]; _c[1]=_a[_b[1]]; _c[2]=_a[_b[2]]; _c[3]=_a[_b[3]];
    }
    template <typename A,typename B,typename C>
    static void assign(C* _c,const A* _a,const B* _b) {
      _c[_b[0]]=_a[0]; _c[_b[1]]=_a[1]; _c[_b[2]]=_a[2]; _c[_b[3]]=_a[3];
     }

    template <typename A,typename B>
    static void isort(B* _b,const A* _a) {
      // 3-sort
      unsigned c[3]={0,1,2};
      if (_a[c[0]]>_a[c[1]]) std::swap(c[0],c[1]);
      if (_a[c[0]]>_a[c[2]]) std::swap(c[0],c[2]);
      if (_a[c[1]]>_a[c[2]]) std::swap(c[1],c[2]);

      // merge 4th in
      if (_a[c[0]]>=_a[3]) {
        _b[0]=B(3); _b[1]=B(c[0]); _b[2]=B(c[1]); _b[3]=B(c[2]);
      }
      else if (_a[c[1]]>=_a[3]) {
        _b[0]=B(c[0]); _b[1]=B(3); _b[2]=B(c[1]); _b[3]=B(c[2]);
      }
      else  if (_a[c[2]]>=_a[3]) {
        _b[0]=B(c[0]); _b[1]=B(c[1]); _b[2]=B(3); _b[3]=B(c[2]);
      }
      else {
        _b[0]=B(c[0]); _b[1]=B(c[1]); _b[2]=B(c[2]); _b[3]=B(3);
      }
    }
  };
# undef FixedAryOps_abcv
# undef FixedAryOps_abcs

# endif // DOXYGEN_SKIP

  /// @} core

# undef _LOOP0
# undef _LOOP1
# undef _LOOP1COND
# undef _LOOP1COND2

} // namespace math
} // namespace VC

#endif // __VC_MATH_ARRAY_HH
