#ifndef VC_VECN_VECN_HH
#define VC_VECN_VECN_HH
//=============================================================================
#include <cstdint>
#include <cassert>
#include <cmath>

#include <iostream>
#include <type_traits>
#include <array>
#include <tuple>
#include <utility>
#include <algorithm>
#include <initializer_list>
#include <functional> // std::hash

#include "vecn_simd.hh"
//=============================================================================

/** \defgroup vc_vecn Linear algebra: small, fixed-size vector
    @{
 */

namespace VC {
/// \brief Linear algebra: small, fixed-size vector \ingroup vc_vecn
namespace vecn {
//=============================================================================

/** \brief A small, fixed-size vector.
    \ingroup vc_vecn

    \tparam D dimension
    \tparam T element type
    \tparam ENABLE_SIMD use SIMD data types *if available*

    The following operations are supported, most are defined as
    external functions.


    ### construction and access

    \code

    T         α,β;           // scalars
    vecn<D,T> x;             // uninitialized vector
    vecn<D,T> y(α);          // (α,α,...)
    vecn<D,T> z(α,β);        // (α,β,β,...)
    // ...
    x=vecn<D,T>(&α);         // initialize from array/pointer

    std::array<T,D> ary = ...;
    x=vecn<D,T>(ary);        // initialize from std::array<>

    auto x=make_vec(x,...);  // make vector, use common_type
    auto x=make_vec(ary);    // make vector from std::array<>

    x.size();                // == std::size_t(D)
    vecn<D,T>::value_type γ; // T
    auto iter=x.begin();     // similarly end()

    α=x[i];                  // get/set element

    y=x.let(i,α);            // copy with i-th element changed to α
    y=x.set(i,α);            // modify i-th element and return reference to self

    auto& ary=x.to_a();      // "cast" to std::array<T,D>&

    \endcode

    ### type traits

    is_vecn<vecn<D,T>>::value // true

    ### conversion and promotion

    \code

    vecn<D,S> xs=to<S>(x);   // convert type
    vecn<D,S> xs=x.to<S>();  // (may need an extra "template", see below)
    vecn<D,S> xs=x.template to<S>();

                             // same but deduce type
    vecn<D,S> xs=x.convert_to(S(0));
    vecn<D,S> ys=x.convert_to(xs);

    auto pr=promote(x,y);    // promote x and y to common type
    auto xp=pr.first;        // (xp and yp have same type)
    auto yp=pr.second;
                             // equivalent
    auto xp=promote_1st(x,y);
    auto yp=promote_2nd(x,y);

    \endcode


    ###  arithmetic

    \code

    x+=y;                    // same as x=x+y and similarly for -=, *=, /=, %=

    z=+x; z=-x;

    z=x+y;
    z=x-y;
    z=x*y;                   // element-wise (x[0]*y[0],x[1]*y[1],...)
    z=x/y;                   // element-wise (x[0]/y[0],x[1]/y[1],...)

    z=x%y;                   // element-wise modulo, for integral T only

    z=x+α; z=α+x;            // add constant vector (x[0]+α,x[1]+α,...)
    z=x-α; z=α-x;
    z=x*α; z=α*x;            // scale
    z=x/α;
    z=α/x;                   // element-wise reciprocal (α/x[0],α/x[1],...)

    z=cross(x,y);            // cross-product for D=3 and T=float|double
    z=x%y;                   // cross-product for D=3 and T=float|double

    z=axpy(α,x,y);           // z=α*x+y
    w=axpy(z,x,y);           // w=z*x+y

    \endcode


    ### dot product and norms

    \code

    α=dot(x,y);              // dot product
    α=(x|y);                 // dot product


    α=sqr(x);                // (x|x)
    α=norm(x);               // 2-norm sqrt(sqr(x)), synonym for norm2()
    α=norm2(x);              // 2-norm
    α=norm1(x);              // 1-norm
    α=norminf(x);            // maximum norm

    \endcode


    ### apply predicates and functions, filter values

    \code

    bool b=isfinitenorm(x);  // Are all entries finite?
    bool b=all(x,λ);         // Does predicate λ(x[i]) hold for all i=0,...,D-1?
    bool b=any(x,λ);         // Does predicate λ(x[i]) hold for any i=0,...,D-1?
    bool b=all(x,y,λ);       // Does predicate λ(x[i],y[i]) hold for all i?
    bool b=any(x,y,λ);       // Does predicate λ(x[i],y[i]) hold for any i?
    bool b=equal(x,y);       // Does x[i]==y[i] hold for all i?
    bool b=allequal(x,α);    // Does x[i]==α hold for all i?

                             // predefined predicates (lambdas):
    bool b=all(x,predicates::positive);
    bool b=all(x,predicates::gt(0)); // same: "greater than 0"
          // similar: negative, lt(x) ("less than x"),
          //          ge(x) "greater or equal x", le(x) "less or equal"


    bool b=in_box(a,x,b);    // Does a[i]<=x[i] && x[i]<=b[i] hold for all i?

    α=sum(x);  α=sum_of(x);   // sum of entries
    α=prod(x); α=prod_of(x);  // product of entries

    α=minimum(x);            // minimal entry
    β=maximum(x);            // maximal entry
    vecn<2,T> αβ=extrema(x); // extrema as 2-vector
    int i=indmin(x);         // get index of minimum
    int j=indmax(x);         // get index of maximum

    vec<D,int> xi=...;
    z=shuffle(x,xi);         // select entries x[xi[i]], xi has integral type
    xi=iota<D>(n0);          // generate vec<D,int>{n0,n0+1,...}
    auto q=invperm(p);       // generate inverse permuation s.th. p[q[i]]==i

    z=min(x,y);              // (min(x[0],y[0],),min(x[1],y[1],...)
    z=max(x,y);              // (max(x[0],y[0],),max(x[1],y[1],...)

    y=reverse(x);            // return (...,x[1],x[0]) with elements revered

    \endcode

    Let λ and ⊕ be functors such that β=λ(α) or γ=λ(α,β) and
    γ=⊕(α,β), which is below writtem  as binary infix operator α ⊕ β.

    \code

    y=map(x,λ);              // get (λ(x[0]),λ(x[1]),...)
    z=zip(x,y,λ);            // get (λ(x[0],y[0]),λ(x[1],y[1]),...)
    α=reduce(x,α₀,λ);        // get λ(...(λ(λ(α₀,x[0]),x[1]))...)
    α=reduce(x,λ);           // reduce(x,T(0),f)
    α=reduce(x,λ,⊕,y);       // get λ(x[0],y[0]) ⊕ λ(x[1],y[1]) ⊕ ... , where
                             // ⊕(a,b)=a⊕b defines a binary operator

    y=vmap(x,λ);             // same as above but use SIMD instructions ...
    z=vzip(x,y,λ);           // ... if available: apply λ on vector types
    α=vreduce(x,λ,⊕,y)

    x=mapind<D>(λ);          // map indices 0,1,... to (λ(0),λ(1),...)

    \endcode

    The functions map(), vmap(), zip(), vzip() and reduce(), vreduce()
    provide the basis for the implementation of many operators and
    functions. For instance, `x+y` is defined as `vzip(x,y,[](a,b){
    return a+b; }`.


    ### split and concatenate, "count"

    \code

    α=head(x);               // α==x[0]
    α=last(x);               // α==x[D-1]
    xr=tail(x);              // vector (x[1],x[2],...) of dimension D-1
    concat(x,α);             // concatenate as Matlab's vcat
    concat(α,x);
    concat(x,y);

    y=drop<K>(x);            // drop dimension K of x and return D-1 vector y

    x=inc(x,bounds);         // "count" x up to bounds for integral T

    α=get<I>(x);             // Returns X[I] for constexpr I.
    α=std::get<I>(x);

    auto r=
    call_unpacked(f,x,...);  // Calls r=f(x[0],x[1],...,x[D-1],...).
                             //  Note: A more general solution exists in
                             //  VC::utility::splat::invoke_splat!

    \endcode

    ### stream I/O

    \code

    cin >> x;
    cout << x;
    cout << pp(x);           // print as Matlab/Julia column vector

    \endcode

    ### Promotion

    `vecn` defines several rules to convert automatically between
     (potential) SIMD and non-SIMD types (vecn::is_simd()).

    `vecn` **doesn't** provide a default or implicit **type
    promotion**, e.g., the sum

    \code
    auto x=vec<2,int>(1,2)+vec<2,float>(4,5);
    \endcode

    is *undefined* unless functions in namespace
    VC::vecn::auto_promote are used.

    You require an explicit

    \code
    using namespace auto_promote
    \endcode

    ("promote all") or, e.g.,

    \code
    using auto_promote:operator+
    auto x=vec<2,int>(1,2)+vec<2,float>(4,5); // x is of type vec<2,float>
    \endcode

    This requires some extra effort but may help to make you aware of
    conversions and to keep code more "type safe".

    ### std

    \code
    using V = vecn<T,D,S>;
    V x {};
    std::get<D>(x)            // == x[D]
    std::tuple_size<V>::value // == std::size_t(D)
    std::rank<V>::value       // == std::size_t(1)
    std::extent<V,I>::value   // == D for I==0, else 0
    \end

 */
template <int D,typename T,bool ENABLE_SIMD=true>
struct vecn {

  /** @name constructors, type and dimension.
      @{
  */

  /// create uninitialized vector
  vecn() {}

  /// create vector from `_list`
  vecn(std::initializer_list<T> _list) {
    assert(_list.size()<=D);
    std::copy(_list.begin(),_list.end(),begin());
    const int n=int(_list.size());
    if (n>0) {
      for (int i=n;i<D;++i)
        data[i]=data[n-1];
    }
  }

  /// create vector `(_x,_x,...)`
  explicit vecn(T _x) {
    for (int i=0;i<D;++i)
      data[i]=_x;
  }
  /// create vector `(_x0,_x1,_xfirst,...,xlast,xlast,...)`
  template <typename... R, std::size_t... I>
  vecn(T _x0, T _x1, R... _xr)
      : self_t {_x0, _x1, static_cast<T>(_xr)...} {
    // TODO: give a more helpful/less annoying error message w/o "no matching..."
    //static_assert(std::conjunction<std::is_convertible<R,T>...>::value, "type mismtach");
    static constexpr int N= int(2 + sizeof...(R));
    static_assert(N <= D, "dimension mismatch");
  }


  /// create vector `(_x[0],_x[1],...)`
  vecn(const T* _x) {
    for (int i=0;i<D;++i)
      data[i]=_x[i];
  }
  /// create vector from `std::array`
  vecn(const std::array<T,std::size_t(D)>& _a) : vecn(&_a[0]) {}

  /// SIMD type *if requested and available* else `T[D]`
  using simd_t = simd_vecn<D,T,ENABLE_SIMD>;

  /// create vector from SIMD type
  vecn(simd_t _x) : vdata(_x) {}

  /// type of this vector
  using self_t = vecn<D,T,ENABLE_SIMD>;
  /// same but don't request SIMD type
  using self_nosimd_t = vecn<D,T,false>;

  /// element type
  using value_type = T;
  /// dimension
  static constexpr std::size_t size() { return std::size_t(D); }
  /// use SIMD type and operations
  static constexpr bool is_simd() { return simd_t::is_simd(); }

  /// create uninitialized vector of same type
  self_t similar() const { return self_t{}; }


  /// @}

  /** @name STL range
      @{
  */
  T* begin() { return data.begin(); }
  T* end() { return data.begin()+D; }
  const T* begin() const { return data.begin(); }
  const T* end() const { return data.begin()+D; }

  /// @}

  /** @name Access entries, operators
      @{
  */

  T operator[](int _i) const { assert(0<=_i && _i<D); return data[_i]; }
  T& operator[](int _i) { assert(0<=_i && _i<D); return data[_i]; }

  template <typename X> auto operator+=(X _x) { return *this=*this+_x; }
  template <typename X> auto operator-=(X _x) { return *this=*this-_x; }
  template <typename X> auto operator*=(X _x) { return *this=*this*_x; }
  template <typename X> auto operator/=(X _x) { return *this=*this/_x; }
  template <typename X> auto operator%=(X _x) { return *this=*this%_x; }

  /** Return copy of vector with `_i`-th component replaced by `_x`.
      `auto y=x.let(k,v);` is short for `auto y=x; y[k]=v;`
   */
  auto let(int _i,T _x) const {
    auto y {*this};
    y[_i]=_x;
    return y;
  }
  /// Set `_i`-th component of to `_x` and return modified vector `*this`.
  auto set(int _i,T _x) {
    (*this)[_i]=_x;
    return *this;
  }

  /// get data as `std::array<T,D>`
  auto& to_a() { return data; }
  /// get data as `std::array<T,D>`
  const auto& to_a() const { return data; }

  /// convert to `D` tuple
  auto to_tuple() const {
    return to_tuple(std::make_integer_sequence<int,D>{});
  }

  /// @}

  /** @name Explicit conversion and promotion.
      @{
  */

  /** Convert to `vecn<D,Tnew>`.
      \tparam Tnew is a scalar type for entries of `vecn`

      This should be used as follows

      \code
      auto xnew=x.template to<Tnew>(); // Mind "template" keyword!
      \endcode

      Calling the "external" version may be more convenient: the
      following us equivalent (but requires qualification by namespace)

      \code
      auto xnew=to<Tnew>(x);
      \endcode
   */
  template <typename Tnew> auto to() const;


  /// same as `to<Tnew>()` but deduct type from argument
  template <typename Tnew> auto to(Tnew) const { return to<Tnew>(); }

  /// @}

  union {
    std::array<T,D> data; //!< data
    simd_t vdata;         //!< data as SIMD type *if requested and available*
  };

 private:
  template <int... I>
  auto to_tuple(std::integer_sequence<int,I...>) const {
    return std::make_tuple(data[I]...);
  }
};

//-----------------------------------------------------------------------------
/** @name VC::vecn::vecn Explicit conversion and promotion.
    \ingroup vc_vecn

    \sa VC::vecn::auto_promote
 */
//-----------------------------------------------------------------------------

/// utility for converting `vecn` see VC::vecn::to(), VC::vecn::convert_to()
template <typename Tnew>
struct convertor {
  template <int D,typename T,bool S>
  auto operator()(const vecn<D,T,S>& _x) const {
    // TODO: __builtin_convertvector for SIMD (clang; only few type combinations)
    vecn<D,Tnew,S> z;
    for (int i=0;i<D;++i)
      z.data[i]=static_cast<Tnew>(_x.data[i]);
    return z;
  }
};

/// convert vector type `vecn<...,T,...>` to `vecn<...,Tnew,...>`
template <typename Tnew>
constexpr convertor<Tnew> to = convertor<Tnew>();

template <int D,typename T,bool S>
template <typename Tnew>
auto vecn<D,T,S>::to() const { return VC::vecn::to<Tnew>(*this); }

/// convert vector type `vecn<...,T,...>` to `vecn<...,Tnew,...>`
template <int D,typename T,typename Tnew,bool S>
auto convert_to(const vecn<D,T,S>& _x,Tnew) { return to<Tnew>(_x); }

template <int D,typename T,typename Tnew,bool S>
auto convert_to(const vecn<D,T,S>& _x,const vecn<D,Tnew,S>&) {
  return to<Tnew>(_x);
}

/** Utility to convert between `vecn` and other "alien" vectors.

    Convert between different vector representations/data structures.
    Assumes that the other, "alien" vector type supports `operator[]`.
*/
template <typename Vec,typename AlienVec>
struct alien_convertor {
  auto operator()(const Vec& _x) const {
    AlienVec x;
    for (int i=0;i<int(_x.size());++i)
      x[i]=_x[i];
    return x;
  }
  auto operator()(const AlienVec& _x) const {
    Vec x;
    for (int i=0;i<int(x.size());++i)
      x[i]=_x[i];
    return x;
  }
};


/// Create a vector from `_xr` with `T=std::common_type_t<R...>`
template <typename ...R>
auto make_vec(R... _xr) {
  using T=std::common_type_t<R...>;
  return vecn<sizeof...(R),T>(static_cast<T>(std::forward<R>(_xr))...);
}
/// Create a vector from `_ary`
template <typename T, std::size_t D>
auto make_vec(const std::array<T,D>& _ary) {
  return vecn<int(D),T> { _ary };
}

// TODO: make_vec from tuple


/// Test if `V` is an instance of `vecn`.
template <typename V>
struct is_vecn : std::integral_constant<bool,false> {};

template <int D,typename T>
struct is_vecn<vecn<D,T>> : std::integral_constant<bool,true> {};

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name functional-style operations
    \ingroup vc_vecn
    @{
 */
//-----------------------------------------------------------------------------

/// return `(_f(_x[0],_y[0]),_f(_x[1],_y[1]),...)`
template <int D,typename T,bool S,typename F>
auto zip(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,F _f) {
  using R=decltype(_f(_x[0],_y[0]));
  vecn<D,R,S> z;
  for (int i=0;i<D;++i)
    z.data[i]=_f(_x[i],_y[i]);
  return z;
}
/// same as VC::vecn::zip() but use SIMD operations if possible
template <int D,typename T,bool S,typename F>
auto vzip(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,F _f) {
  return vecn<D,T,S>(_x.vdata.zip(_f,_y.vdata));
}
/// similar to VC::vecn::zip() on 3 vectors, use SIMD operations if possible
template <int D,typename T,bool S,typename F>
auto
vzip(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,const vecn<D,T,S>& _a,F _f) {
  return vecn<D,T,S>(_x.vdata.zip(_f,_y.vdata,_a.vdata));
}

/// return `(_f(_x[0]),_f(_x[1]),...)`
template <int D,typename T,bool S,typename F>
auto map(const vecn<D,T,S>& _x,F _f) {
  using R=decltype(_f(_x[0]));
  vecn<D,R,S> z;
  for (int i=0;i<D;++i)
    z.data[i]=_f(_x[i]);
  return z;
}
/// same as VC::vecn::map() but use SIMD operations if possible
template <int D,typename T,bool S,typename F>
auto vmap(const vecn<D,T,S>& _x,F _f) {
  return vecn<D,T,S>(_x.vdata.map(_f));
}

/// map indices 0,1,... to `(_f(0),_f(1),...)`
template <int D,typename F>
auto mapind(F _f) {
  using T=decltype(_f(0));
  vecn<D,T> x;
  for (int i=0;i<D;++i)
    x[i]=_f(i);
  return x;
}

/// reduce entries to scalar `_op(..._op(_op(_initial,_x[0]),_x[1]),...)`
template <int D,typename T,bool S,typename R,typename OP>
auto reduce(R _initial,const vecn<D,T,S>& _x,OP _op) {
  R r=_initial;
  for (int i=0;i<D;++i)
    r=_op(r,_x[i]);
  return r;
}
/// reduce to `_f(_x[0],_y[0])⊕_f(_x[1],_y[1])⊕...`, where `a⊕b≡_op(a,b)`
template <int D,typename T,bool S,typename F,typename OP>
auto reduce(const vecn<D,T,S>& _x,F _f,OP _op,const vecn<D,T,S>& _y) {
  auto r=_f(_x[0],_y[0]);
  for (int i=1;i<D;++i)
    r=_op(r,_f(_x[i],_y[i]));
  return r;
}
/// reduce to `_f(_x[0])⊕_f(_x[1])⊕...`, where `a⊕b≡_op(a,b)`
template <int D,typename T,bool S,typename F,typename OP>
auto reduce(const vecn<D,T,S>& _x,F _f,OP _op) {
  return reduce(_x,[_f](auto x,auto) { return _f(x); },_op,_x);
}
/// same as VC::vecn::reduce() but use SIMD operations if possible
template <int D,typename T,bool S,typename F,typename OP>
auto vreduce(const vecn<D,T,S>& _x,F _f,OP _op,const vecn<D,T,S>& _y) {
  return _x.vdata.reduce(_f,_op,_y.vdata);
}
/// same as VC::vecn::reduce() but use SIMD operations if possible
template <int D,typename T,bool S,typename F,typename OP>
auto vreduce(const vecn<D,T,S>& _x,F _f,OP _op) {
  return vreduce(_x,[_f](auto x,auto) { return _f(x); },_op,_x);
}

/// test if *all* entries `_x[i]` satisfy predicate `_f(_x[i])`
template <int D,typename T,bool S,typename F>
bool all(const vecn<D,T,S>& _x,F _f) {
  for (int i=0;i<D;++i)
    if (!_f(_x[i]))
      return false;
  return true;
}
/// test if *any* entry `_x[i]` satisfies predicate `_f(_x[i])`
template <int D,typename T,bool S,typename F>
bool any(const vecn<D,T,S>& _x,F _f) {
  return !all(_x,[_f](auto x) { return !_f(x); });
}

//-----------------------------------------------------------------------------
namespace predicates {
//-----------------------------------------------------------------------------

/// lambda that tests if `>=x`
inline constexpr auto gt =
    [](const auto& x) { return [x](auto const& y) { return y>x; }; };

/// lambda that tests if `<=x`
inline constexpr auto lt =
    [](const auto& x) { return [x](auto const& y) { return y<x; }; };

/// lambda that tests if `_>=x`
inline constexpr auto ge =
    [](const auto& x) { return [x](auto const& y) { return y>=x; }; };

/// lambda that tests if `_<=x`
inline constexpr auto le =
    [](const auto& x) { return [x](auto const& y) { return y<=x; }; };

/// lambda that tests if `_==x`
inline constexpr auto eq =
    [](const auto& x) { return [x](auto const& y) { return y<=x; }; };

/// lambda that tests if `_!=x`
inline constexpr auto neq =
    [](const auto& x) { return [x](auto const& y) { return y!=x; }; };

/// lambda that tests if `x==0`
inline constexpr auto zero = eq(0);

/// lambda that tests if `x==0`
inline constexpr auto nonzero = neq(0);

/// lambda that tests if `x>0`
inline constexpr auto positive = gt(0);

/// lambda that tests if `x<0`
inline constexpr auto negative = lt(0);

/// lambda that tests if `x>=0`
inline constexpr auto nonnegative = ge(0);

/// lambda that tests if `x<=0`
inline constexpr auto nonpositive = le(0);

//-----------------------------------------------------------------------------
} // namespace predicates
//-----------------------------------------------------------------------------

/// test if *all* pairs of entries satisfy predicate `_f(_x[i],_y[i])`
template <int D,typename T,bool S,typename F>
bool all(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,F _f) {
  for (int i=0;i<D;++i)
    if (!_f(_x[i],_y[i]))
      return false;
  return true;
}
/// test if *any* pair of entries satisfies predicate `_f(_x[i],_y[i])`
template <int D,typename T,bool S,typename F>
bool any(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,F _f) {
  return !all(_x,_y,[_f](auto x,auto y) { return !_f(x,y); });
}

/// short for `all(_x,_y,==)`
template <int D,typename T,bool S>
bool equal(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return all(_x,_y,[](auto x,auto y) { return x==y; });
}
/// short for `all(_x,_y,==)`
template <int D,typename T,bool S>
bool allequal(const vecn<D,T,S>& _x,T _y) {
  return all(_x,[_y](auto x) { return x==_y; });
}
template <int D,typename T,bool S>
bool
in_box(const vecn<D,T,S>& _a,const vecn<D,T,S>& _x,const vecn<D,T,S>& _b) {
  for (int i=0;i<D;++i)
    if (!(_a[i]<=_x[i] && _x[i]<=_b[i]))
      return false;
  return true;
}



/// scan indices and select new index `i` if `_f(_x[selected],_x[i])`
template <int D,typename T,bool S,typename F>
int indscan(const vecn<D,T,S>& _x,F _f) {
  int is=0;
  for (int i=1;i<D;++i)
    is=_f(_x[is],_x[i]) ? i : is;
  return is;
}

/// get sum of entries
template <int D,typename T,bool S>
T sum(const vecn<D,T,S>& _x) {
  return vreduce(_x,[](auto x){ return x; },[](auto s,auto x){ return s+x; });
}
/// get product of entries
template <int D,typename T,bool S>
T prod(const vecn<D,T,S>& _x) {
  return vreduce(_x,[](auto x){ return x; },[](auto p,auto x){ return p*x; });
}

//
// The following are defined for CONVENIENCE: In case someone names a
// variable `sum=sum(x)` you don't require to add the namespace
// `sum=VC::vecn::sum(x)` but use `sum=sum_of(x)`.
//

/// get sum of entries (*synonym* for `sum`)
template <int D,typename T,bool S>
T sum_of(const vecn<D,T,S>& _x) { return sum(_x); }

/// get product of entries (*synonym* for `prod`)
template <int D,typename T,bool S>
T prod_of(const vecn<D,T,S>& _x) { return prod(_x); }

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name VC::vecn::vecn algebra
    \ingroup vc_vecn
    @{

    Let `x`,`y` be vectors of type `vecn<D,T>`, and `a` is a scalar of
    type `T`.

    - unary `-x`,`+x`
    - sum `x+y`, `x-y`
    - sum `x+a ≡ x+vecn(a)` and similarly `a+x`, `x-a`, `a-x`
    - scale `x*a`, `a*x`, `x/a`, `a/x ≡ (a/x₀,a/x₁,...)`
    - element-wise scale `x*y ≡ (x₀*y₀,x₁*y₁,...)` (**not** a dot product),
      and similarly `x/y` (**not** a pseudo-inverse)

    Modulo for *integral* types

    - `x%y ≡ (x[0]%y[0],x[1]%y[1],...)`, and similarly
    - `x%a ≡ x%vecn(a)`
*/
//-----------------------------------------------------------------------------

//
// +x, -x
//

template <int D,typename T,bool S>
auto operator+(const vecn<D,T,S>& _x) { return _x; }

template <int D,typename T,bool S>
auto operator-(const vecn<D,T,S>& _x) {
  return vmap(_x,[](auto x) { return -x; });
}

//
// x+y, x-y
//

template <int D,typename T,bool S>
auto operator+(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return vzip(_x,_y,[](auto a,auto b) { return a+b; });
}

template <int D,typename T,bool S>
auto operator-(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return vzip(_x,_y,[](auto a,auto b) { return a-b; });
}

//
// x+a, a+x, x-a, a-x
//

template <int D,typename T,bool S,typename Ta>
auto operator+(const vecn<D,T,S>& _x,Ta _a) {
  auto a=_x.vdata.vbroadcast(T(_a));
  return vmap(_x,[a](auto x) { return x+a; });
}
template <int D,typename T,bool S,typename Ta>
auto operator+(Ta _a,const vecn<D,T,S>& _x) { return _x+_a; }

template <int D,typename T,bool S,typename Ta>
auto operator-(const vecn<D,T,S>& _x,Ta _a) { return _x+(-_a); }

template <int D,typename T,bool S,typename Ta>
auto operator-(Ta _a,const vecn<D,T,S>& _x) { return (-_x)+_a; }

//
// x*a, a*x
//

template <int D,typename T,bool S,typename Ta>
auto operator*(const vecn<D,T,S>& _x,Ta _a) {
  // TODO: should check that Ta is a scalar type! (similarly below...)
  auto a=_x.vdata.vbroadcast(T(_a));
  return vmap(_x,[a](auto x) { return x*a; });
}
template <int D,typename T,bool S,typename Ta>
auto operator*(Ta _a,const vecn<D,T,S>& _x) { return _x*_a; }

//
// x/a, a/x
//

template <int D,typename T,bool S,typename Ta>
auto operator/(const vecn<D,T,S>& _x,Ta _a) {
  auto a=_x.vdata.vbroadcast(T(_a));
  return vmap(_x,[a](auto x) { return x/a; });
}
template <int D,typename T,bool S,typename Ta>
auto operator/(Ta _a,const vecn<D,T,S>& _x) {
  auto a=_x.vdata.vbroadcast(T(_a));
  return vmap(_x,[a](auto x) { return a/x; });
}

//
// element-wise x*y, x/y
//

template <int D,typename T,bool S>
auto operator*(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return vzip(_x,_y,[](auto a,auto b) { return a*b; });
}
template <int D,typename T,bool S>
auto operator/(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return vzip(_x,_y,[](auto a,auto b) { return a/b; });
}

//
// element-wise x%y (only integral types)
//

template <int D,typename T,bool S>
auto operator%(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  static_assert(std::is_integral<T>(),"require integral type");
  return vzip(_x,_y,[](auto a,auto b) { return a%b; });
}
template <int D,typename T,bool S>
auto operator%(const vecn<D,T,S>& _x,T _a) {
  static_assert(std::is_integral<T>(),"require integral type");
  auto a=_x.vdata.vbroadcast(T(_a));
  return vmap(_x,[a](auto x) { return x%a; });
}

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name Euclidean vector space, norms, etc.
    \ingroup vc_vecn
    @{
*/
//-----------------------------------------------------------------------------

//
// z=a*x+y
//

/// get `_a*_x+y` possibly using SIMD and FMA
template <int D,typename T,typename R,bool S>
auto axpy(R _a,const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  auto a=_x.vdata.vbroadcast(T(_a));
  return vzip(_x,_y,[a](auto x,auto y) { return a*x+y; });
}
/// get `_a*_x+y` possibly using SIMD and FMA
template <int D,typename T,bool S>
auto axpy(const vecn<D,T,S>& _a,const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return vzip(_x,_y,_a,[](auto x,auto y,auto a) { return a*x+y; });
}

//
// dot(x,y), (x|y), norm(x)
//

/// get dot product
template <int D,typename T,bool S>
T dot(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  // T sum=T(0);
  // for (int i=0;i<D;++i)
  //   sum+=_x[i]*_y[i];
  // return sum;
  return vreduce(_x,
                 [](auto x,auto y) { return x*y; },
                 [](auto sum,auto xy) { return sum+xy; },
                 _y);
}
/// same as VC::vecn::dot(), **always** use parenthesis `(_x|_y)`
template <int D,typename T,bool S>
T operator|(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) { return dot(_x,_y); }
/// same as `dot(_x,_x)`
template <int D,typename T,bool S>
T sqr(const vecn<D,T,S>& _x) { return dot(_x,_x); }

/// get Euclidean norm, synonym for VC::vecn::norm2()
template <int D,typename T,bool S>
T norm(const vecn<D,T,S>& _x) {
  static_assert(std::is_floating_point<T>(),"require floating point");
  return sqrt(sqr(_x));
}
/// get 2-norm, same as VC::vecn::norm()
template <int D,typename T,bool S>
T norm2(const vecn<D,T,S>& _x) { return norm(_x); }
/// get 1-norm
template <int D,typename T,bool S>
T norm1(const vecn<D,T,S>& _x) {
  return reduce(_x, // no vreduce: clang would not support x<T(0)
                [](auto x){ return std::abs(x); },
                [](auto sum,auto x){ return sum+x; });
}
/// get infinity norm (maximum norm)
template <int D,typename T,bool S>
T norminf(const vecn<D,T,S>& _x) {
  return reduce(_x,
                [](auto x){return std::abs(x);},
                [](auto r,auto x) { return std::max(r,x); });
}

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/** @name cross product of 3-vectors
    \ingroup vc_vecn

    - For 3-vectors (`D=3`) cross product as `cross(x,y)`.
    - If, in addition, `T` is a floating point type (`T=double|float`)
      the following is equivalent `x%y` (read `%` as operator `×`).

    @{
*/
//-----------------------------------------------------------------------------

/// return cross product `_x × _y`
template <typename T,bool S>
auto cross(const vecn<3,T,S>& _x,const vecn<3,T,S>& _y) {
  return vecn<3,T,S>(_x[1]*_y[2]-_x[2]*_y[1],
                     _x[2]*_y[0]-_x[0]*_y[2],
                     _x[0]*_y[1]-_x[1]*_y[0]);
}

/// return cross product `_x × _y`
template <bool S>
auto operator%(const vecn<3,float,S>& _x,const vecn<3,float,S>& _y) {
  return cross(_x,_y);
}
/// return cross product `_x × _y`
template <bool S>
auto operator%(const vecn<3,double,S>& _x,const vecn<3,double,S>& _y) {
  return cross(_x,_y);
}

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name predicates and "vectorized" functions
    \ingroup vc_vecn
    @{
*/
//-----------------------------------------------------------------------------


//
// isfinite(x)
//

/// test if all entries are finite
template <int D,typename T,bool S>
auto isfinitenorm(const vecn<D,T,S>& _x) {
  return all(_x,[](auto x) { return std::isfinite(x); });
}

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name shuffling and sorting
    \ingroup vc_vecn
    @{
*/
//-----------------------------------------------------------------------------

//
// shuffle
//

/// get `(_x[_mask[0]],_x[_mask[1],...)`
template <int D,typename T,bool S,typename TI,bool SI>
auto shuffle(const vecn<D,T,S>& _x,const vecn<D,TI,SI>& _mask) {
  return vecn<D,T,S>(shuffle(_x.vdata,_mask.vdata));
}

/// generate `(n0,n0+1,...)`, e.g., identity permutation for shuffle()
template <int D>
auto iota(int n0=0) {
  vecn<D,int> xi {};
  for (int i=0;i<D;++i)
    xi[i]=i+n0;
  return xi;
}

/// generate inverse of permutation vector `p`
template <int D,bool S>
auto invperm(const vecn<D,int,S>& _permutation, int _idefault = 0) {
  vecn<D,int,S> ip { _idefault };
  for (int i=0;i<D;++i)
    ip[_permutation[i]]=i;
  return ip;
}

//
// extrema
//

/// get smallest entry
template <int D,typename T,bool S>
auto minimum(const vecn<D,T,S>& _x) {
  return reduce(_x,
                [](auto x){return x;},
                [](auto r,auto x) { return std::min(r,x); });
}
/// get largest entry
template <int D,typename T,bool S>
auto maximum(const vecn<D,T,S>& _x) {
  return reduce(_x,
                [](auto x){return x;},
                [](auto r,auto x) { return std::max(r,x); });
}

/// get minimum and maximum
template <int D,typename T,bool S>
auto extrema(const vecn<D,T,S>& _x) {
  T low=_x[0],high=_x[0];
  for (int i=1;i<D;++i) {
    low =std::min(low,_x[i]);
    high=std::max(low,_x[i]);
  }
  return vecn<2,T,S>(low,high);
}

/// get index of smallest entry
template <int D,typename T,bool S>
int indmin(const vecn<D,T,S>& _x) {
  return indscan(_x,[](auto xs,auto x) { return x<xs; });
}
/// get index of smallest entry
template <int D,typename T,bool S>
int indmax(const vecn<D,T,S>& _x) {
  return indscan(_x,[](auto xs,auto x) { return x>xs; });
}

/// get vector of minima
template <int D,typename T,bool S>
auto min(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return zip(_x,_y,[](auto x,auto y) { return std::min(x,y); }); // no vzip
}
/// get vector of maxima
template <int D,typename T,bool S>
auto max(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  return zip(_x,_y,[](auto x,auto y) { return std::max(x,y); });  // no vzip
}

template <int D,typename T,bool S>
auto reverse(const vecn<D,T,S>& _x) {
  vecn<D,T,S> r;
  for (int i=0;i<D;++i)
    r[i]=_x[D-i-1];
  return r;
}

// TODO: construct from reverse

// TODO: sort, psort

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name split and concatenate
    \ingroup vc_vecn
    @{

    Let `x=(x0,x1,x2,...,xend)` and `y=(y0,y1,...,yend)`. And let `a`
    be a scalar.

    ### split

    - `x0 == VC::vecn::head(x)`
    - `xend == VC::vecn::last(x)`
    - `(x1,x2,...,xend) == tail(x)`
    - `VC::vecn::first<N>(x)` gives vector of `N` first entries `(x0,x1,...)`
    - `VC::vecn::last<N>(x)` gives vector of `N` last entries `(,...,xend)`

    ### concatenate

    - `(x0,x1,...,xend,y1,y2,...,yend) == concat(x,y)`
    - `(a,x0,x1,...,xend) == concat(a,x)`
    - `(x0,x1,...,xend,a) == concat(x,a)`

    The namespace vc::vecn::concatenate provides concatenation with
    `operator,`. Then `(x,y)` is short for `concat(x,y)`.

    ### increment: count indices

    Let `i=(i0,i1,...)` be an vector of an integer type and same for
    `bounds=(b0,b1,...)`.

    - `inc(i,bounds) == (i0+1,i1,...)` for `i0+1<b0`, and
    - `inc(i,bounds) == concat(0,inc(tail(i),tail(bounds))` else

    Then VC::vecn::inc() would count up each dimension to `bounds`
    with a possible "overflow", i.e., `last(x)>last(bounds)` when the
    least dimension `D-1` is reached.
*/

//-----------------------------------------------------------------------------

//
// first and last
//

/// get first `N` entries `(_x[0],_x[1],...,_x[N-1])`
template <int N,int D,typename T,bool S>
auto first(const vecn<D,T,S>& _x) {
  static_assert(D>=N,"too short");
  vecn<N,T,S> y;
  for (int i=0;i<N;++i)
    y.data[i]=_x.data[i];
  return y;
}
/// get last `N` entries `(_x[D-N],_x[D-N+1],...,_x[D-1])`
template <int N,int D,typename T,bool S>
auto last(const vecn<D,T,S>& _x) {
  static_assert(D>=N,"too short");
  vecn<N,T,S> y;
  for (int i=0;i<N;++i)
    y.data[i]=_x.data[D-N+i];
  return y;
}

//
// head, last, tail, concat
//

/// get first entry `_x[0]`, use with VC::vecn::tail()
template <int D,typename T,bool S>
auto head(const vecn<D,T,S>& _x) { return _x[0]; }
/// get last entry `_x[_x.size()-1]`
template <int D,typename T,bool S>
auto last(const vecn<D,T,S>& _x) { return _x[D-1]; }

/// get vector with VC::vecn::head() removed
template <int D,typename T,bool S>
auto tail(const vecn<D,T,S>& _x) {
  static_assert(D>1,"tail must not be empty");
  vecn<D-1,T,S> y;
  for (int i=1;i<D;++i)
    y.data[i-1]=_x.data[i];
  return y;
}

/// get concatenation
template <int Dx,int Dy,typename T,bool Sx,bool Sy>
auto concat(const vecn<Dx,T,Sx>& _x,const vecn<Dy,T,Sy>& _y) {
  vecn<Dx+Dy,T,Sx || Sy> z;
  for (int i=0;i<Dx;++i)
    z.data[i]=_x.data[i];
  for (int i=0;i<Dy;++i)
    z.data[i+Dx]=_y.data[i];
  return z;
}
template <int D,typename T,bool S>
auto concat(const vecn<D,T,S>& _x,T _y) {
  vecn<D+1,T,S> z;
  for (int i=0;i<D;++i)
    z.data[i]=_x.data[i];
  z.data[D]=_y;
  return z;
}
template <int D,typename T,bool S>
auto concat(T _y,const vecn<D,T,S>& _x) {
  vecn<D+1,T,S> z;
  z.data[0]=_y;
  for (int i=0;i<D;++i)
    z.data[i+1]=_x.data[i];
  return z;
}

//
// drop dimension
//

/// drop dimension `K-1`
template <int D,int K>
struct drop_dimension_t {
  template <typename T,bool S>
  auto operator()(const vecn<D,T,S>& _x) const {
    static_assert(1<=K && K<=D,"invalid dimension (mind 1-based K)");
    return concat(first<K-1>(_x),last<D-K>(_x));
  }
};

// Note: K is in 1,...,D,
// because specialization drop_dimension_t<D,D-1> is not allowed!

template <int D>
struct drop_dimension_t<D,D> {
  template <typename T,bool S>
  auto operator()(const vecn<D,T,S>& _x) const { return first<D-1>(_x); }
};
template <int D>
struct drop_dimension_t<D,1> {
  template <typename T,bool S>
  auto operator()(const vecn<D,T,S>& _x) const { return last<D-1>(_x); }
};

template <int K,int D,typename T,bool S>
auto drop(const vecn<D,T,S>& _x) {
  return (drop_dimension_t<D,K+1>{})(_x);
}

template <int K,int D,typename T,bool S>
constexpr const T& get(const vecn<D,T,S>& _x) noexcept {
  static_assert(0<=K && K<D,"invalid index");
  return _x.begin()[K];
}
template <int K,int D,typename T,bool S>
constexpr const T&& get(const vecn<D,T,S>&& _x) noexcept {
  static_assert(0<=K && K<D,"invalid index");
  return _x.begin()[K];
}

template <int K,int D,typename T,bool S>
constexpr T& get(vecn<D,T,S>& _x) noexcept {
  static_assert(0<=K && K<D,"invalid index");
  return _x.begin()[K];
}
template <int K,int D,typename T,bool S>
constexpr T&& get(vecn<D,T,S>&& _x) noexcept {
  static_assert(0<=K && K<D,"invalid index");
  return _x.begin()[K];
}

// TODO: vecn to std::tuple, to std::array


template <typename F,
          int D,typename T,bool S,typename... Args,std::size_t... I>
auto _call_unpacked(F f,
                    const vecn<D,T,S>& _x,std::index_sequence<I...>,Args... args) {
  return f(std::forward<T>(_x[I])...,std::forward<Args>(args)...);
}

/** Returns `f(_x[0],_x[1],...,_x[D-1], args...)`.

    If `f` has additional arguments *before* `_x[0],...` (e.g.,
    `f(a,b,x0,x1,...,c,...)` invoke, for instance, as

    ~~~~
    call_unpacked([f,a,b](auto&& ...args) { f(a,b,args...); }, x);
    ~~~

    **Note:** VC::utility::splat defines a more general solution that
    enables splatting of multiple arguments and uses a nicer syntax:

    ~~~
    using VC::utility::splat::invoke_splat;
    using VC::utility::splat::__;

    invoke_splat(f,__[_x], args...);
    ~~~
 */

template <typename F,
          int D,typename T,bool S,typename... Args>
auto call_unpacked(F f,const vecn<D,T,S>& _x,Args... args) {
  return _call_unpacked(f,_x,
                        std::make_index_sequence<D>{}, std::forward<Args>(args)...);
}

//-----------------------------------------------------------------------------
/// \brief Concatenate `vecn` by `operator,`. \ingroup vc_vecn
namespace concatenate {
//-----------------------------------------------------------------------------

/// short for VC::vecn::concat(_x,_y)
template <int Dx,int Dy,typename T,bool Sx,bool Sy>
auto operator,(const vecn<Dx,T,Sx>& _x,const vecn<Dy,T,Sy>& _y) {
  return concat(_x,_y);
}
template <int D,typename T,bool S>
auto operator,(const vecn<D,T,S>& _x,T _y) { return concat(_x,_y); }

template <int D,typename T,bool S>
auto operator,(T _x,const vecn<D,T,S>& _y) { return concat(_x,_y); }

//-----------------------------------------------------------------------------
} // namespace concatenate
//-----------------------------------------------------------------------------

/// functor used by VC::vecn::inc()
template <int D,typename T>
struct inc_t {
  auto operator()(const vecn<D,T>& _i,const vecn<D,T>& _bounds) {
    static_assert(std::is_integral<T>(),"require integral type");

    if (_i[0]+1<_bounds[0]) {
      vecn<D,T > next=_i;
      next[0]=_i[0]+1;
      return next;
    }
    else {
      auto next=concat(T(0),inc_t<D-1,T>()(tail(_i),tail(_bounds)));
      return next;
    }
  }
};

# ifndef DOXYGEN_SKIP

template <typename T>
struct inc_t<1,T> {
  auto operator()(const vecn<1,T>& _i,const vecn<1,T>&) {
    static_assert(std::is_integral<T>(),"require integral type");
    return vecn<1,T>(_i[0]+1);
  }
};

# endif

/// count `_i` up by one and respect `_bounds` to "switch" dimensions
template <int D,typename T>
vecn<D,T> inc(const vecn<D,T>& _i,const vecn<D,T>& _bounds) {
  return inc_t<D,T>()(_i,_bounds);
}


//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name stream I/O
    \ingroup vc_vecn
    @{
*/
//-----------------------------------------------------------------------------

//
// out << x, in >> x
//

/// print `_x` as text to `_out
template <int D,typename T,bool S>
std::ostream& operator<<(std::ostream& _out,const vecn<D,T,S>& _x) {
  static_assert(D>0,"dimension");
  for (int i=0;i<D-1;++i)
    _out << _x[i] << ' ';
  return _out << _x[D-1];
}

/// pretty print tag generated by VC::vecn::pp()
template <int D,typename T,bool S>
struct pp_t {
  const vecn<D,T,S>& x;
};

/// tag `_x` for pretty printing as Matlab/Julia column vector
template <int D,typename T,bool S>
auto pp(const vecn<D,T,S>& _x) { return pp_t<D,T,S> { _x }; }

/// print `_x` as text to `_out
template <int D,typename T,bool S>
std::ostream& operator<<(std::ostream& _out,const pp_t<D,T,S>& _x) {
  static_assert(D>0,"dimension");
  _out << '[';
  for (int i=0;i<D-1;++i)
    _out << _x.x[i] << ';';
  return _out << _x.x[D-1] << ']';
}

/// get `_x` from text read from `_in`
template <int D,typename T,bool S>
std::istream& operator>>(std::istream& _in, vecn<D,T,S>& _x) {
  static_assert(D>0,"dimension");
  for (int i=0;i<D && _in;++i)
    _in >> _x[i];
  return _in;
}

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
/** @name Relational operators: lexicographic order
    \ingroup vc_vecn

    Equality testing `x == y` is only supported for integer types.

    @{
*/
//-----------------------------------------------------------------------------

template <int D,typename T,bool S>
constexpr bool operator<(const vecn<D,T,S>& _x, const vecn<D,T,S>& _y) {
  for (int i=0;i<D;++i) {
    if (_x[i] < _y[i])
      return true;
    else if (_x[i] > _y[i])
      return false;
  }
  return false;
}

template <int D,typename T,bool S>
constexpr bool operator<=(const vecn<D,T,S>& _x, const vecn<D,T,S>& _y) {
  for (int i=0;i<D;++i) {
    if (_x[i] < _y[i])
      return true;
    else if (_x[i] > _y[i])
      return false;
  }
  return true;
}

template <int D,typename T,bool S>
constexpr bool operator>(const vecn<D,T,S>& _x, const vecn<D,T,S>& _y) {
  for (int i=0;i<D;++i) {
    if (_x[i] > _y[i])
      return true;
    else if (_x[i] < _y[i])
      return false;
  }
  return false;
}

template <int D,typename T,bool S>
constexpr bool operator>=(const vecn<D,T,S>& _x, const vecn<D,T,S>& _y) {
  for (int i=0;i<D;++i) {
    if (_x[i] > _y[i])
      return true;
    else if (_x[i] < _y[i])
      return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

template <int D,typename T,bool S,
          typename std::enable_if_t<!std::is_floating_point<T>::value, int> = 0>
constexpr bool operator==(const vecn<D,T,S>& _x, const vecn<D,T,S>& _y) {
  for (int i=0;i<D;++i) {
    if (_x[i] != _y[i])
      return false;
  }
  return true;
}

template <int D,typename T,bool S,
          typename std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
constexpr bool operator==(const vecn<D,T,S>&, const vecn<D,T,S>&) {
  static_assert(!std::is_same<T,T>::value,
                "Equality test x == y s not supported for floating point types!");
  return false;
}

template <int D,typename T,bool S>
constexpr bool operator!=(const vecn<D,T,S>& _x, const vecn<D,T,S>& _y) {
  return !(_x == _y);
}


// TODO: is_floating_point -- error message

//-----------------------------------------------------------------------------
/// @}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

// TODO: rand (?)

// TODO: scale(2) --> "factor" axpy
// TODO: specialize small dimensions

// TODO: geometry (clamp dot, angle, ...)

//-----------------------------------------------------------------------------

//=============================================================================
} //namespace vecn
} //namespace VC
//=============================================================================
// # ifndef DOXYGEN_SKIP
namespace std {
//-----------------------------------------------------------------------------
template <int D,typename T,bool S>
struct tuple_size<VC::vecn::vecn<D,T,S>>
    : public std::integral_constant<std::size_t,D> {};

template <std::size_t I,int N,typename T,bool S>
constexpr const auto& get(const VC::vecn::vecn<N,T,S>& _x) noexcept {
  return VC::vecn::get<I>(_x);
}
template <std::size_t I,int N,typename T,bool S>
constexpr const auto&& get(const VC::vecn::vecn<N,T,S>&& _x) noexcept {
  return VC::vecn::get<I>(_x);
}
template <std::size_t I,int N,typename T,bool S>
constexpr auto& get(VC::vecn::vecn<N,T,S>& _x) noexcept {
  return VC::vecn::get<I>(_x);
}
template <std::size_t I,int N,typename T,bool S>
constexpr auto&& get(VC::vecn::vecn<N,T,S>&& _x) noexcept {
  return VC::vecn::get<I>(_x);
}

//-----------------------------------------------------------------------------

template<int D,typename T,bool S>
struct rank<VC::vecn::vecn<D,T,S>>
    : public std::integral_constant<std::size_t, 1> {};

template<int D,typename T,bool S>
struct extent<VC::vecn::vecn<D,T,S>, 0>
    : public std::integral_constant<std::size_t, D> {};

//-----------------------------------------------------------------------------

template<int D,typename T,bool S>
struct hash<VC::vecn::vecn<D,T,S>> {
  using argument_type = VC::vecn::vecn<D,T,S>;
  using result_type = std::size_t;
  result_type operator()(const argument_type& _x) {
    // https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
    // https://www.nullptr.me/2018/01/15/hashing-stdpair-and-stdtuple/
    constexpr std::hash<T> hasher{};
    return VC::vecn::reduce(std::size_t{0}, _x,
             [hasher](const result_type seed, const T& xi) {
               auto val = hasher(xi) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
               return seed ^ val;
             });
  }
};

template<int D,typename T,bool S>
struct less<VC::vecn::vecn<D,T,S>> {
  using argument_type = VC::vecn::vecn<D,T,S>;

  constexpr bool
  operator()(const argument_type& _x, const argument_type& _y) {
    return _x < _y;
  }
};

//-----------------------------------------------------------------------------
} // namespace std
// # endif
//=============================================================================
/// @}
# include "vecn_promotion.hh"
# include "vecn_specialization.hh"
//=============================================================================
#endif // VC_VECN_VECN_HH defined
