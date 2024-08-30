#ifndef VC_VECN_VECN_PROMOTION_HH
#define VC_VECN_VECN_PROMOTION_HH
//=============================================================================
namespace VC {
namespace vecn {
//=============================================================================
# ifndef DOXYGEN_SKIP
//-----------------------------------------------------------------------------

template <int D,typename T>
auto simd(const vecn<D,T,true>& _x) { return _x; }

template <int D,typename T>
auto simd(const vecn<D,T,false>& _x) {
  vecn<D,T,true> x; // Note that we don't have guarantees on alignment!
  for (size_t i=0;i<D;++i)
    x[i]=_x[i];
  return x;
}
// TODO: specialize

template <int D,typename T>
auto no_simd(const vecn<D,T,false>& _x) { return _x; }

template <int D,typename T>
auto no_simd(const vecn<D,T,true>& _x) { return vecn<D,T,false>(_x.begin()); }

#  define PROMOTE_SIMD_OP(op)                                       \
  template <int D,typename T,bool S>                              \
  auto operator op(const vecn<D,T,S>& _x,const vecn<D,T,!S> _y) { \
    return simd(_x) op simd(_y);                                  \
  }


PROMOTE_SIMD_OP(+)
PROMOTE_SIMD_OP(-)
PROMOTE_SIMD_OP(*)
PROMOTE_SIMD_OP(/)
PROMOTE_SIMD_OP(%)

PROMOTE_SIMD_OP(|)

#  undef PROMOTE_SIMD_OP

template <int D,typename T,bool S,typename F>
auto zip(const vecn<D,T,S>& _x,const vecn<D,T,!S>& _y,F _f) {
  return zip(no_simd(_x),no_simd(_y),_f);
}
template <int D,typename T,bool S,typename F>
auto vzip(const vecn<D,T,S>& _x,const vecn<D,T,!S>& _y,F _f) {
  return vzip(simd(_x),simd(_y),_f);
}

template <int D,typename T,bool S,typename F,typename OP>
auto reduce(const vecn<D,T,S>& _x,F _f,OP _op,const vecn<D,T,!S>& _y) {
  return reduce(no_simd(_x),_f,_op,no_simd(_y));
}
template <int D,typename T,bool S,typename F,typename OP>
auto vreduce(const vecn<D,T,S>& _x,F _f,OP _op,const vecn<D,T,!S>& _y) {
  return reduce(simd(_x),_f,_op,simd(_y));
}

template <int D,typename T,bool S,typename F>
auto dot(const vecn<D,T,S>& _x,const vecn<D,T,!S>& _y) {
  return dot(simd(_x),simd(_y));
}

// TODO: axpy (?)

//-----------------------------------------------------------------------------
# else // DOXYGEN_SKIP
//-----------------------------------------------------------------------------

/// promote `_x` to (potential) SIMD type
template <int D,typename T,bool SIMD> auto simd(const vecn<D,T,true>& _x);

// promote `_x` to non-SIMD type (e.g., to drop assumption on alignment)
template <int D,typename T,bool SIMD> auto no_simd(const vecn<D,T,false>& _x);

//-----------------------------------------------------------------------------
# endif // DOXYGEN_SKIP
//-----------------------------------------------------------------------------

/** \brief Automatic type promotion for VC::vecn::vecn.
    \ingroup vc_vecn

    The namespace defines automatic type promotion for

    - VC::vecn::operator+
    - VC::vecn::operator-
    - VC::vecn::operator*
    - VC::vecn::operator/
    - VC::vecn::operator%
    - VC::vecn::operator|
    - VC::vecn::zip
    - VC::vecn::vzip
    - VC::vecn::reduce (version with with two vector values parameters)
    - VC::vecn::vreduce (version with with two vector values parameters)
    - VC::vecn::dot()

    This means that `op(x,y)` will first convert `x` and `y` to a
    common type defined by VC::vecn::auto_promote::promote_type()
    before doing the computation.

    This enables, e.g.,

    \code
    vecn<2,int>   x(1,2);
    vecn<2,float> y(3,4);

    using auto_promote::operator+;

    auto z=x+y;   // type of z is vecn<2,float>
    \endcode
 */
namespace auto_promote {
//-----------------------------------------------------------------------------

/// get *undefined* value of common type, use return value with `decltype`
template <typename Tx,typename Ty>
auto promote_type(Tx _x,Ty _y) {
  static_assert(std::is_scalar<Tx>() && std::is_scalar<Ty>(),
                "predefined promotion rules only for scalars");
  return _x+_y;
}

/// get `std::pair<>` of `(_x,_y)` promoted to common promote_type()
template <int D,typename Tx,typename Ty,bool Sx,bool Sy>
auto promote(const vecn<D,Tx,Sx>& _x,const vecn<D,Ty,Sy>& _y) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return std::make_pair(C(_x),C(_y));
}

/// get first argument `_x` promoted to common promote_type() of `_x,_y`
template <int D,typename Tx,typename Ty,bool Sx,bool Sy>
auto promote_1st(const vecn<D,Tx,Sx>& _x,const vecn<D,Ty,Sy>& _y) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return C(_x);
}
/// get second argument `_x` promoted to common promote_type() of `_x,_y`
template <int D,typename Tx,typename Ty,bool Sx,bool Sy>
auto promote_2nd(const vecn<D,Tx,Sx>& _x,const vecn<D,Ty,Sy>& _y) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return C(_y);
}

//-----------------------------------------------------------------------------
# ifndef DOXYGEN_SKIP
//-----------------------------------------------------------------------------

# define AUTO_PROMOTE_RULE(op)                                          \
  template <int D,typename Tx,typename Ty,bool Sx,bool Sy>              \
  auto operator op(const vecn<D,Tx,Sx>& _x,const vecn<D,Ty,Sy>& _y) {   \
    auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();            \
    return C(_x) op C(_y);                                              \
  }                                                                     \

AUTO_PROMOTE_RULE(+)
AUTO_PROMOTE_RULE(-)
AUTO_PROMOTE_RULE(*)
AUTO_PROMOTE_RULE(/)
AUTO_PROMOTE_RULE(%)

AUTO_PROMOTE_RULE(|)

template <int D,typename T,bool S,typename F>
auto zip(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,F _f) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return zip(C(_x),C(_y),_f);
}
template <int D,typename T,bool S,typename F>
auto vzip(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y,F _f) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return vzip(C(_x),C(_y),_f);
}

template <int D,typename T,bool S,typename F,typename OP>
auto reduce(const vecn<D,T,S>& _x,F _f,OP _op,const vecn<D,T,S>& _y) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return reduce(C(_x),_f,_op,C(_y));
}
template <int D,typename T,bool S,typename F,typename OP>
auto vreduce(const vecn<D,T,S>& _x,F _f,OP _op,const vecn<D,T,S>& _y) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return reduce(C(_x),_f,_op,C(_y));
}

template <int D,typename T,bool S,typename F>
auto dot(const vecn<D,T,S>& _x,const vecn<D,T,S>& _y) {
  auto C=convertor<decltype(promote_type(_x[0],_y[0]))>();
  return dot(C(_x),C(_y));
}

// TODO: axpy (?)

//-----------------------------------------------------------------------------
# endif // DOXYGEN_SKIP
//-----------------------------------------------------------------------------
} // namespace auto_promote
//-----------------------------------------------------------------------------

//=============================================================================
} //namespace vecn
} //namespace VC
//=============================================================================
#endif // VC_VECN_VECN_PROMOTION_HH defined
