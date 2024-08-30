#pragma once
//--------------------------------------------------------------------------//
#include <vc/vecn.hh>
//--------------------------------------------------------------------------//
namespace VC::vecn
{
  //--------------------------------------------------------------------------//
  template <int D, typename T, bool ENABLE_SIMD = true>
  vecn<D, T, ENABLE_SIMD> normalized(const vecn<D, T, ENABLE_SIMD> &_a)
  {
    return _a / norm(_a);
  }
  //--------------------------------------------------------------------------//
  template <int D, typename T, bool ENABLE_SIMD = true>
  T angle_u(const vecn<D, T, ENABLE_SIMD> &_a, const vecn<D, T, ENABLE_SIMD> &_b)
  {
    double d(_a | _b);
    assert(fabs(d) < 1.0 + std::numeric_limits<double>::epsilon() * 1000.0);
    if (d < -1.0)
      d = -1.0;
    else if (d > +1.0)
      d = +1.0;
    return acos(d);
  }
  //--------------------------------------------------------------------------//
  template <int D, typename T, bool ENABLE_SIMD = true>
  T angle(const vecn<D, T, ENABLE_SIMD> &_a, const vecn<D, T, ENABLE_SIMD> &_b)
  {
    return angle_u(normalized(_a), normalized(_b));
  }
  //--------------------------------------------------------------------------//
  template <int D, typename T, bool ENABLE_SIMD = true>
  bool same(const vecn<D, T, ENABLE_SIMD> &_a, const vecn<D, T, ENABLE_SIMD> &_b)
  {
    for (int i = 0; i < D; ++i)
      if (_a[i] != _b[i])
        return false;
    return true;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//