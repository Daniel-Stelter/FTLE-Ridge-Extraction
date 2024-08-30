//--------------------------------------------------------------------------//
#include "particlesystem.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  real energyProfileDerivative(real r)
  {
    constexpr real w = Globals::ENERGY_W;
    constexpr real d = Globals::ENERGY_D;
    if (r < w)
    {
      // (3*(d-1) / w) - (6*(d-1) * r/w^2) + (3*(d-1) * r^2/w^3)
      // = (3*(d-1) / w) * (1 - 2r/w + r^2/w^2)
      real r_w = r / w;
      return 3 * (d - 1) / w * (1 - (2 * r_w) + (r_w * r_w));
    }
    else if (r < 1)
    {
      // - (6*d * (r-w)^2 / (w-1)^3) - (6*d * (r-w) / (w-1)^2)
      // = (-6*d * (r-w) / (w-1)^2) * ((r-w)/(w-1) + 1)
      real rw_w1 = (r - w) / (w - 1);
      return (-6 * d * rw_w1 / (w - 1)) * (rw_w1 + 1);
    }
    else
      return 0;
  }
  //--------------------------------------------------------------------------//
  real energyProfile(real r)
  {
    constexpr real w = Globals::ENERGY_W;
    constexpr real d = Globals::ENERGY_D;
    if (r < w)
    {
      // 1 + (3 * (d-1) * r/w) - (3 * (d-1) * r^2/w^2) + ((d-1) * r^3/w^3)
      real r_w = r / w;
      real temp = (d - 1) * r_w;
      real result = 1 + 3 * temp;
      temp *= r_w;
      result -= 3 * temp;
      temp *= r_w;
      return result + temp;
    }
    else if (r < 1)
    {
      // d - (3*d * (r-w)^2 / (w-1)^2) - (2*d * (r-w)^3 / (w-1)^3)
      real rw_w1 = (r - w) / (w - 1);
      real temp = rw_w1 * rw_w1;
      real result = d - 3 * d * temp;
      temp *= rw_w1;
      return result - 2 * d * temp;
    }
    else
      return 0;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//