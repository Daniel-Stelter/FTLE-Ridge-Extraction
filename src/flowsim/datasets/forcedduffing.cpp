#include "forcedduffing.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  ForcedDuffing2D::ForcedDuffing2D()
      : ForcedDuffing2D({Range(-2, 2), Range(-1.5, 1.5)}) {}
  //--------------------------------------------------------------------------//
  ForcedDuffing2D::ForcedDuffing2D(const Domain<2> &space_domain_used)
      : Flow<2>({Range(), Range()},
                space_domain_used,
                Range(),
                "ForcedDuffing2D") {}
  //--------------------------------------------------------------------------//
  MaybeVec<2>
  ForcedDuffing2D::v(const VecR<2> &pos, real t) const
  {
    real x = pos[0];
    real y = pos[1];
    VecR<2> v;
    v[0] = y;
    v[1] = x - x * x * x + 0.1 * sin(t);
    return v;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//
