#include "tornado.hpp"
//--------------------------------------------------------------------------//
#include "globals.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  Tornado3D::Tornado3D()
      : Tornado3D({Range(-2.5, 2.5), Range(-2.5, 2.5), Range(-5, 5)}) {}
  //--------------------------------------------------------------------------//
  Tornado3D::Tornado3D(const Domain<3> &space_domain_used)
      : Flow<3>({Range(), Range(), Range()},
                space_domain_used,
                Range(),
                "Tornado3D") {}
  //--------------------------------------------------------------------------//
  MaybeVec<3>
  Tornado3D::v(const VecR<3> &pos, real) const
  {
    const real &x = pos[0];
    const real &y = pos[1];
    const real &z = pos[2];
    VecR<3> v;

    const real xc = 0.5 + 0.1 * sin(10.0 * z);                 // For each z-slice, determine the spiral circle.
    const real yc = 0.5 + 0.1 * cos(3.0 * z);                  // (xc,yc) determine the center of the circle.
    const real r = 0.1 + 0.4 * z * z + 0.1 * z * sin(8.0 * z); // The radius also changes at each z-slice.
    const real r2 = 0.2 + 0.1 * z;                             // r is the center radius, r2 is for damping

    real temp = sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc));
    real scale = fabs(r - temp);

    if (scale > r2)
      scale = 0.8 - scale;
    else
      scale = 1.0;
    const real z0 = std::max(0.0, 0.1 * (0.1 - temp * z));
    temp = sqrt(temp * temp + z0 * z0);
    scale = (r + r2 - temp) * scale / (temp + Globals::SMALL);
    scale = scale / (1 + z);

    v[0] = scale * (y - yc) + 0.1 * (x - xc);
    v[1] = scale * -(x - xc) + 0.1 * (y - yc);
    v[2] = scale * z0;

    return v;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//
