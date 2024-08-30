#include "doublegyre.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  DoubleGyre3D::DoubleGyre3D()
      : DoubleGyre3D({Range(0, 2), Range(0, 10), Range(0, 1)}) {}
  //--------------------------------------------------------------------------//
  DoubleGyre3D::DoubleGyre3D(const Domain<3> &space_domain_used)
      : Flow<3>({Range(), Range(), Range()},
                space_domain_used,
                Range(),
                "DoubleGyre3D") {}
  //--------------------------------------------------------------------------//
  MaybeVec<3>
  DoubleGyre3D::v(const VecR<3> &pos, real t) const
  {
    const real &x = pos[0];
    const real &y = pos[1];
    const real &z = pos[2];
    VecR<3> v;

    // use y-axis for an interval of [t0, t0 + 10]
    t += y;

    const real &A = m_A;
    const real &eps = m_eps;
    const real &omega = m_omega;

    const real a = eps * sin(omega * t);
    const real b = 1.0 - 2.0 * a;
    const real f = a * x * x + b * x;

    v[0] = -M_PI * A * sin(M_PI * f) * cos(M_PI * z);
    v[1] = 0;
    v[2] = M_PI * A * cos(M_PI * f) * sin(M_PI * z) * (2.0 * a * x + b);

    return v;
  }
  //--------------------------------------------------------------------------//
  void
  DoubleGyre3D::setParameters(real A, real omega, real eps)
  {
    m_A = A;
    m_omega = omega;
    m_eps = eps;
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  DoubleGyre2D::DoubleGyre2D()
      : DoubleGyre2D({Range(0, 2), Range(0, 1)}) {}
  //--------------------------------------------------------------------------//
  DoubleGyre2D::DoubleGyre2D(const Domain<2> &space_domain_used)
      : Flow<2>({Range(), Range()},
                space_domain_used,
                Range(),
                "DoubleGyre2D") {}
  //--------------------------------------------------------------------------//
  MaybeVec<2>
  DoubleGyre2D::v(const VecR<2> &pos, real t) const
  {
    real x = pos[0];
    real y = pos[1];
    VecR<2> v;

    const real &A = m_A;
    const real &eps = m_eps;
    const real &omega = m_omega;

    const real a = eps * sin(omega * t);
    const real b = 1.0 - 2.0 * a;
    const real f = a * x * x + b * x;

    v[0] = -M_PI * A * sin(M_PI * f) * cos(M_PI * y);
    v[1] = M_PI * A * cos(M_PI * f) * sin(M_PI * y) * (2.0 * a * x + b);

    return v;
  }
  //--------------------------------------------------------------------------//
  void
  DoubleGyre2D::setParameters(real A, real omega, real eps)
  {
    m_A = A;
    m_omega = omega;
    m_eps = eps;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//
