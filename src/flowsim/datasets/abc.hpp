#pragma once
//--------------------------------------------------------------------------//
#include "../flow.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  class ABC3D : public Flow<3>
  //--------------------------------------------------------------------------//
  {
  public:
    //--------------------------------------------------------------------------//
    ABC3D();
    //--------------------------------------------------------------------------//
    ABC3D(const Domain<3> &space_domain_standard);
    //--------------------------------------------------------------------------//
    virtual ~ABC3D() {}
    //--------------------------------------------------------------------------//
    virtual MaybeVec<3> v(const VecR<3> &pos, real t) const override;
    //--------------------------------------------------------------------------//
private:
    //--------------------------------------------------------------------------//
    real a = 3;
    real b = 2;
    real c = 1;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//