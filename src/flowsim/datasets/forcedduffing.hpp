#pragma once
//--------------------------------------------------------------------------//
#include "../flow.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  class ForcedDuffing2D : public Flow<2>
  //--------------------------------------------------------------------------//
  {
  public:
    //--------------------------------------------------------------------------//
    ForcedDuffing2D();
    //--------------------------------------------------------------------------//
    ForcedDuffing2D(const Domain<2> &space_domain_used);
    //--------------------------------------------------------------------------//
    virtual ~ForcedDuffing2D() {}
    //--------------------------------------------------------------------------//
    virtual MaybeVec<2> v(const VecR<2> &pos, real t) const override;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//