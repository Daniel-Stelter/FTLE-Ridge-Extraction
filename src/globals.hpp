#pragma once
//--------------------------------------------------------------------------//
#include "types.hpp"
//--------------------------------------------------------------------------//
namespace dst
{
  //--------------------------------------------------------------------------//
  class Globals
  {
    //--------------------------------------------------------------------------//
  private:
    Globals() {}
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    // smallest possible value difference
    constexpr static real EPS = std::numeric_limits<real>::epsilon();
    // quasi-zero (very small)
    constexpr static real ZERO = 1000 * Globals::EPS;
    // small value, greater than ZERO
    constexpr static real SMALL = 10000000 * Globals::EPS;
    // infinity
    constexpr static real INF = std::numeric_limits<real>::infinity();
    //--------------------------------------------------------------------------//
    // maximum steps when constraining particles
    constexpr static int MAX_CONSTRAIN_STEPS = 100;
    // position of energy potential well
    constexpr static real ENERGY_W = 0.6;
    // depth of energy potential well
    constexpr static real ENERGY_D = -0.02;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//