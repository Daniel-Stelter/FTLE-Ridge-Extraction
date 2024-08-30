#pragma once
//--------------------------------------------------------------------------//
#include "math.hpp"
//--------------------------------------------------------------------------//
#include <vc/odeint.hh>
//--------------------------------------------------------------------------//
#include <Eigen/Dense>
//--------------------------------------------------------------------------//
namespace dst
{
  //--------------------------------------------------------------------------//
  using real = double;
  //--------------------------------------------------------------------------//
  template <uint n>
  using VecR = VC::vecn::vecn<n, real>;
  template <uint n>
  using MaybeVec = VC::odeint::maybe<VecR<n>>;
  template <uint n>
  using VecI = VC::vecn::vecn<n, uint>;
  //--------------------------------------------------------------------------//
  template <uint n>
  using ODE = VC::odeint::ode_t<n, real>;
  template <uint n>
  using Pathline = typename ODE<n>::spline_t;
  //--------------------------------------------------------------------------//
  template <uint n, uint m>
  using EMat = Eigen::Matrix<real, n, m>;
  template <uint n>
  using EVec = Eigen::Vector<real, n>;
  //--------------------------------------------------------------------------//
  using color = VecR<3>;
  //--------------------------------------------------------------------------//
  using VC::odeint::OutOfDomain;
  //--------------------------------------------------------------------------//

  // //--------------------------------------------------------------------------//
  // /** \brief Converts a vector of the Eigen library to a vector of the VC library */
  // template <uint n>
  // VecR<n> cvt(const EVec<n> &v) { return VecR<n>(v.data()); }
  // //--------------------------------------------------------------------------//
  // /** \brief Converts a vector of the VC library to a vector of the Eigen library */
  // template <uint n>
  // EVec<n> cvt(const VecR<n> &v) { return EVec<n>(v.data.data()); }
  // //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  struct Range
  {
    //--------------------------------------------------------------------------//
    Range(real min_val, real max_val) : min(min_val), max(max_val) {}
    //--------------------------------------------------------------------------//
    Range() : Range(-std::numeric_limits<real>::infinity(), std::numeric_limits<real>::infinity()) {}
    //--------------------------------------------------------------------------//
    real min;
    real max;
    real dif() const { return max - min; }
    //--------------------------------------------------------------------------//
    bool isInside(real val) const { return min <= val && val <= max; }
    //--------------------------------------------------------------------------//
    bool operator==(const Range &other) const
    {
      return min == other.min && max == other.max;
    }
    //--------------------------------------------------------------------------//
    operator std::string() const
    {
      return std::string("[") + std::to_string(min) + std::string(", ") +
             std::to_string(max) + std::string("]");
    }
    //--------------------------------------------------------------------------//
    friend std::ostream &operator<<(std::ostream &os, const Range &r)
    {
      return os << "[" << r.min << ", " << r.max << "]";
    }
    //--------------------------------------------------------------------------//
    friend std::istream &operator>>(std::istream &is, Range &r)
    {
      std::string str_min, str_max;
      is >> str_min >> str_max;
      assert(str_min.front() == '[' && str_min.back() == ',' && str_max.back() == ']');
      r.min = std::stod(str_min.substr(1, str_min.length() - 2));
      r.max = std::stod(str_max.substr(0, str_max.length() - 1));
      return is;
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  struct SteppedRange : public Range
  {
    //--------------------------------------------------------------------------//
    uint steps;
    //--------------------------------------------------------------------------//
    SteppedRange() : Range(), steps(0) {}
    //--------------------------------------------------------------------------//
    SteppedRange(real min, real max, uint steps) : Range(min, max), steps(steps) {}
    //--------------------------------------------------------------------------//
    /**
     * Automatically creates a stepped range which goes from zero to max. The
     * step of 0 can be skipped by setting include_zero to false.
    */
    SteppedRange(real max, uint steps, bool include_zero)
        : Range(include_zero ? 0.0 : max / steps, max), steps(steps) {}
    //--------------------------------------------------------------------------//
    real value(uint step) const
    {
      if (steps < 2)
        return this->min;
      return this->min + this->dif() * ((real)step / (steps - 1));
    }
    //--------------------------------------------------------------------------//
    real delta() const
    {
      if (steps < 2)
        return 0.0;
      return this->dif() / (steps - 1);
    }
    //--------------------------------------------------------------------------//
    bool operator==(const SteppedRange &other) const
    {
      return this->min == other.min && this->max == other.max && steps == other.steps;
    }
    //--------------------------------------------------------------------------//
  };

  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  struct Domain : std::array<Range, n>
  {
    //--------------------------------------------------------------------------//
    bool isInside(const VecR<n> &pos) const
    {
      for (uint i = 0; i < n; ++i)
        if (!this->at(i).isInside(pos[i]))
          return false;
      return true;
    }
    //--------------------------------------------------------------------------//
    std::array<real, n> difs() const
    {
      std::array<real, n> a;
      std::transform(this->begin(), this->end(), a.begin(), [](const Range &r)
                     { return r.dif(); });
      return a;
    }
    //--------------------------------------------------------------------------//
    friend std::ostream &operator<<(std::ostream &os, const Domain &d)
    {
      if (n > 0)
      {
        os << d[0];
        for (uint i = 1; i < n; ++i)
          os << " x " << d[i];
      }
      return os;
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//