#pragma once
//--------------------------------------------------------------------------//
#include "../types.hpp"
//--------------------------------------------------------------------------//
#include <array>
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  template <uint n>
  class Flow
  {
  public:
    //--------------------------------------------------------------------------//
    Flow(const Domain<n> &space_domain,
         const Range &time_domain,
         const std::string &name = "Flow")
        : m_space_domain_defined(space_domain),
          m_space_domain_standard(space_domain),
          m_time_domain(time_domain),
          m_name(name) {}
    //--------------------------------------------------------------------------//
    Flow(const Domain<n> &space_domain_defined,
         const Domain<n> &space_domain_standard,
         const Range &time_domain,
         const std::string &name = "Flow")
        : m_space_domain_defined(space_domain_defined),
          m_space_domain_standard(space_domain_standard),
          m_time_domain(time_domain),
          m_name(name) {}
    //--------------------------------------------------------------------------//
    virtual ~Flow() {}
    //--------------------------------------------------------------------------//
    virtual MaybeVec<n> v(const VecR<n> &pos, real t) const = 0;
    //--------------------------------------------------------------------------//
    VecR<n> integrate(const VecR<n> &pos, real t, real tau) const
    {
      auto solver = ODE<n>::solver(VC::odeint::RK43, ODE<n>::make_options(VC::odeint::AbsTol = 1e-4,
                                                                          VC::odeint::RelTol = 1e-4,
                                                                          VC::odeint::InitialStep = 1e-2,
                                                                          VC::odeint::MaxStep = 1e-1));
      solver.integrate([this](real t, const VecR<n> &pos) -> MaybeVec<n>
                       { return this->v(pos, t); },
                       t, t + tau, pos);
      return solver.state == VC::odeint::AcceptStep ? solver.y : VecR<n>(NAN);
    }
    //--------------------------------------------------------------------------//
    const std::string &name(void) const { return m_name; }
    //--------------------------------------------------------------------------//
    const Domain<n> &spaceDomainDefined() const { return m_space_domain_defined; }
    const Range &spaceDomainDefined(uint i) const { return m_space_domain_defined[i]; }
    //--------------------------------------------------------------------------//
    const Domain<n> &spaceDomainStandard() const { return m_space_domain_standard; }
    const Range &spaceDomainStandard(uint i) const { return m_space_domain_standard[i]; }
    //--------------------------------------------------------------------------//
    const Range &timeDomain() const { return m_time_domain; }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    Domain<n> m_space_domain_defined;
    Domain<n> m_space_domain_standard;
    Range m_time_domain;
    std::string m_name;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//