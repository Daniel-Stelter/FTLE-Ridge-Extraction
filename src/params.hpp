#pragma once
//--------------------------------------------------------------------------//
#include "types.hpp"
#include "utils.hpp"
//--------------------------------------------------------------------------//
#include <optional>
//--------------------------------------------------------------------------//
namespace dst::params
{
  //--------------------------------------------------------------------------//
  template <uint n>
  struct ParamsGeneral
  {
    //--------------------------------------------------------------------------//
    /**
     * Spatial domain for the scalar field. If no domain is provided, it uses
     * the standard domain of the scalar field. For FTLE computation, this refers
     * to the standard domain of the flow.
     */
    std::optional<Domain<n>> domain;
    /**
     * Start time for the series. In case of FTLE computation, this is the start
     * time of the integration of path lines.
     */
    real t;
    /**
     * Time period for the series, i.e. the time range will be [t, t + tau]. For
     * FTLE computation, the start time t will stay the same. Instead, tau refers
     * to the integration time of the flow in the interval [tau/steps, tau]. This
     * prevents the first step to be tau=0.0.
     */
    real tau;
    /**
     * Number of steps for tau for which ridges will be searched.
     */
    uint steps;
    /**
     * Resolution for the interpolation grid of the scalar field.
     */
    VecI<n> space_res;
    //--------------------------------------------------------------------------//
    friend std::ostream &operator<<(std::ostream &os, const ParamsGeneral &p)
    {
      os << "GENERAL PARAMS:\n";
      if (p.domain.has_value())
        utils::prettyPrintListItem("Domain", *p.domain, 26);
      else
        utils::prettyPrintListItem("Domain", "standard", 26);
      utils::prettyPrintListItem("t", p.t, 26);
      utils::prettyPrintListItem("tau", p.tau, 26);
      utils::prettyPrintListItem("steps", p.steps, 26);
      utils::prettyPrintListItem("Space resolution", p.space_res, 26, false);
      return os;
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  struct ParamsStelter
  {
    //--------------------------------------------------------------------------//
    /**
     * Resolution for the initialization grid of new particles.
     */
    VecI<n> par_init_res;
    /**
     * Resolution of the voxel grid which is responsible for efficiently finding
     * particles in the execution.
     */
    VecI<n> voxel_res;
    /**
     * Threshold for the ridge strength of a particle.
     */
    real min_ridge_strength;
    /**
     * Periodicity for full initialization of new particles to find new emerging
     * ridges. Value 1 would lead to initialization in each iteration.
     */
    uint full_init_period;
    /**
     * Factor for elliptic / anisotropic distance in ridge tangent direction(s).
     * Value 1.0 leads to Euclidean distance.
     */
    real elliptic_dis_fac;
    //--------------------------------------------------------------------------//
    friend std::ostream &operator<<(std::ostream &os, const ParamsStelter &p)
    {
      os << "PARAMS STELTER:\n";
      utils::prettyPrintListItem("Particle init resolution", p.par_init_res, 26);
      utils::prettyPrintListItem("Voxel resolution", p.voxel_res, 26);
      utils::prettyPrintListItem("Min ridge strength", p.min_ridge_strength, 26);
      utils::prettyPrintListItem("Elliptic distance factor", p.elliptic_dis_fac, 26);
      utils::prettyPrintListItem("Full init period", p.full_init_period, 26, false);
      return os;
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  struct ParamsKindlmann
  {
    //--------------------------------------------------------------------------//
    /**
     * Resolution for the initialization grid of new particles.
     */
    VecI<n> par_init_res;
    /**
     * Resolution of the voxel grid which is responsible for efficiently finding
     * particles in the execution.
     */
    VecI<n> voxel_res;
    /**
     * Threshold for the ridge strength of a particle.
     */
    real min_ridge_strength;
    /**
     * Radius for the neighborhood space of a particle.
     */
    real sigma;
    /**
     * Difference for numerical computation of gradients and Hessians.
     */
    real numeric_diff_delta;
    //--------------------------------------------------------------------------//
    friend std::ostream &operator<<(std::ostream &os, const ParamsKindlmann &p)
    {
      os << "PARAMS KINDLMANN:\n";
      utils::prettyPrintListItem("Particle init resolution", p.par_init_res, 26);
      utils::prettyPrintListItem("Voxel resolution", p.voxel_res, 26);
      utils::prettyPrintListItem("Min ridge strength", p.min_ridge_strength, 26);
      utils::prettyPrintListItem("Sigma", p.sigma, 26);
      utils::prettyPrintListItem("NumericDiff delta", p.numeric_diff_delta, 26, false);
      return os;
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  struct ParamsMarRidges
  {
    //--------------------------------------------------------------------------//
    /**
     * Threshold for the ridge strength of a particle.
     */
    real min_ridge_strength;
    //--------------------------------------------------------------------------//
    friend std::ostream &operator<<(std::ostream &os, const ParamsMarRidges &p)
    {
      os << "PARAMS MARCHING RIDGES:\n";
      utils::prettyPrintListItem("Min ridge strength", p.min_ridge_strength, 26, false);
      return os;
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//