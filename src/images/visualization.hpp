#pragma once
//--------------------------------------------------------------------------//
#include "texture.hpp"
#include "colormap.hpp"
#include "../flowsim/scalarfield.hpp"
//---------------------------------------------------------------------------//
#include <functional>
#include <optional>
//---------------------------------------------------------------------------//
namespace dst::images
{
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(
      const flowsim::ScalarFieldBase<2> &scalar_field,
      const VecI<2> &resolution,
      colormap::CMap col_map = colormap::VIRIDIS);
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(
      const flowsim::ScalarFieldBase<2> &scalar_field,
      const VecI<2> &resolution,
      const Domain<2> &render_domain,
      colormap::CMap col_map = colormap::VIRIDIS);
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(
      const flowsim::ScalarFieldBase<2> &scalar_field,
      const VecI<2> &resolution,
      const Range &map_range,
      colormap::CMap col_map = colormap::VIRIDIS);
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(
      const flowsim::ScalarFieldBase<2> &scalar_field,
      const VecI<2> &resolution,
      const Domain<2> &render_domain,
      const Range &map_range,
      colormap::CMap col_map = colormap::VIRIDIS);
  //---------------------------------------------------------------------------//
  void drawParticles(
      Texture &texture,
      const std::vector<VecR<2>> &particles,
      const Domain<2> &render_domain,
      const color &particle_color = black,
      const std::optional<color> &opt_border_color = {});
  //---------------------------------------------------------------------------//
  void drawParticlesWithValue(
      Texture &texture,
      const std::vector<VecR<2>> &particles,
      const std::vector<real> &values,
      const Domain<2> &render_domain,
      colormap::CMap col_map = colormap::INFERNO,
      const std::optional<color> &opt_border_color = {},
      const std::function<real(real)> &func = [](real x)
      { return x; });
  //---------------------------------------------------------------------------//
  void drawPath(
      Texture &texture,
      const std::vector<VecR<2>> &points,
      const Domain<2> &render_domain,
      const color &particle_color = black);
  //---------------------------------------------------------------------------//
  void drawCircle(
      Texture &texture,
      const VecR<2> &center,
      real radius,
      const Domain<2> &render_domain,
      const color &color = black);
  //---------------------------------------------------------------------------//
}
//---------------------------------------------------------------------------//