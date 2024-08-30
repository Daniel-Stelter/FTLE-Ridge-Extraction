#include "visualization.hpp"
//--------------------------------------------------------------------------//
using namespace std;
//---------------------------------------------------------------------------//
namespace dst::images
{
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(const flowsim::ScalarFieldBase<2> &scalar_field,
                                 const VecI<2> &resolution,
                                 colormap::CMap col_map)
  {
    return createScalarFieldImage(scalar_field, resolution, scalar_field.domain(), col_map);
  }
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(const flowsim::ScalarFieldBase<2> &scalar_field,
                                 const VecI<2> &resolution,
                                 const Domain<2> &render_domain,
                                 colormap::CMap col_map)
  {
    // prepare texture
    const color col_NAN = white;
    Texture texture(resolution[0], resolution[1], col_NAN);
    if (resolution[0] < 1 || resolution[1] < 1)
      return texture;

    // calc all values
    flowsim::CellDataGrid<2, real> pixel_grid(resolution, render_domain);
#pragma omp parallel for schedule(dynamic)
    for (uint cell_idx = 0; cell_idx < pixel_grid.numCells(); ++cell_idx)
    {
      VecI<2> cell_coord = pixel_grid.cellIndexToCoord(cell_idx);
      VecR<2> pos = pixel_grid.cellCenterPos(cell_coord);
      pixel_grid.cellData(cell_idx) = scalar_field.value(pos);
    }

    // find map range, but circumvent bug of minmax_element: leave out NAN values
    auto iter_begin = pixel_grid.cellData().begin();
    while (iter_begin != pixel_grid.cellData().end() && *iter_begin != *iter_begin)
      ++iter_begin;
    if (iter_begin == pixel_grid.cellData().end()) // if all values are NAN, just return
      return texture;
    // find min/max values and create range
    auto iters_min_max = minmax_element(iter_begin, pixel_grid.cellData().end());
    Range map_range(*iters_min_max.first, *iters_min_max.second);

    // transform all saved scalar values to corresponding colors
#pragma omp parallel for schedule(dynamic)
    for (uint cell_idx = 0; cell_idx < pixel_grid.numCells(); ++cell_idx)
    {
      VecI<2> cell_coord = pixel_grid.cellIndexToCoord(cell_idx);
      const real sc_val = pixel_grid.cellData(cell_idx);
      if (sc_val == sc_val) // check NAN
        texture.pixel(cell_coord[0], cell_coord[1]) = colormap::apply(col_map, map_range, sc_val);
    }
    return texture;
  }
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(const flowsim::ScalarFieldBase<2> &scalar_field,
                                 const VecI<2> &resolution,
                                 const Range &map_range,
                                 colormap::CMap col_map)
  {
    return createScalarFieldImage(scalar_field, resolution, scalar_field.domain(), map_range, col_map);
  }
  //---------------------------------------------------------------------------//
  Texture createScalarFieldImage(const flowsim::ScalarFieldBase<2> &scalar_field,
                                 const VecI<2> &resolution,
                                 const Domain<2> &render_domain,
                                 const Range &map_range,
                                 colormap::CMap col_map)
  {
    if (resolution[0] < 1 || resolution[1] < 1)
      return Texture(resolution[0], resolution[1]);

    // prepare texture
    const color col_NAN = white;
    Texture texture(resolution[0], resolution[1], col_NAN);

    // calc all values
    flowsim::CellDataGrid<2, real> pixel_grid(resolution, render_domain);
#pragma omp parallel for schedule(dynamic)
    for (uint cell_idx = 0; cell_idx < pixel_grid.numCells(); ++cell_idx)
    {
      VecI<2> cell_coord = pixel_grid.cellIndexToCoord(cell_idx);
      VecR<2> pos = pixel_grid.cellCenterPos(cell_coord);
      pixel_grid.cellData(cell_idx) = scalar_field.value(pos);
    }

    // transform all saved scalar values to corresponding colors
#pragma omp parallel for schedule(dynamic)
    for (uint cell_idx = 0; cell_idx < pixel_grid.numCells(); ++cell_idx)
    {
      VecI<2> cell_coord = pixel_grid.cellIndexToCoord(cell_idx);
      const real sc_val = pixel_grid.cellData(cell_idx);
      if (sc_val == sc_val) // check NAN
        texture.pixel(cell_coord[0], cell_coord[1]) = colormap::apply(col_map, map_range, sc_val);
    }
    return texture;
  }
  //---------------------------------------------------------------------------//
  void drawParticleBorder(Texture &texture,
                          const VecI<2> center,
                          const color &border_color)
  {
    // min statements already handle out of domain
    for (uint u = min(center[0], center[0] - 1); u <= min<uint>(center[0] + 1, texture.resU() - 1); ++u)
      for (uint v = min(center[1], center[1] - 1); v <= min<uint>(center[1] + 1, texture.resV() - 1); ++v)
        if (u != center[0] || v != center[1])
          texture.pixel(u, v) = border_color;
  }
  //---------------------------------------------------------------------------//
  VecI<2> toPixelCoord(const Texture &texture,
                       const VecR<2> &pos,
                       const Domain<2> &render_domain)
  {
    return {min<uint>(texture.resU() - 1,
                      max(0.0, round((pos[0] - render_domain[0].min) / render_domain[0].dif() * texture.resU()))),
            min<uint>(texture.resV() - 1,
                      max(0.0, round((pos[1] - render_domain[1].min) / render_domain[1].dif() * texture.resV())))};
  }
  //---------------------------------------------------------------------------//
  void drawParticles(Texture &texture,
                     const vector<VecR<2>> &particles,
                     const Domain<2> &render_domain,
                     const color &particle_color,
                     const optional<color> &opt_border_color)
  {
    for (const VecR<2> &par : particles)
    {
      auto uv = toPixelCoord(texture, par, render_domain);
      texture.pixel(uv[0], uv[1]) = particle_color;
      if (opt_border_color)
        drawParticleBorder(texture, uv, *opt_border_color);
    }
  }
  //---------------------------------------------------------------------------//
  void drawParticlesWithValue(Texture &texture,
                              const vector<VecR<2>> &particles,
                              const vector<real> &values,
                              const Domain<2> &render_domain,
                              colormap::CMap col_map,
                              const optional<color> &opt_border_color,
                              const function<real(real)> &func)
  {
    assert(particles.size() == values.size());
    if (particles.size() == 0)
      return;

    // find max value, apply func to it and create its inverse
    real max_val = func(*max_element(values.begin(), values.end()));
    real inv_max_val = 1.0 / max_val;

#pragma omp parallel for
    for (uint i = 0; i < particles.size(); ++i)
    {
      VecI<2> uv = toPixelCoord(texture, particles[i], render_domain);
      real val = func(values[i]);

      texture.pixel(uv[0], uv[1]) = colormap::apply(col_map, val / inv_max_val);
      if (opt_border_color)
        drawParticleBorder(texture, uv, *opt_border_color);
    }
  }
  //---------------------------------------------------------------------------//

  //---------------------------------------------------------------------------//
  vector<VecI<2>> rasterize(const VecI<2> &p1, const VecI<2> &p2)
  {
    // trivial case
    if (p1 == p2)
      return {p1};
    // setup
    int dir_x = p1[0] <= p2[0] ? 1 : -1;
    int dir_y = p1[1] <= p2[1] ? 1 : -1;
    int delta_x = abs((int)p2[0] - (int)p1[0]);
    int delta_y = abs((int)p2[1] - (int)p1[1]);
    int dim_x = 0;
    int dim_y = 1;
    // switch if necessary
    if (delta_x < delta_y)
    {
      swap(dir_x, dir_y);
      swap(delta_x, delta_y);
      swap(dim_x, dim_y);
    }
    // create output
    vector<VecI<2>> path;
    path.reserve(abs(delta_x) + 1);
    path.push_back(p1);
    // actual execution
    int x = p1[dim_x];
    int y = p1[dim_y];
    int d = 2 * delta_y - delta_x;
    while (x != (int)p2[dim_x])
    {
      x += dir_x;
      if (d > 0)
      {
        y += dir_y;
        d += 2 * delta_y - 2 * delta_x;
      }
      else if (d <= 0)
        d += 2 * delta_y;
      path.push_back(dim_x == 0 ? VecI<2>(x, y) : VecI<2>(y, x));
    }
    return path;
  }
  //---------------------------------------------------------------------------//
  void drawPath(Texture &texture,
                const vector<VecR<2>> &points,
                const Domain<2> &render_domain,
                const color &particle_color)
  {
    for (uint j = 0; j < points.size() - 1; ++j)
    {
      VecI<2> p1 = toPixelCoord(texture, points[j], render_domain);
      VecI<2> p2 = toPixelCoord(texture, points[j + 1], render_domain);
      vector<VecI<2>> path = rasterize(p1, p2);
      for (const VecI<2> pixel : path)
        texture.pixel(pixel[0], pixel[1]) = particle_color;
    }
  }
  //---------------------------------------------------------------------------//
  void drawCircle(Texture &texture,
                  const VecR<2> &real_center,
                  real real_radius,
                  const Domain<2> &render_domain,
                  const color &color)
  {
    VecI<2> discrete_center = toPixelCoord(texture, real_center, render_domain);
    const int u_center = discrete_center[0];
    const int v_center = discrete_center[1];
    const int radius = (real_radius / render_domain[0].dif()) * texture.resU();

    auto setPixel = [&](int u, int v)
    {
      if (0 <= u && u < (int)texture.resU() && 0 <= v && v < (int)texture.resV())
        texture.pixel(u, v) = color;
    };

    int u = radius;
    int v = 0;
    int p = 1 - radius; // initial decision parameter for the circle

    // printing initial points (if radius is zero, only a single point will be printed)
    setPixel(u_center + radius, v_center);
    if (radius > 0)
    {
      setPixel(u_center - radius, v_center);
      setPixel(u_center, v_center + radius);
      setPixel(u_center, v_center - radius);
    }
    // iterate over whole
    while (u > v++)
    {
      // mid-point is inside or on the perimeter of the circle
      if (p <= 0)
        p = p + 2 * v + 1;
      // mid-point is outside the perimeter of the circle
      else
      {
        --u;
        p = p + 2 * v - 2 * u + 1;
      }
      // all the perimeter points have already been printed
      if (u < v)
        break;

      // draw generated point on all eight parts
      setPixel(u_center + u, v_center - v);
      setPixel(u_center - u, v_center - v);
      setPixel(u_center + u, v_center + v);
      setPixel(u_center - u, v_center + v);
      // check if already drawn
      if (u != v)
      {
        setPixel(u_center + v, v_center - u);
        setPixel(u_center - v, v_center - u);
        setPixel(u_center + v, v_center + u);
        setPixel(u_center - v, v_center + u);
      }
    }
  }
  //---------------------------------------------------------------------------//
}
//---------------------------------------------------------------------------//