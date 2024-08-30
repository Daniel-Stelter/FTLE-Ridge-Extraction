#pragma once
//--------------------------------------------------------------------------//
#include "colors.hpp"
//--------------------------------------------------------------------------//
#include <array>
#include <string>
#include <cmath>
#include <vector>
//--------------------------------------------------------------------------//
namespace dst::images
{
  //--------------------------------------------------------------------------//
  class Texture
  {
    //--------------------------------------------------------------------------//
    std::array<size_t, 2> m_resolution;
    std::vector<color> m_pixel_data;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    Texture(const std::string &filepath);
    Texture(size_t res_u, size_t res_v, const color &col = white);
    Texture(size_t res_u, size_t res_v, const std::vector<color> &data_pixel);
    //--------------------------------------------------------------------------//
    Texture(const Texture &) = default;
    Texture(Texture &&) = default;
    //--------------------------------------------------------------------------//
    Texture &operator=(const Texture &) = default;
    Texture &operator=(Texture &&) = default;
    //--------------------------------------------------------------------------//
    size_t resU() const { return m_resolution[0]; }
    size_t resV() const { return m_resolution[1]; }
    //--------------------------------------------------------------------------//
    /// Samples texture with bilinear interpolation.
    color sample(double u, double v) const;
    //--------------------------------------------------------------------------//
    /// returns pixel at position [u,v] with read/write access
    color &pixel(size_t u, size_t v);
    //--------------------------------------------------------------------------//
    /// returns pixel at position [u,v] with read-only access
    const color &pixel(size_t u, size_t v) const;
    //--------------------------------------------------------------------------//
    /// reads image data data structure
    void readPPM(const std::string &filepath);
    //--------------------------------------------------------------------------//
    /// writes image data into a file with ppm format
    void writePPM(const std::string &filepath) const;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//