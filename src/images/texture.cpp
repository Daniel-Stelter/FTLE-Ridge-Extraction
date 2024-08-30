#include "texture.hpp"
//--------------------------------------------------------------------------//
#include "../types.hpp"
//--------------------------------------------------------------------------//
#include <cassert>
#include <filesystem>
#include <fstream>
//-------------------------------------------------------------------------//
namespace dst::images
{
  //-------------------------------------------------------------------------//
  Texture::Texture(const std::string &filepath) : m_resolution{0, 0}
  {
    readPPM(filepath);
  }

  //-------------------------------------------------------------------------//
  Texture::Texture(size_t res_u, size_t res_v, const color &col)
      : m_resolution{res_u, res_v}, m_pixel_data(res_u * res_v, col) {}

  //-------------------------------------------------------------------------//
  Texture::Texture(size_t res_u, size_t res_v,
                   const std::vector<color> &pixel_data)
      : m_resolution{res_u, res_v}, m_pixel_data(pixel_data)
  {
    assert(res_u * res_v == pixel_data.size());
  }

  //-------------------------------------------------------------------------//
  color Texture::sample(double u, double v) const
  {
    // scale position to have u in range of [0, res_u - 1]
    u *= resU() - 1;

    // scale position to have v in range of [0, res_v - 1]
    v *= resV() - 1;

    //             ut
    //      lu      |     ru
    //   tv []---- p1 ----[] tv
    //      |       .      |
    //      |       .      |
    //  vt -|......p2......|- vt
    //      |       .      |
    //      |       .      |
    //   bv []---- p0 ----[] bv
    //      lu      |     ru
    //             ut
    //
    // [] are given pixel values stored in m_pixel_data.
    // lu, ru, bv, tv are left, right, bottom and top pixel indices.
    // ut and vt are interpolation factors.
    // p0 and p1 are linearly interpolated pixel colors.
    // p2 is final bilinearly interpolated pixel color.

    // get left and right pixel indices.
    // using modulo we get a texture with repeat behavior
    size_t lu =
        std::abs(static_cast<int>(std::floor(u)) % static_cast<int>(resU()));
    size_t ru = std::min(lu + 1, resU() - 1);

    // get bottom and top pixel indices
    // using modulo we get a texture with repeat behavior
    size_t bv =
        std::abs(static_cast<int>(std::floor(v)) % static_cast<int>(resV()));
    size_t tv = std::min(bv + 1, resV() - 1);

    // interpolation factor for first two linear interpolations
    double ut = u - std::floor(u);
    // interpolation factor for final bilinear interpolations
    double vt = v - std::floor(v);

    // bottom linear interpolation
    color p0 = pixel(lu, bv) * (1 - ut) + pixel(ru, bv) * ut;

    // top linear interpolation
    color p1 = pixel(lu, tv) * (1 - ut) + pixel(ru, tv) * ut;

    // final biliniear interpolation
    color p2 = p0 * (1 - vt) + p1 * vt;

    return p2;
  }

  //-------------------------------------------------------------------------//
  color &Texture::pixel(size_t u, size_t v)
  {
    return m_pixel_data[u + v * m_resolution[0]];
  }

  //-------------------------------------------------------------------------//
  const color &Texture::pixel(size_t u, size_t v) const
  {
    return m_pixel_data[u + v * m_resolution[0]];
  }

  //-------------------------------------------------------------------------//
  void Texture::readPPM(const std::string &file_path)
  {
    using namespace std;
    ifstream in(file_path, ios::binary);
    if (!in)
      throw std::runtime_error("Texture: could not read from file '" + file_path + "'");

    string P6;
    size_t max_val;
    in >> P6 >> m_resolution[0] >> m_resolution[1] >> max_val;
    in.ignore(1);
    if ("P6" != P6)
      throw std::runtime_error("Texture: file does not start with 'P6' (for PPM files)");
    m_pixel_data.resize(m_resolution[0] * m_resolution[1]);

    uint8_t vals[3];
    for (size_t v = 0; v < resV(); ++v)
    {
      for (size_t u = 0; u < resU(); ++u)
      {
        in.read((char *)&vals[0], sizeof(vals));
        pixel(u, resV() - v - 1)[0] = (real)vals[0] / max_val;
        pixel(u, resV() - v - 1)[1] = (real)vals[1] / max_val;
        pixel(u, resV() - v - 1)[2] = (real)vals[2] / max_val;
      }
    }
    in.close();
  }

  //-------------------------------------------------------------------------//
  void Texture::writePPM(const std::string &file_path) const
  {
    using namespace std;
    filesystem::create_directories(filesystem::absolute(file_path).parent_path());
    ofstream out(file_path, ios::binary);
    if (!out)
      throw runtime_error("Texture: could not write to file '" + file_path + "'");

    const size_t max_val = 255;
    auto remap = [&](auto u, auto v, auto i)
    {
      const auto &comp = pixel(u, resV() - v - 1)[i];
      return (uint8_t)abs(floor(max(0.0, min(1.0, comp)) * max_val));
    };
    if (out)
    {
      out << "P6 " << resU() << " " << resV() << " " << max_val << "\n";
      for (size_t v = 0; v < resV(); ++v)
      {
        for (size_t u = 0; u < resU(); ++u)
        {
          uint8_t vals[3] = {remap(u, v, 0), remap(u, v, 1), remap(u, v, 2)};
          out.write((char *)&vals[0], sizeof(vals));
        }
      }
      out.close();
    }
  }
  //-------------------------------------------------------------------------//
}
//-------------------------------------------------------------------------//