#pragma once
//--------------------------------------------------------------------------//
#include "texture.hpp"
#include "../types.hpp"
//--------------------------------------------------------------------------//
#include <algorithm>
//--------------------------------------------------------------------------//
namespace dst::images::colormap
{
  //--------------------------------------------------------------------------//
  enum CMap
  {
    CHROMAJS,
    INFERNO,
    VIRIDIS,
    RED_BLUE,
    RING
  };
  //--------------------------------------------------------------------------//
  /// Get the mapped color for a given value between 0 and 1
  color apply(CMap map, double percent);
  //--------------------------------------------------------------------------//
  /// Get the mapped color for a given range and value
  color apply(CMap map, const Range &range, double value);
  //--------------------------------------------------------------------------//
  /// Create a texture of the color mapping with defined width, height and border size
  Texture createTextureRing(size_t width);
  //--------------------------------------------------------------------------//
  /// Create a texture of the color mapping with defined width, height and border size
  Texture createTexture(CMap map, size_t width, size_t height, size_t border_size);
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//