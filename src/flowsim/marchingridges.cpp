#include "marchingridges.hpp"
//--------------------------------------------------------------------------//
#include "ridgeinfo.hpp"
//--------------------------------------------------------------------------//
#include <numeric>
//--------------------------------------------------------------------------//
using namespace std;
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  // binary encoding of required edges for a given Marching Cubes case
  const uint EDGE_LOOKUP_TABLE[256]{
      0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
      0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
      0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
      0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
      0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
      0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
      0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
      0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
      0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
      0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
      0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
      0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
      0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
      0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
      0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
      0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
      0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
      0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
      0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
      0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
      0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
      0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
      0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
      0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
      0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
      0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
      0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
      0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
      0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
      0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
      0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
      0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0};

  //--------------------------------------------------------------------------//
  // edge indices for triangulation for a given Marching Cubes case
  const vector<uint> TRIANGLE_LOOKUP_TABLE[256]{
      {},
      {0, 8, 3},
      {0, 1, 9},
      {1, 8, 3, 9, 8, 1},
      {1, 2, 10},
      {0, 8, 3, 1, 2, 10},
      {9, 2, 10, 0, 2, 9},
      {2, 8, 3, 2, 10, 8, 10, 9, 8},
      {3, 11, 2},
      {0, 11, 2, 8, 11, 0},
      {1, 9, 0, 2, 3, 11},
      {1, 11, 2, 1, 9, 11, 9, 8, 11},
      {3, 10, 1, 11, 10, 3},
      {0, 10, 1, 0, 8, 10, 8, 11, 10},
      {3, 9, 0, 3, 11, 9, 11, 10, 9},
      {9, 8, 10, 10, 8, 11},
      {4, 7, 8},
      {4, 3, 0, 7, 3, 4},
      {0, 1, 9, 8, 4, 7},
      {4, 1, 9, 4, 7, 1, 7, 3, 1},
      {1, 2, 10, 8, 4, 7},
      {3, 4, 7, 3, 0, 4, 1, 2, 10},
      {9, 2, 10, 9, 0, 2, 8, 4, 7},
      {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4},
      {8, 4, 7, 3, 11, 2},
      {11, 4, 7, 11, 2, 4, 2, 0, 4},
      {9, 0, 1, 8, 4, 7, 2, 3, 11},
      {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1},
      {3, 10, 1, 3, 11, 10, 7, 8, 4},
      {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4},
      {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3},
      {4, 7, 11, 4, 11, 9, 9, 11, 10},
      {9, 5, 4},
      {9, 5, 4, 0, 8, 3},
      {0, 5, 4, 1, 5, 0},
      {8, 5, 4, 8, 3, 5, 3, 1, 5},
      {1, 2, 10, 9, 5, 4},
      {3, 0, 8, 1, 2, 10, 4, 9, 5},
      {5, 2, 10, 5, 4, 2, 4, 0, 2},
      {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8},
      {9, 5, 4, 2, 3, 11},
      {0, 11, 2, 0, 8, 11, 4, 9, 5},
      {0, 5, 4, 0, 1, 5, 2, 3, 11},
      {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5},
      {10, 3, 11, 10, 1, 3, 9, 5, 4},
      {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10},
      {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3},
      {5, 4, 8, 5, 8, 10, 10, 8, 11},
      {9, 7, 8, 5, 7, 9},
      {9, 3, 0, 9, 5, 3, 5, 7, 3},
      {0, 7, 8, 0, 1, 7, 1, 5, 7},
      {1, 5, 3, 3, 5, 7},
      {9, 7, 8, 9, 5, 7, 10, 1, 2},
      {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3},
      {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2},
      {2, 10, 5, 2, 5, 3, 3, 5, 7},
      {7, 9, 5, 7, 8, 9, 3, 11, 2},
      {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11},
      {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7},
      {11, 2, 1, 11, 1, 7, 7, 1, 5},
      {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11},
      {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0},
      {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0},
      {11, 10, 5, 7, 11, 5},
      {10, 6, 5},
      {0, 8, 3, 5, 10, 6},
      {9, 0, 1, 5, 10, 6},
      {1, 8, 3, 1, 9, 8, 5, 10, 6},
      {1, 6, 5, 2, 6, 1},
      {1, 6, 5, 1, 2, 6, 3, 0, 8},
      {9, 6, 5, 9, 0, 6, 0, 2, 6},
      {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8},
      {2, 3, 11, 10, 6, 5},
      {11, 0, 8, 11, 2, 0, 10, 6, 5},
      {0, 1, 9, 2, 3, 11, 5, 10, 6},
      {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11},
      {6, 3, 11, 6, 5, 3, 5, 1, 3},
      {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6},
      {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9},
      {6, 5, 9, 6, 9, 11, 11, 9, 8},
      {5, 10, 6, 4, 7, 8},
      {4, 3, 0, 4, 7, 3, 6, 5, 10},
      {1, 9, 0, 5, 10, 6, 8, 4, 7},
      {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4},
      {6, 1, 2, 6, 5, 1, 4, 7, 8},
      {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7},
      {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6},
      {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9},
      {3, 11, 2, 7, 8, 4, 10, 6, 5},
      {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11},
      {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6},
      {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6},
      {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6},
      {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11},
      {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7},
      {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9},
      {10, 4, 9, 6, 4, 10},
      {4, 10, 6, 4, 9, 10, 0, 8, 3},
      {10, 0, 1, 10, 6, 0, 6, 4, 0},
      {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10},
      {1, 4, 9, 1, 2, 4, 2, 6, 4},
      {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4},
      {0, 2, 4, 4, 2, 6},
      {8, 3, 2, 8, 2, 4, 4, 2, 6},
      {10, 4, 9, 10, 6, 4, 11, 2, 3},
      {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6},
      {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10},
      {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1},
      {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3},
      {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1},
      {3, 11, 6, 3, 6, 0, 0, 6, 4},
      {6, 4, 8, 11, 6, 8},
      {7, 10, 6, 7, 8, 10, 8, 9, 10},
      {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10},
      {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0},
      {10, 6, 7, 10, 7, 1, 1, 7, 3},
      {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7},
      {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9},
      {7, 8, 0, 7, 0, 6, 6, 0, 2},
      {7, 3, 2, 6, 7, 2},
      {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7},
      {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7},
      {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11},
      {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1},
      {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6},
      {0, 9, 1, 11, 6, 7},
      {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0},
      {7, 11, 6},
      {7, 6, 11},
      {3, 0, 8, 11, 7, 6},
      {0, 1, 9, 11, 7, 6},
      {8, 1, 9, 8, 3, 1, 11, 7, 6},
      {10, 1, 2, 6, 11, 7},
      {1, 2, 10, 3, 0, 8, 6, 11, 7},
      {2, 9, 0, 2, 10, 9, 6, 11, 7},
      {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8},
      {7, 2, 3, 6, 2, 7},
      {7, 0, 8, 7, 6, 0, 6, 2, 0},
      {2, 7, 6, 2, 3, 7, 0, 1, 9},
      {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6},
      {10, 7, 6, 10, 1, 7, 1, 3, 7},
      {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8},
      {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7},
      {7, 6, 10, 7, 10, 8, 8, 10, 9},
      {6, 8, 4, 11, 8, 6},
      {3, 6, 11, 3, 0, 6, 0, 4, 6},
      {8, 6, 11, 8, 4, 6, 9, 0, 1},
      {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6},
      {6, 8, 4, 6, 11, 8, 2, 10, 1},
      {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6},
      {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9},
      {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3},
      {8, 2, 3, 8, 4, 2, 4, 6, 2},
      {0, 4, 2, 4, 6, 2},
      {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8},
      {1, 9, 4, 1, 4, 2, 2, 4, 6},
      {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1},
      {10, 1, 0, 10, 0, 6, 6, 0, 4},
      {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3},
      {10, 9, 4, 6, 10, 4},
      {4, 9, 5, 7, 6, 11},
      {0, 8, 3, 4, 9, 5, 11, 7, 6},
      {5, 0, 1, 5, 4, 0, 7, 6, 11},
      {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5},
      {9, 5, 4, 10, 1, 2, 7, 6, 11},
      {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5},
      {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2},
      {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6},
      {7, 2, 3, 7, 6, 2, 5, 4, 9},
      {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7},
      {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0},
      {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8},
      {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7},
      {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4},
      {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10},
      {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10},
      {6, 9, 5, 6, 11, 9, 11, 8, 9},
      {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5},
      {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11},
      {6, 11, 3, 6, 3, 5, 5, 3, 1},
      {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6},
      {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10},
      {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5},
      {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3},
      {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2},
      {9, 5, 6, 9, 6, 0, 0, 6, 2},
      {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8},
      {1, 5, 6, 2, 1, 6},
      {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6},
      {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0},
      {0, 3, 8, 5, 6, 10},
      {10, 5, 6},
      {11, 5, 10, 7, 5, 11},
      {11, 5, 10, 11, 7, 5, 8, 3, 0},
      {5, 11, 7, 5, 10, 11, 1, 9, 0},
      {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1},
      {11, 1, 2, 11, 7, 1, 7, 5, 1},
      {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11},
      {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7},
      {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2},
      {2, 5, 10, 2, 3, 5, 3, 7, 5},
      {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5},
      {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2},
      {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2},
      {1, 3, 5, 3, 7, 5},
      {0, 8, 7, 0, 7, 1, 1, 7, 5},
      {9, 0, 3, 9, 3, 5, 5, 3, 7},
      {9, 8, 7, 5, 9, 7},
      {5, 8, 4, 5, 10, 8, 10, 11, 8},
      {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0},
      {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5},
      {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4},
      {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8},
      {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11},
      {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5},
      {9, 4, 5, 2, 11, 3},
      {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4},
      {5, 10, 2, 5, 2, 4, 4, 2, 0},
      {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9},
      {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2},
      {8, 4, 5, 8, 5, 3, 3, 5, 1},
      {0, 4, 5, 1, 0, 5},
      {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5},
      {9, 4, 5},
      {4, 11, 7, 4, 9, 11, 9, 10, 11},
      {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11},
      {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11},
      {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4},
      {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2},
      {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3},
      {11, 7, 4, 11, 4, 2, 2, 4, 0},
      {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4},
      {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9},
      {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7},
      {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10},
      {1, 10, 2, 8, 7, 4},
      {4, 9, 1, 4, 1, 7, 7, 1, 3},
      {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1},
      {4, 0, 3, 7, 4, 3},
      {4, 8, 7},
      {9, 10, 8, 10, 11, 8},
      {3, 0, 9, 3, 9, 11, 11, 9, 10},
      {0, 1, 10, 0, 10, 8, 8, 10, 11},
      {3, 1, 10, 11, 3, 10},
      {1, 2, 11, 1, 11, 9, 9, 11, 8},
      {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9},
      {0, 2, 11, 8, 0, 11},
      {3, 2, 11},
      {2, 3, 8, 2, 8, 10, 10, 8, 9},
      {9, 10, 2, 0, 9, 2},
      {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8},
      {1, 10, 2},
      {1, 3, 8, 9, 1, 8},
      {0, 9, 1},
      {0, 3, 8},
      {}};

  //--------------------------------------------------------------------------//
  MarchingRidges3D::MarchingRidges3D(ItplScField_ptr_t p_itpl_sc_field,
                                     real min_ridge_strength)
      : mp_sc_field_base(p_itpl_sc_field),
        m_grid(p_itpl_sc_field->cellResolution(), p_itpl_sc_field->domain()),
        m_num_diff(p_itpl_sc_field, p_itpl_sc_field->cellSideLengths()),
        m_min_ridge_strength(min_ridge_strength) {}

  //--------------------------------------------------------------------------//
  MarchingRidges3D::MarchingRidges3D(ItplScField_ptr_t p_itpl_sc_field,
                                     const Params_t &params)
      : MarchingRidges3D(p_itpl_sc_field, params.min_ridge_strength) {}

  //--------------------------------------------------------------------------//
  MarchingRidges3D::MarchingRidges3D(ScFieldBase_ptr_t p_sc_field_base,
                                     real min_ridge_strength,
                                     const VecI<3> &cell_res)
      : mp_sc_field_base(p_sc_field_base),
        m_grid(cell_res, p_sc_field_base->domain()),
        m_num_diff(p_sc_field_base, m_grid.cellSideLengths()),
        m_min_ridge_strength(min_ridge_strength) {}

  //--------------------------------------------------------------------------//
  void MarchingRidges3D::calcRidgeSurface()
  {
    m_vertices.clear();
    m_edge_to_vertex_map.clear();
    m_triangles.clear();
#pragma omp parallel for schedule(dynamic)
    for (uint cell_idx = 0; cell_idx < m_grid.numCells(); ++cell_idx)
      findCellTriangles(m_grid.cellIndexToCoord(cell_idx));
  }

  //--------------------------------------------------------------------------//
  void MarchingRidges3D::saveRidgeSurface(const string &file_path) const
  {
    filesystem::create_directories(filesystem::absolute(file_path).parent_path());

    using UIntBlender_t = u_int32_t;
    using FloatBlender_t = float;
    using VecBlender_t = VC::vecn::vecn<3, FloatBlender_t>;
    using TriIndicesBlender_t = VC::vecn::vecn<3, UIntBlender_t>;

    ofstream out(file_path, ios::binary);
    if (!out)
      throw runtime_error("MarchingRidges3D: could not write to file '" + file_path + "'");

    // write number of vertices and triangles
    UIntBlender_t num_vertices = m_vertices.size();
    out.write((char *)&num_vertices, sizeof(UIntBlender_t));
    UIntBlender_t num_triangles = m_triangles.size();
    out.write((char *)&num_triangles, sizeof(UIntBlender_t));

    // precalculate all scalar values
    vector<FloatBlender_t> vals(num_vertices);
#pragma omp parallel for schedule(dynamic)
    for (VertIndex vert_idx = 0; vert_idx < num_vertices; ++vert_idx)
      vals[vert_idx] = mp_sc_field_base->value(m_vertices[vert_idx]);

    // find and write min and max values
    auto iters_minmax = minmax_element(vals.begin(), vals.end());
    FloatBlender_t val_min = iters_minmax.first == vals.end() ? -Globals::INF : *iters_minmax.first;
    FloatBlender_t val_max = iters_minmax.second == vals.end() ? Globals::INF : *iters_minmax.second;
    out.write((char *)&val_min, sizeof(FloatBlender_t));
    out.write((char *)&val_max, sizeof(FloatBlender_t));

    // write all vertex data with scalar values
    for (VertIndex vert_idx = 0; vert_idx < num_vertices; ++vert_idx)
    {
      VecBlender_t pos(m_vertices[vert_idx][0], m_vertices[vert_idx][1], m_vertices[vert_idx][2]);
      out.write((char *)&pos, sizeof(VecBlender_t));
      out.write((char *)&vals[vert_idx], sizeof(FloatBlender_t));
    }

    // write all triangle data
    for (TriangleIndex tri_idx = 0; tri_idx < num_triangles; ++tri_idx)
    {
      TriIndicesBlender_t tri_indices{m_triangles[tri_idx][0], m_triangles[tri_idx][1], m_triangles[tri_idx][2]};
      out.write((char *)&tri_indices, sizeof(TriIndicesBlender_t));
    }

    out.close();
  }

  //--------------------------------------------------------------------------//
  optional<VecR<3>> MarchingRidges3D::calcLineIntersection(const VecR<3> &pos1,
                                                           real grad1,
                                                           const VecR<3> &pos2,
                                                           real grad2,
                                                           const VecR<3> &avg_transverse) const
  {
    // there can only be an intersection if there are different signs
    if ((grad1 > 0 && grad2 > 0) || (grad1 < 0 && grad2 < 0))
      return {};
    // assume a linear behavior and find relative position
    real rel_pos = grad1 == 0 && grad2 == 0 ? 0.5 : (grad2 / (grad2 - grad1));
    VecR<3> pos_ridge = rel_pos * pos1 + (1 - rel_pos) * pos2;
    // check if second derivative is negative
    VecR<3> transverse_ridge = VC::vecn::normalized(ridgeinfo::stelter::direction<3>(m_num_diff, pos_ridge));
    EMat<3, 3> hessian_ridge = m_num_diff.hessian(pos_ridge);
    if ((transverse_ridge | avg_transverse) < 0)
      transverse_ridge *= -1;
    EVec<3> temp(transverse_ridge.data.data());
    real curvature_ridge = (temp.transpose() * hessian_ridge * temp);
    // ridge strength is the negative of the curvature, thus the minus
    if (-curvature_ridge >= m_min_ridge_strength)
      return pos_ridge;
    return {};
  }

  //--------------------------------------------------------------------------//
  VecR<3> calcAverageTransverseDirection(const array<VecR<3>, 8> &transverse_vecs)
  {
    assert(!transverse_vecs.empty());
    // create average transverse direction
    EMat<3, 3> avg_mat = Eigen::MatrixXd::Zero(3, 3);
    for (const VecR<3> &t_vec : transverse_vecs)
    {
      EVec<3> temp(t_vec.data.data());
      avg_mat += temp * temp.transpose();
    }
    avg_mat /= transverse_vecs.size();
    // find average as greatest eigenvector of matrix
    Eigen::SelfAdjointEigenSolver<EMat<3, 3>> solver;
    solver.compute(avg_mat);
    return VecR<3>(solver.eigenvectors().col(2).data());
  }

  //--------------------------------------------------------------------------//
  /**
   * Note: The used lookup table and the implemented grid use different indexing
   * for cells and edges. Especially, the lookup table only considers local
   * indices which are structured in the following way:
   *     7-----------6-----------6
   *    /|                      /|
   *   7 |                     5 |
   *  /  |                    /  |
   * 4-----------4-----------5   |
   * |   11                  |   10
   * |   |                   |   |
   * |   |                   |   |
   * |   |                   |   |
   * 8   |                   9   |
   * |   3-----------2-------|---2
   * |  /                    |  /
   * | 3                     | 1
   * |/                      |/
   * 0-----------0-----------1
   * The following functions help to convert between both structures.
   */

  const array<VecI<3>, 8> LOCAL_POS_COORDS{
      VecI<3>(0, 0, 0),
      VecI<3>(1, 0, 0),
      VecI<3>(1, 1, 0),
      VecI<3>(0, 1, 0),
      VecI<3>(0, 0, 1),
      VecI<3>(1, 0, 1),
      VecI<3>(1, 1, 1),
      VecI<3>(0, 1, 1)};
  const array<VecI<4>, 12> LOCAL_EDGE_COORDS{
      VC::vecn::concat(LOCAL_POS_COORDS[0], 0u),
      VC::vecn::concat(LOCAL_POS_COORDS[1], 1u),
      VC::vecn::concat(LOCAL_POS_COORDS[3], 0u),
      VC::vecn::concat(LOCAL_POS_COORDS[0], 1u),
      VC::vecn::concat(LOCAL_POS_COORDS[4], 0u),
      VC::vecn::concat(LOCAL_POS_COORDS[5], 1u),
      VC::vecn::concat(LOCAL_POS_COORDS[7], 0u),
      VC::vecn::concat(LOCAL_POS_COORDS[4], 1u),
      VC::vecn::concat(LOCAL_POS_COORDS[0], 2u),
      VC::vecn::concat(LOCAL_POS_COORDS[1], 2u),
      VC::vecn::concat(LOCAL_POS_COORDS[2], 2u),
      VC::vecn::concat(LOCAL_POS_COORDS[3], 2u)};
  VecI<3> globalCornerCoord(const VecI<3> &cell_coord, uint lookup_corner_idx)
  {
    return cell_coord + LOCAL_POS_COORDS[lookup_corner_idx];
  }
  VecI<4> globalEdgeCoord(const VecI<3> &cell_coord, uint lookup_edge_idx)
  {
    return VC::vecn::concat(cell_coord, 0u) + LOCAL_EDGE_COORDS[lookup_edge_idx];
  }
  uint lookupEdgeMinCornerIndex(uint lookup_edge_idx)
  {
    // 0 1 2 3 4 5 6 7 0 1 2 3
    return lookup_edge_idx % 8;
  }
  uint lookupEdgeMaxCornerIndex(uint lookup_edge_idx)
  {
    // 1 2 3 0  5 6 7 4  4 5 6 7
    if (lookup_edge_idx < 4)
      return (lookup_edge_idx + 1) % 4;
    else if (lookup_edge_idx < 8)
      return 4 + (lookup_edge_idx + 1) % 4;
    else
      return lookup_edge_idx - 4;
  }

  //--------------------------------------------------------------------------//
  void MarchingRidges3D::findCellTriangles(const VecI<3> &cell_coord)
  {
    // find position, gradient and transverse for each cell corner
    array<VecR<3>, 8> positions;
    array<VecR<3>, 8> gradients;
    array<VecR<3>, 8> transverses;
    for (uint lookup_corner_idx = 0; lookup_corner_idx < 8; ++lookup_corner_idx)
    {
      VecR<3> pos = m_grid.cornerCoordToPos(globalCornerCoord(cell_coord, lookup_corner_idx));
      positions[lookup_corner_idx] = pos;
      gradients[lookup_corner_idx] = m_num_diff.grad(pos);
      transverses[lookup_corner_idx] = VC::vecn::normalized(ridgeinfo::stelter::direction<3>(m_num_diff, pos));
    }
    // orient all transverse directions to the average and calculate the directional first derivative (corner values)
    VecR<3> avg_transverse = calcAverageTransverseDirection(transverses);
    array<real, 8> corner_vals;
    for (uint lookup_corner_idx = 0; lookup_corner_idx < 8; ++lookup_corner_idx)
    {
      if ((transverses[lookup_corner_idx] | avg_transverse) < 0)
        transverses[lookup_corner_idx] *= -1;
      corner_vals[lookup_corner_idx] = transverses[lookup_corner_idx] | gradients[lookup_corner_idx];
    }

    // there are 256 cases with positive and negative corner value signs
    uint lookup_case_idx = 0;
    for (uint lookup_corner_idx = 0; lookup_corner_idx < 8; ++lookup_corner_idx)
      lookup_case_idx |= (corner_vals[lookup_corner_idx] < 0.0) << lookup_corner_idx;

    // use the case index to calculate all required intersections
    array<optional<VecR<3>>, 12> intersections;
    for (uint lookup_edge_idx = 0; lookup_edge_idx < 12; ++lookup_edge_idx)
    {
      if (EDGE_LOOKUP_TABLE[lookup_case_idx] & (1 << lookup_edge_idx))
      {
        uint min_corner = lookupEdgeMinCornerIndex(lookup_edge_idx);
        uint max_corner = lookupEdgeMaxCornerIndex(lookup_edge_idx);
        intersections[lookup_edge_idx] = calcLineIntersection(positions[min_corner],
                                                              corner_vals[min_corner],
                                                              positions[max_corner],
                                                              corner_vals[max_corner],
                                                              avg_transverse);
      }
    }

    // use the case index to lookup the correct triangulation
    const vector<uint> lookup_triangle_indices = TRIANGLE_LOOKUP_TABLE[lookup_case_idx];
    // step by step, create all triangles if the corresponding triangles exist
    for (uint lookup_tri_idx = 0; lookup_tri_idx < lookup_triangle_indices.size(); lookup_tri_idx += 3)
    {
      // obtain lookup and global edge indices and intersections
      uint lookup_edge_indices[3];
      uint glob_edge_indices[3];

      for (uint i = 0; i < 3; ++i)
      {
        lookup_edge_indices[i] = lookup_triangle_indices[lookup_tri_idx + i];
        glob_edge_indices[i] = m_grid.edgeCoordToIndex(globalEdgeCoord(cell_coord, lookup_edge_indices[i]));
      }
      // check whether all triangle vertices exist
      if (intersections[lookup_edge_indices[0]] && intersections[lookup_edge_indices[1]] && intersections[lookup_edge_indices[2]])
      {
        // insert to containers in isolation!
#pragma omp critical
        {
          array<VertIndex, 3> vert_indices;
          for (uint i = 0; i < 3; ++i)
          {
            // try to insert global edge index
            // if success: push back vertex
            // if fail:    already found an intersection
            const auto [it, success] = m_edge_to_vertex_map.insert({glob_edge_indices[i], m_vertices.size()});
            if (success)
              m_vertices.push_back(*intersections[lookup_edge_indices[i]]);
            // in any case, use it->second which is the vertex index
            vert_indices[i] = it->second;
          }
          m_triangles.push_back(vert_indices);
        }
      }
    }
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void transformMarRidgesScValues(ScalarFieldSeriesBase<3> &sc_series_base,
                                  const string &input_dir,
                                  const string &output_dir)
  {
    namespace fs = filesystem;
    fs::create_directories(output_dir);

    using UIntBlender_t = u_int32_t;
    using FloatBlender_t = float;
    using VecBlender_t = VC::vecn::vecn<3, FloatBlender_t>;
    using TriIndicesBlender_t = VC::vecn::vecn<3, UIntBlender_t>;

    real global_min = Globals::INF;
    real global_max = -Globals::INF;
    real global_sum = 0.0;
    uint global_num_particles = 0;

    for (uint step = 0; step < sc_series_base.steppedTime().steps; ++step)
    {
      sc_series_base.loadStep(step);

      // prepare files
      string in_path = input_dir + "/" + to_string(step) + ".data";
      string out_path = output_dir + "/" + to_string(step) + ".data";
      ifstream in(in_path, ios::binary);
      ofstream out(out_path, ios::binary);
      if (!in)
        throw runtime_error("transformMarRidgesScValues: could not read from file '" + in_path + "'");
      if (!out)
        throw runtime_error("transformMarRidgesScValues: could not write to file '" + out_path + "'");

      // read / write number of vertices and triangles
      UIntBlender_t num_vertices, num_triangles;
      in.read((char *)&num_vertices, sizeof(UIntBlender_t));
      in.read((char *)&num_triangles, sizeof(UIntBlender_t));
      out.write((char *)&num_vertices, sizeof(UIntBlender_t));
      out.write((char *)&num_triangles, sizeof(UIntBlender_t));

      // skip min and max values
      in.ignore(sizeof(FloatBlender_t) * 2);

      // read all vertices and triangle indices
      vector<VecBlender_t> vertices(num_vertices);
      vector<TriIndicesBlender_t> triangles(num_triangles);
      for (uint vert_idx = 0; vert_idx < num_vertices; ++vert_idx)
      {
        in.read((char *)&vertices[vert_idx], sizeof(VecBlender_t));
        in.ignore(sizeof(FloatBlender_t)); // ignore old scalar value
      }
      in.read((char *)&triangles[0], sizeof(TriIndicesBlender_t) * num_triangles);

      // precalculate all scalar values
      vector<FloatBlender_t> vals(num_vertices);
#pragma omp parallel for schedule(dynamic)
      for (uint vert_idx = 0; vert_idx < num_vertices; ++vert_idx)
      {
        VecR<3> real_pos(vertices[vert_idx][0], vertices[vert_idx][1], vertices[vert_idx][2]);
        vals[vert_idx] = sc_series_base.value(real_pos);
      }

      // find and write min and max values
      auto iters_minmax = minmax_element(vals.begin(), vals.end());
      FloatBlender_t val_min = iters_minmax.first == vals.end() ? -Globals::INF : *iters_minmax.first;
      FloatBlender_t val_max = iters_minmax.second == vals.end() ? Globals::INF : *iters_minmax.second;
      out.write((char *)&val_min, sizeof(FloatBlender_t));
      out.write((char *)&val_max, sizeof(FloatBlender_t));

      // write all vertex data with scalar values
      for (uint vert_idx = 0; vert_idx < num_vertices; ++vert_idx)
      {
        out.write((char *)&vertices[vert_idx], sizeof(VecBlender_t));
        out.write((char *)&vals[vert_idx], sizeof(FloatBlender_t));
      }
      // write all triangle data
      out.write((char *)&triangles[0], sizeof(TriIndicesBlender_t) * num_triangles);

      in.close();
      out.close();

      // handle statistics
      if (val_min < global_min && val_min != -Globals::INF)
        global_min = val_min;
      if (val_max > global_max && val_max != Globals::INF)
        global_max = val_max;
      global_sum += std::accumulate(vals.begin(), vals.end(), 0.0);
      global_num_particles += num_vertices;
    }
    string stats_save_path = output_dir + "/stats.txt";
    ofstream out_stats(stats_save_path);
    out_stats << "Min:  " << global_min
              << "\nMax:  " << global_max
              << "\nMean: " << (global_sum / global_num_particles);
    out_stats.close();
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//