#pragma once
//--------------------------------------------------------------------------//
#include "../scalarfield.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  /**
   * @brief A 3D example scalar field. It represents a 2D scalar field with a
   * third component for the time, i.e. it is a time-dependent scalar field. The
   * spatial domain is [0, 1] x [0, 1]. Circular ridges come closer to (but
   * never touch) the unit circle while becoming sharper the closer they are.
   * This emulates the typical behavior of FTLE ridges.
   */
  class MovingCircleRidges2D : public ScalarField<3>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    MovingCircleRidges2D();
    //--------------------------------------------------------------------------//
    MovingCircleRidges2D(const Domain<2> &domain);
    //--------------------------------------------------------------------------//
    virtual ~MovingCircleRidges2D() {}
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
  /**
   * @brief A scalar field which represents the minimum distance to the ridges
   * of MovingCircleRidges2D for any position.
   */
  class MovingCircleRidgesMinDistance2D : public ScalarField<3>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    MovingCircleRidgesMinDistance2D();
    //--------------------------------------------------------------------------//
    MovingCircleRidgesMinDistance2D(const Domain<2> &domain);
    //--------------------------------------------------------------------------//
    virtual ~MovingCircleRidgesMinDistance2D() {}
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  /**
   * @brief A 4D example scalar field. It represents a 3D scalar field with a
   * fourth component for the time, i.e. it is a time-dependent scalar field. The
   * spatial domain is [0, 1] x [0, 1] x [0, 1]. Similar to
   * MovingCircleRidges2D where values along the third dimension are
   * constant.
   */
  class MovingCircleRidges3D : public ScalarField<4>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    MovingCircleRidges3D();
    //--------------------------------------------------------------------------//
    MovingCircleRidges3D(const Domain<3> &domain);
    //--------------------------------------------------------------------------//
    virtual ~MovingCircleRidges3D() {}
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
  /**
   * @brief A scalar field which represents the minimum distance to the ridges
   * of MovingCircleRidges3D for any position.
   */
  class MovingCircleRidgesMinDistance3D : public ScalarField<4>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    MovingCircleRidgesMinDistance3D();
    //--------------------------------------------------------------------------//
    MovingCircleRidgesMinDistance3D(const Domain<3> &domain);
    //--------------------------------------------------------------------------//
    virtual ~MovingCircleRidgesMinDistance3D() {}
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//