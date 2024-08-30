#include "movingcircleridges.hpp"
//--------------------------------------------------------------------------//
using namespace std;
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  real waveMu(real t, real start_val)
  {
    return 1.0 - exp(-start_val - t); // for increasing t mu comes exponentially closer to unit circle
  }
  //--------------------------------------------------------------------------//
  real waveSigma(real t, real start_val)
  {
    return (1.01 - waveMu(t, start_val)) * 0.1; // for increasing t sigma decreases exponentially
  }
  //--------------------------------------------------------------------------//
  real waveValue(real distance, real t, real start_val)
  {
    const real mu = waveMu(t, start_val);
    if (mu <= 0.0)
      return 0.0;
    const real sigma = waveSigma(t, start_val);
    return exp(-0.5 * pow((distance - mu) / sigma, 2)) * sqrt(mu);
  }
  //--------------------------------------------------------------------------//
  constexpr real WAVE_START_0 = 0.1;          // start parameter for first wave
  constexpr real WAVE_START_DIFFERENCE = 1.0; // difference between two waves
  //--------------------------------------------------------------------------//
  real timeDepMovingCircleRidges2D(const VecR<2> &pos, real t)
  {
    const real distance = VC::vecn::norm(VecR<2>(pos[0], pos[1]));
    const real start_0 = WAVE_START_0;
    const real start_dif = WAVE_START_DIFFERENCE;

    // in case of larger distance than mu of first wave: full influence of this wave
    if (distance >= waveMu(t, start_0))
      return waveValue(distance, t, start_0);

    // search correct interval
    real start_left = start_0 - start_dif;
    while (distance < waveMu(t, start_left))
      start_left -= start_dif;
    real start_right = start_left + start_dif;
    // linear combination of left and right wave
    real fac_right = (distance - waveMu(t, start_left)) / (waveMu(t, start_right) - waveMu(t, start_left));
    return (1.0 - fac_right) * waveValue(distance, t, start_left) +
           fac_right * waveValue(distance, t, start_right);
  };
  //--------------------------------------------------------------------------//
  real timeDepMovingCircleRidgesMinDistance2D(const VecR<2> &pos, real t)
  {
    const real distance = VC::vecn::norm(VecR<2>(pos[0], pos[1]));

    real start = WAVE_START_0;
    // search as long as distance to current waveMu becomes smaller
    real min_dist = Globals::INF, cur_dist = Globals::INF;
    do
    {
      min_dist = cur_dist;
      cur_dist = abs(distance - waveMu(t, start));
      start -= WAVE_START_DIFFERENCE;
    } while (cur_dist < min_dist);
    return min_dist;
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  MovingCircleRidges2D::MovingCircleRidges2D()
      : MovingCircleRidges2D({Range(0, 1), Range(0, 1)}) {}
  //--------------------------------------------------------------------------//
  MovingCircleRidges2D::MovingCircleRidges2D(const Domain<2> &domain)
      : ScalarField<3>([](const VecR<3> &pos)
                       { return timeDepMovingCircleRidges2D({pos[0], pos[1]}, pos[2]); },
                       {domain[0], domain[1], Range()}) {}
  //--------------------------------------------------------------------------//
  MovingCircleRidgesMinDistance2D::MovingCircleRidgesMinDistance2D()
      : MovingCircleRidgesMinDistance2D({Range(0, 1), Range(0, 1)}) {}
  //--------------------------------------------------------------------------//
  MovingCircleRidgesMinDistance2D::MovingCircleRidgesMinDistance2D(const Domain<2> &domain)
      : ScalarField<3>([](const VecR<3> &pos)
                       { return timeDepMovingCircleRidgesMinDistance2D({pos[0], pos[1]}, pos[2]); },
                       {domain[0], domain[1], Range()}) {}
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  MovingCircleRidges3D::MovingCircleRidges3D()
      : MovingCircleRidges3D({Range(0, 1), Range(0, 1), Range(0, 1)}) {}
  //--------------------------------------------------------------------------//
  MovingCircleRidges3D::MovingCircleRidges3D(const Domain<3> &domain)
      : ScalarField<4>([](const VecR<4> &pos)
                       { return timeDepMovingCircleRidges2D({pos[0], pos[1]}, pos[3]); },
                       {domain[0], domain[1], domain[2], Range()}) {}
  //--------------------------------------------------------------------------//
  MovingCircleRidgesMinDistance3D::MovingCircleRidgesMinDistance3D()
      : MovingCircleRidgesMinDistance3D({Range(0, 1), Range(0, 1), Range(0, 1)}) {}
  //--------------------------------------------------------------------------//
  MovingCircleRidgesMinDistance3D::MovingCircleRidgesMinDistance3D(const Domain<3> &domain)
      : ScalarField<4>([](const VecR<4> &pos)
                       { return timeDepMovingCircleRidgesMinDistance2D({pos[0], pos[1]}, pos[3]); },
                       {domain[0], domain[1], domain[2], Range()}) {}
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//