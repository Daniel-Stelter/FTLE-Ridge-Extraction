#ifndef VC_VECN_VECN_SPECIALIZATION_HH
#define VC_VECN_VECN_SPECIALIZATION_HH
//-----------------------------------------------------------------------------
# ifndef DOXYGEN_SKIP
//=============================================================================
namespace VC {
namespace vecn {
//=============================================================================
#  ifndef VC_VECN_NO_SPECIALIZATION
//-----------------------------------------------------------------------------
template <typename T>
struct inc_t<2,T> {
  auto operator()(const vecn<2,T>& _ii,const vecn<2,T>& _d) {
    if (_ii[0]+1<_d[0])
      return vecn<2,T>(_ii[0]+1,_ii[1]);
    else
      return vecn<2,T>(T(0),_ii[1]+1);
  }
};

template <typename T>
struct inc_t<3,T> {
  auto operator()(const vecn<3,T>& _ii,const vecn<3,T>& _d) {
    if (_ii[0]+1<_d[0])
      return vecn<3,T>(_ii[0]+1,_ii[1],_ii[2]);
    else if (_ii[1]+1<_d[1])
      return vecn<3,T>(T(0),_ii[1]+1,_ii[2]);
    else
      return  vecn<3,T>(T(0),T(0),_ii[2]+1);
  }
};

template <typename T>
struct inc_t<4,T> {
  auto operator()(const vecn<4,T>& _ii,const vecn<4,T>& _d) {
    if (_ii[0]+1<_d[0])
      return vecn<4,T>(_ii[0]+1,_ii[1],_ii[2],_ii[3]);
    else if (_ii[1]+1<_d[1])
      return vecn<4,T>(T(0),_ii[1]+1,_ii[2],_ii[3]);
    else if (_ii[2]+1<_d[2])
      return vecn<4,T>(T(0),T(0),_ii[2]+1,_ii[3]);
    else
      return vecn<4,T>(T(0),T(0),T(0),_ii[3]+1);
  }
};

//-----------------------------------------------------------------------------
#  endif // VC_VECN_NO_SPECIALIZATION
//=============================================================================
} //namespace vecn
} //namespace VC
//=============================================================================
# endif // DOXYGEN_SKIP
//-----------------------------------------------------------------------------
#endif // VC_VECN_VECN_SPECIALIZATION_HH defined
