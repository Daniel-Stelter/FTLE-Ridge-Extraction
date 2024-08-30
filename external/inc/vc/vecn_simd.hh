#ifndef VC_VECN_VECN_SIMD_HH
#define VC_VECN_VECN_SIMD_HH
//=============================================================================
namespace VC {
namespace vecn {
//=============================================================================

# ifdef DOXYGEN_SKIP

/// defines a vector type if available
template <int D,typename T,bool ENABLE=true>
struct simd_vecn {};

# else

//-----------------------------------------------------------------------------

// TODO: __has_builtin(__bulitin_convertvector) clang only,
//       applies only to few types/combinations

// TODO: __builtin_shuffle (gcc), __builtin_shufflevector (clang)

template <int D,typename T,bool ENABLE=true>
struct simd_vecn {
  using value_type = T;
  using type = T[D];
  static constexpr bool is_simd() { return false; }

  template <typename F>
  auto map(F _f) const {
    using R=decltype(_f(data[0]));
    simd_vecn<D,R,ENABLE> z;
    for (size_t i=0;i<D;++i)
      z.data[i]=_f(data[i]);
    return z;
  }

  template <typename F>
  auto zip(F _f,const simd_vecn<D,T,ENABLE>& _y) const {
    using R=decltype(_f(data[0],data[0]));
    simd_vecn<D,R,ENABLE> z;
    for (size_t i=0;i<D;++i)
      z.data[i]=_f(data[i],_y.data[i]);
    return z;
  }

  template <typename F>
  auto zip(F _f,
           const simd_vecn<D,T,ENABLE>& _y,
           const simd_vecn<D,T,ENABLE>& _a) const {
    using R=decltype(_f(data[0],data[0],data[0]));
    simd_vecn<D,R,ENABLE> z;
    for (size_t i=0;i<D;++i)
      z.data[i]=_f(data[i],_y.data[i],_a.data[i]);
    return z;
  }

  template <typename F,typename Op>
  auto reduce(F _f,Op _op,const simd_vecn<D,T,ENABLE>& _y) const {
    auto r=_f(data[0],_y.data[0]);
    for (size_t i=1;i<D;++i)
      r=_op(r,_f(data[i],_y.data[i]));
    return r;
  }

  auto vbroadcast(T _a) const { return _a; }

  type data;
};

template <int D,typename T,bool S,typename TI,bool SI>
auto shuffle(const simd_vecn<D,T,S>& _x,const simd_vecn<D,TI,SI>& _mask)  {
  static_assert(std::is_integral<TI>(),"mask must be an integer vector");
  simd_vecn<D,T,S> z;
  for (size_t i=0;i<D;++i)
    z.data[i]=_x.data[_mask.data[i]];
  return z;
}

// Note: {_a} yields (_a,0,...); {0}+_a does not work for clang

# define VBROADCAST2(x) type {x,x}
# define VBROADCAST4(x) type {x,x,x,x}
# define VBROADCAST8(x) type {x,x,x,x, x,x,x,x}
# define VBROADCAST16(x) type {x,x,x,x, x,x,x,x, x,x,x,x, x,x,x,x}

# define SPECIALIZE_SIMD(D,T)                                           \
  template <>                                                           \
  struct simd_vecn<D,T,true> {                                          \
    using value_type = T;                                               \
    using type = T __attribute__ ((__vector_size__ (sizeof(T)*D)));     \
    static constexpr bool is_simd() { return true; }                    \
                                                                        \
    type data;                                                          \
                                                                        \
    template <typename F>                                               \
    auto map(F _f) const { return simd_vecn<D,T,true> { _f(data) }; }   \
    template <typename F>                                               \
    auto zip(F _f,const simd_vecn<D,T,true>& _y) const {                \
      return simd_vecn<D,T,true> { _f(data,_y.data) };                  \
    }                                                                   \
    template <typename F>                                               \
    auto zip(F _f,                                                      \
             const simd_vecn<D,T,true>& _y,                             \
             const simd_vecn<D,T,true>& _a) const {                     \
      return simd_vecn<D,T,true> { _f(data,_y.data,_a.data) };          \
    }                                                                   \
    template <typename F,typename Op>                                   \
    auto reduce(F _f,Op _op,const simd_vecn<D,T,true>& _y) const {      \
      auto z=_f(data,_y.data);                                          \
      auto r=z[0];                                                      \
      for (size_t i=1;i<D;++i)                                          \
        r=_op(r,z[i]);                                                  \
      return r;                                                         \
    }                                                                   \
    type vbroadcast(T _a) const { return VBROADCAST ## D(_a); }         \
  }                                                                     \

# define SPECIALIZE_SHUFFLE(D,T,I)                                      \
  template <>                                                           \
  inline auto shuffle(const simd_vecn<D,T,true>& _x,                    \
               const simd_vecn<D,I,true>& _mask) {                      \
    return simd_vecn<D,T,true> { BULTIIN_SHUFFLE(_x.data,_mask.data) }; \
  }                                                                     \

# if !defined(VC_VECN_NO_SIMD)

# if (defined(__GNUC__) || defined(__clang__)) && !defined(__CUDACC__)

# if defined(__clang__)
#  define BULTIIN_SHUFFLE(x,m) __builtin_shufflevector(x,m)
# else
#  define BULTIIN_SHUFFLE(x,m) __builtin_shuffle(x,m)
# endif

#  if defined(__SSE__)

// template <>
// auto simd_vecn<2,float,true>::vbroadcast(float _a) const {
//   return simd_vecn<2,float,true>::type {_a,_a};
// }

SPECIALIZE_SIMD(2,float);
SPECIALIZE_SIMD(4,float);
SPECIALIZE_SIMD(2,std::int32_t);
SPECIALIZE_SIMD(4,std::int32_t);
SPECIALIZE_SIMD(2,double);

SPECIALIZE_SHUFFLE(4,float,std::int32_t)
SPECIALIZE_SHUFFLE(4,int32_t,std::int32_t)

SPECIALIZE_SIMD(4,std::int8_t);
SPECIALIZE_SIMD(4,std::uint8_t);
SPECIALIZE_SIMD(8,std::int8_t);
SPECIALIZE_SIMD(8,std::uint8_t);
SPECIALIZE_SIMD(16,std::int8_t);
SPECIALIZE_SIMD(16,std::uint8_t);

// https://godbolt.org/ suggests that st of the following are useless?!

SPECIALIZE_SHUFFLE(4,std::int8_t,std::int8_t)
SPECIALIZE_SHUFFLE(4,std::uint8_t,std::int8_t)
SPECIALIZE_SHUFFLE(8,std::int8_t,std::int8_t)
SPECIALIZE_SHUFFLE(8,std::uint8_t,std::int8_t)
SPECIALIZE_SHUFFLE(16,std::int8_t,std::int8_t)
SPECIALIZE_SHUFFLE(16,std::uint8_t,std::int8_t)

# if defined(__clang__)
  SPECIALIZE_SHUFFLE(2,double,std::int32_t)
#endif

#  endif

#  if defined(__AVX__)

SPECIALIZE_SIMD(8,float);
SPECIALIZE_SIMD(8,std::int32_t);
SPECIALIZE_SIMD(4,double);

SPECIALIZE_SHUFFLE(8,float,std::int32_t)
SPECIALIZE_SHUFFLE(8,int32_t,std::int32_t)

# if defined(__clang__)
  SPECIALIZE_SHUFFLE(4,double,std::int32_t)
#endif

#  endif

# endif

# undef VBROADCAST2
# undef VBROADCAST4
# undef VBROADCAST8

# undef BUILTIN_SHUFFLE

# endif // VC_VECN_NO_SIMD

# endif // DOXYGEN_SKIP

//=============================================================================
} //namespace vecn
} //namespace VC
//=============================================================================
#endif // VC_VECN_VECN_SIMD_HH defined
