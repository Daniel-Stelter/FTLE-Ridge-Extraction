//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MVFIELDS_GRID_D1_HH
#define VC_MVFIELDS_GRID_D1_HH

#include <cassert>
#include <climits> // _MSC_VER
#include <cmath>
#include <algorithm>

namespace VC {
namespace mvfields {

/// univariate reconstruction [\ref vc_mvfields_grid] \ingroup vc_mvfields_grid
namespace d1 {
//-----------------------------------------------------------------------------

/// number of dimensions
const int DIMS=1;

/// vector type (a scalar for 1d)
template <typename T> using vec_t = T;

typedef int index_t;        //!< linear index type

typedef vec_t<int> mindex_t; //!< multi-dimensional index type

typedef mindex_t msize_t;   //!< multi-dimensional size type

//-----------------------------------------------------------------------------

/// get product of indices \ingroup vc_mvfields_grid
inline static index_t prod(const mindex_t& _mi) {
  return _mi;
}
/// zero index? \ingroup vc_mvfields_grid
inline static bool is_zero(const mindex_t& _mi) {
  return _mi==0;
}

//-----------------------------------------------------------------------------

/** \class marray
    \brief Describe 1d array.
    \ingroup vc_mvfields_grid

    \tparam T `value_type`

    This class provides
    \arg indexed access to data,
    \arg copying of blocks.

    Note that marray does *not* manage memory for data!
 */
template <typename T>
class marray {
public:
  typedef T value_type; //!< value type

  /// define array
  marray(T* _data,msize_t _size)
    : m_data(_data), m_size(_size) {
  }
  /// copy array (does _not_ copy data)
  marray(const marray& _other)
    : m_data(_other.m_data), m_size(_other.m_size) {
  }

  /// number of dimensions
  static constexpr unsigned ndims() { return 1; }
  /// array size per dimension (trivial in 1d)
  msize_t size() const { return m_size; }
  /// total number of elements (same as size() 1d)
  unsigned numel() const { return unsigned(m_size); }

  /// get index from multi-index (trivial in 1d)
  index_t index(mindex_t _idx) const { return _idx; }
  /// get multi-index from linear index
  mindex_t mindex(index_t _i) const { return _i; }

  /// get data
  T* data() const { return m_data; }
  /// same as data()
  T* begin() const { return m_data; }
  /// end of data
  T* end() const { return m_data+numel(); }

  /// copy a block of `_size` starting at `_idx` to `dst`
  template <typename TT>
  void cp_block(TT* _dst,mindex_t _idx,msize_t _size) const {
    for (int i=0;i<_size;++i)
      _dst[i]=TT(m_data[_idx+i]);
  }

  /** Get "overlap" of block  with data bounds.

      The method assumes that size() `>_size`, i.e., there may be _either_
      left overlap _or_ right overlap.

      \param _idx index
      \param _size block size
      \return number of overlapping elements, sign denotes
      "left" and "right" boundary. A return value `0` indicates that
      the requested block is inside the domain (no overlap).
   */
  msize_t overlap(mindex_t _idx,msize_t _size) const {
    assert(_size<=m_size);

    if (_idx<0)                  return _idx;
    else if (_idx+_size>=m_size) return (_idx+_size)-m_size;
    return 0;
  }

protected:
  T*      m_data; //!< pointer to data
  msize_t m_size; //!< size of data
};

//-----------------------------------------------------------------------------

/** \class grid
    \brief Describe 1d regular grid.
    \ingroup vc_mvfields_grid

    \tparam R real numbers
    \tparam T `value_type`

    This class provides in addition to marray
    \arg description of the rectangular domain,
    \arg mapping between indices and domain points.

    There are two types of coordinates:
    \arg The domain is defined in _world_ coordinates with cell size h().
    \arg Reconstruction is is done in _local_ coordinates with cell size 1.

    Note that grid does *not* manage memory for data!
 */
template <typename R,typename T>
class grid : public marray<T> {
public:
  typedef R real_t;          //!< real
  typedef R coord_t;         //!< coordinate (in 1d same as real_t)
  typedef T value_type;      //!< value type
  typedef marray<T> array_t; //!< array type

  /// undefined grid
  grid() : marray<T>(nullptr,mindex_t(0)) {}

  /// define grid from `_data` and `_size` in `[0,1]`
  grid(T* _data,msize_t _size)
    : marray<T>(_data,_size) {
    assert(_size>1);
    set_origin(coord_t(0));
    set_extent(coord_t(1));
  }
  /// copy grid (does _not_ copy data!)
  grid(const grid& _other)
    : marray<T>(_other),
      m_origin(_other.m_origin), m_extent(_other.m_extent), m_h(_other.m_h) {
  }

  /// origin of domain
  const coord_t& origin() const { return this->m_origin; }
  /// extent of rectangular domain
  const coord_t& extent() const { return this->m_extent; }
  /// cell size
  const coord_t& h() const { return this->m_h; }

  /// set origin
  void set_origin(const coord_t& _origin) { this->m_origin=_origin; }
  /// set extent of rectangular domain
  void set_extent(const coord_t& _extent) {
    assert(_extent>real_t(0));
    this->m_extent=_extent;
    this->m_h=this->m_extent/coord_t(this->size()-1);
  }
  /// set bounding box
  template <typename RR>
  void set_bbox(const RR* _bbox) {
    R bbox[2];
    bbox[0]=std::min(_bbox[0],_bbox[1]);
    bbox[1]=std::max(_bbox[0],_bbox[1]);
    set_origin(coord_t(bbox[0]));
    set_extent(coord_t(bbox[1])-this->m_origin);
  }
  /// Is `_x` inside?
  bool inside(const coord_t& _x) const {
    return
      this->m_origin<=_x && _x<=this->m_origin+this->m_extent;
  }

  /// map local coordinates `_x` to world coordinates
  coord_t to_world(const coord_t& _x) const {
    return this->m_origin+_x*m_h;
  }
  /// map world coordinates `_x` to local coordinates
  coord_t to_local(const coord_t& _x) const {
    return (_x-this->m_origin)/m_h;
  }
  /// truncate local coordinates to integers
  static coord_t floor(const coord_t& _x) {
    return coord_t(std::floor(_x));
  }

protected:
  coord_t m_origin; //!< origin
  coord_t m_extent; //!< extent of rectangular domain
  coord_t m_h;      //!< cell size
};

//-----------------------------------------------------------------------------
/// reconstruction kernels \ingroup vc_mvfields_grid
namespace kernel {
//-----------------------------------------------------------------------------

/** \class linear_interpolatory_c0
    \brief Linear interpolatory C^0 reconstruction kernel.
    \ingroup vc_mvfields_grid

    \tparam GRID type must represent a grid
    \tparam T `value_type` used for reconstruction (default: same as in grid)
    \tparam R real numbers (default: same is in grid)

    Note that derivatives evaluated fx_local() are generally
    _discontinuous_! (fx_local() does _not_ apply a finite difference
    scheme over multiple grid cells!)
 */
template <typename GRID,
          typename T=typename GRID::value_type,
          typename R=typename GRID::real_t>
class linear_interpolatory_c0 {
public:
  typedef GRID grid_t;      //!< the grid
  typedef T value_type;     //!< value type
  typedef R real_t;         //!< real numbers

# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t; //!< coordinates (real numbers in 1d)
# else
  typedef R coord_t;
# endif

  enum { SUPPORT_BUFFER_SIZE=2 }; // _MSC_VER fix: should be support_buffer_size()

  /// size of support
  static constexpr msize_t support_size() { return 2; }
  /// dimensions of of support (same as support_size() for 1d)
  static constexpr size_t support_buffer_size() { return 2; }
  /// offset to index returned by get_cell()
  static constexpr int offset() { return 0; }

  /** Copy support_size() coefficients from `_grid` at `_idx`.
      Eventually apply boundary conditions.
      \param _grid the grid: data source
      \param _idx source index
      \param _c destination buffer of size support_buffer_size()
      \return `false` if `_idx` is (too far) outside grid
   */
  static bool get_coeffs(const grid_t& _grid,mindex_t _idx,value_type* _c) {

    if (!is_zero(_grid.overlap(_idx,support_size())))
      return false;

    _grid.cp_block(_c,_idx,support_size());

    return true;
  }

  /** Evaluate function in local coordinates `_s`.

      \param _s local coordinates as from get_cell(), i.e., they are
      truncated such that each coordinate is in `[0,1]`
      \param _c coefficients as from get_coeffs, they define the function
      \return function value
   */
  static value_type f_local(real_t _s,const value_type* _c) {
    return _c[0]*(real_t(1)-_s)+_c[1]*_s;
  }

  /// evaluate first derivative to f_local() (independent of grid size)
  static value_type fx_local(real_t /*_s*/,const value_type* _c) {
    return _c[1]-_c[0];
  }


  /** Get cell index and local coordinates.
      \param _grid the grid
      \param _x domain point in world coordinates
      \return grid index and local coordinate of `_x`
   */
  static std::pair<mindex_t,coord_t>
  get_cell(const grid_t& _grid,const typename grid_t::coord_t& _x) {

    auto x =_grid.to_local(_x);
    auto xk=_grid.floor(x);

    real_t   s=x-xk;
    mindex_t k=mindex_t(xk);
    k+=offset();

    if (k>=1 && s==real_t(0)) {
      // prevent invalid evaluation at right boundary
      --k;
      s=real_t(1);
    }

    return std::make_pair(k,s);
  }
};

//-----------------------------------------------------------------------------

/** \class cubic_convolution_interpolatory_c1
    \brief Cubic convolution interpolatory C^1 reconstruction kernel.
    \ingroup vc_mvfields_grid

    \tparam GRID type must represent a grid
    \tparam T `value_type` used for reconstruction (default: same as in grid)
    \tparam R real numbers (default: same is in grid)

    The reconstruction is based on

    <pre>
    Robert G. Keys.
    Cubic Convolution Interpolation for Digital Image Processing.
    IEEE Trans. on Acoustics, Speech, and Signal Processing, 29(6), 1981
    </pre>

    The kernel reproduces quadratic functions and provides a cubic C^1
    continuous reconstruction that interpolates data at grid
    points. The implementation provides appropriate boundary
    conditions.
 */
template <typename GRID,
          typename T=typename GRID::value_type,
          typename R=typename GRID::real_t>
class cubic_convolution_interpolatory_c1 {
public:
  typedef GRID grid_t;  //!< the grid
  typedef T value_type; //!< value type
  typedef R real_t;     //!< real numbers

# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t; //!< coordinates
# else
  typedef real_t coord_t;
# endif

  enum { SUPPORT_BUFFER_SIZE=4 }; // _MSC_VER fix: should be support_buffer_size()

  /// size of support
  static constexpr msize_t support_size() { return 4; }
  /// dimensions of of support (same as support_size() for 1d)
  static constexpr size_t support_buffer_size() { return 4; }
  /// offset to index returned by get_cell()
  static constexpr size_t offset() { return -1; }

  /** Enforce "left" boundary conditions.
      \param _c coefficients where `_c[0]` is to be determined.
      \param _inc increment to next coefficient
   */
  static void set_left_bcond(value_type* _c,msize_t _inc=1) {
    _c[0]=(_c[_inc]-_c[2*_inc])*real_t(3)+_c[3*_inc];
  }
  /** Enforce "right" boundary conditions.
      \param _c coefficients where `_c[3*_inc]` is to be determined.
      \param _inc increment to next coefficient
   */
  static void set_right_bcond(value_type* _c,msize_t _inc=1) {
    _c[3*_inc]=_c[0]-(_c[_inc]-_c[2*_inc])*real_t(3);
  }

  /** Copy support_size() coefficients from `grid` at `_idx`.
      Eventually apply boundary conditions.
      \param _grid the grid: data source
      \param _idx source index
      \param _c destination buffer of size support_buffer_size()
      \return `false` if `_idx` is (too far) outside grid
   */
  static bool get_coeffs(const grid_t& _grid,mindex_t _idx,value_type* _c) {

    mindex_t overlap=_grid.overlap(_idx,support_size());

    if (is_zero(overlap))
      _grid.cp_block(_c,_idx,support_size());
    else {
      // Get overlap, adjust origin and size, and copy
      if (overlap<0) {
        if (-overlap>1) return false;
        _grid.cp_block(_c+1,_idx+1,3);
        set_left_bcond(_c);
      }
      else if (overlap>0) {
        if (+overlap>1) return false;
        _grid.cp_block(_c,_idx,3);
        set_right_bcond(_c);
      }
    }
    return true;
  }

  /** Evaluate function in local coordinates `_s`.

      \param _s local coordinates as from get_cell(), i.e., they are
      truncated such that each coordinate is in `[0,1]`
      \param _c coefficients as from get_coeffs, they define the function
      \return function value
   */
  static value_type f_local(real_t _s,const value_type* _c) {

    // real_t t0=((real_t(2)-_s)*_s-real_t(1))*_s;
    // real_t t1=(real_t(3)*_s-real_t(5))*_s*_s+real_t(2);
    // real_t t2=((real_t(4)-real_t(3)*_s)*_s+real_t(1))*_s;
    // real_t t3=(_s-real_t(1))*_s*_s;

    // return (_c[0]*t0+_c[1]*t1+_c[2]*t2+_c[3]*t3)*real_t(0.5);

    real_t s2=_s*_s;
    real_t s3=s2*_s;
    return ( _c[0]*(          -s3+real_t(2)*s2 -_s           )+
             _c[1]*( real_t(3)*s3-real_t(5)*s2     +real_t(2))+
             _c[2]*(-real_t(3)*s3+real_t(4)*s2 +_s           )+
             _c[3]*(           s3          -s2)
             )*real_t(0.5);
  }

  /// evaluate first derivative to f_local() (independent of grid size)
  static value_type fx_local(real_t _s,const value_type* _c) {
    real_t s2=_s*_s;
    return (_c[0]*(-real_t(3)*s2+real_t(4) *_s -real_t(1))+
            _c[1]*( real_t(9)*s2-real_t(10)*_s)+
            _c[2]*(-real_t(9)*s2+real_t(8) *_s +real_t(1))+
            _c[3]*( real_t(3)*s2-real_t(2) *_s)
           )*real_t(0.5);
  }

  /// evaluate second derivative to f_local() (independent of grid size)
  static value_type fxx_local(real_t _s,const value_type* _c) {
    return (_c[0]*(-real_t(6)*_s +real_t(4))+
            _c[1]*( real_t(18)*_s-real_t(10))+
            _c[2]*(-real_t(18)*_s+real_t(8))+
            _c[3]*( real_t(6)*_s -real_t(2))
           )*real_t(0.5);
  }

  /** Get cell index and local coordinates.
      \param _grid the grid
      \param _x domain point in world coordinates
      \return grid index and local coordinate of `_x`
   */
  static std::pair<mindex_t,real_t>
  get_cell(const grid_t& _grid,const typename grid_t::coord_t& _x) {

    auto x =_grid.to_local(_x);
    auto xk=_grid.floor(x);

    real_t   s=x-xk;
    mindex_t k=mindex_t(xk)+offset();

    if (k>=1 && s==real_t(0)) {
      // prevent invalid evaluation at right boundary
      --k;
      s=real_t(1);
    }

    return std::make_pair(k,s);
  }
};

//-----------------------------------------------------------------------------
} // namespace kernel
//-----------------------------------------------------------------------------

/** \class evaluator
    \brief Evaluate function defined by gridded values.
    \ingroup vc_mvfields_grid

    \tparam KERNEL type must represent a reconstruction kernel such as
    VC::mvfields::d1::kernel::cubic_convolution_interpolatory_c1.

    Note that the evaluator has a state (in contrast to kernels),
    i.e., multiple instances are required for parallel evaluation!

    \todo Check spurious compiler warning `maybe-uninitialized` (`g++`
    with optimization only) due to uninitialized coefficients `c`.
 */
template <typename KERNEL>
class evaluator {
public:
  typedef KERNEL kernel_t; //!< the kernel
  typedef typename kernel_t::grid_t grid_t; //!< the grid
  typedef typename kernel_t::value_type value_type; //!< value type
  typedef typename kernel_t::real_t real_t; //!< real numbers
  typedef typename kernel_t::coord_t coord_t; //!< coordinate


# if !defined(_MSC_VER)  /// tagging cell invalid forces loadc() to call `kernel_t::get_cell()`
  static const index_t INVALIDCELL=std::numeric_limits<index_t>::min()+1;
# else
  enum { INVALIDCELL=INT_MIN+1 };
# endif

  /// create evaluator for `_grid`
  evaluator(const grid_t& _grid) : grid(_grid), cell(INVALIDCELL) {}

  /** Load coefficients that are required to evaluate at `_x` to c.
      \param _x domain point for evaluation
      \return `(false,???)` is `_x` is not in domain, `(true,y)` otherwise
      with local coordinate `y` (as used by `kernel_t::f_local()`).
   */
  std::pair<bool,coord_t> loadc(const typename grid_t::coord_t& _x) {
    auto p=kernel_t::get_cell(grid,_x);
    index_t newcell=grid.index(p.first);
    if (newcell!=cell) {
      if (!kernel_t::get_coeffs(grid,p.first,c))
        return std::make_pair(false,p.second);
      this->cell=newcell;
      assert(newcell!=INVALIDCELL);
    }
    return std::make_pair(true,p.second);
  }

  /** Evaluate function value `_f` at `_x`.
      \param[out] _f value on success
      \param _x domain point
      \return `true` is `_x` is in domain (success)
   */
  bool value(value_type& _f,const typename grid_t::coord_t& _x) {
    auto p=loadc(_x);
    if (p.first)
      _f=kernel_t::f_local(p.second,c);
    return p.first;
  }
  /** Evaluate function gradient `_fx` at `_x`.
      \param[out] _fx gradient on success
      \param _x domain point
      \return `true` is `_x` is in domain (success)
   */
  bool gradient(value_type* _fx,const typename grid_t::coord_t& _x) {
    auto p=loadc(_x);
    if (p.first)
      _fx[0]=kernel_t::fx_local(p.second,c)/grid.h();
    return p.first;
  }

  value_type    c[kernel_t::SUPPORT_BUFFER_SIZE]; //!< "local" coefficients
  const grid_t& grid; //!< the grid
  index_t       cell; //!< current (linear) cell index
};

# ifndef _MSC_VER
#  undef constexpr // TODO: remove this ASAP!
# endif

//-----------------------------------------------------------------------------
} // namespace d1
} // namespace mvfields
} // namespace VC

#endif // VC_MVFIELDS_GRID_D1_HH
