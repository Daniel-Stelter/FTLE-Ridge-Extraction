//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MVFIELDS_GRID_D2_HH
#define VC_MVFIELDS_GRID_D2_HH

#include <cassert>
#include <climits> // _MSC_VER

#include <algorithm>
#include <iostream>
#include <utility>

#include "VecN.hh"

#include "grid_d1.hh"

namespace VC {
namespace mvfields {

/// bivariate reconstruction [\ref vc_mvfields_grid] \ingroup vc_mvfields_grid
namespace d2 {
//-----------------------------------------------------------------------------

/// number of dimensions \ingroup vc_mvfields_grid
const int DIMS=2;

/// vector type for indices and coordinates \ingroup vc_mvfields_grid
# if !defined(_MSC_VER)
template <typename T> using vec_t = VC::math::VecN<T,DIMS>;
# else
// TODO: remove this ASAP !!!
#  define constexpr
# endif

typedef int index_t;        //!< linear index type \ingroup vc_mvfields_grid

# if !defined(_MSC_VER)
typedef vec_t<int> mindex_t; //!< multi-dimensional index type \ingroup vc_mvfields_grid
# else
typedef VC::math::VecN<int,DIMS> mindex_t;
# endif

typedef mindex_t msize_t;   //!< multi-dimensional size type \ingroup vc_mvfields_grid

//-----------------------------------------------------------------------------

/// get product of indices \ingroup vc_mvfields_grid
inline static index_t prod(const mindex_t& _mi) {
  return _mi[0]*_mi[1];
}
/// zero index? \ingroup vc_mvfields_grid
inline static bool is_zero(const mindex_t& _mi) {
  return _mi[0]==0 && _mi[1]==0;
}

//-----------------------------------------------------------------------------

/** \class marray
    \brief Describe 2d array (first index fastest).
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
  static constexpr unsigned ndims() { return DIMS; }
  /// array size per dimension
  const msize_t& size() const { return m_size; }
  /// total number of elements
  unsigned numel() const { return prod(m_size); }

  /// get data
  T* data() const { return m_data; }
  /// same as data()
  T* begin() const { return m_data; }
  /// end of data
  T* end() const { return m_data+numel(); }

  /// get index from multi-index
  index_t index(const mindex_t& _mi) const {
    return _mi[1]*m_size[0]+_mi[0]; // x (0-th dim) fastest
  }
  /// get multi-index from linear index
  mindex_t mindex(index_t _i) const {
    return mindex_t(_i%m_size[0],_i/m_size[0]);
  }

  /// copy a block of `_size` starting at `_idx` to `dst`
  template <typename TT>
  void cp_block(TT* _dst,const mindex_t& _mi,const msize_t& _size) const {
    cp_block(_dst,index(_mi),_size);
  }
  /// copy a block of `_size` starting at `_idx` to `dst`
  template <typename TT>
  void cp_block(TT* _dst,index_t _i,const msize_t& _size) const {
    unsigned n=0;
    for (int j=0;j<_size[1];++j)
      for (int i=0;i<_size[0];++i)
        _dst[n++]=TT(m_data[_i+j*m_size[0]+i]); // let compiler do optimizations!
  }

  /** Get "overlap" of block `_mi -- _mi+_size` with data bounds.

      The method assumes that size() `>_size`, i.e., for each
      dimension there may be only one part of the requested block
      crossing the boundary.

      \param _mi "lower left" index
      \param _size block size
      \return number of overlapping elements per dimension, sign denotes
      "left" and "right" boundary. A return value `(0,0)` indicates that
      the requested block is inside the domain (no overlap).
   */
  mindex_t overlap(const mindex_t& _mi,const msize_t& _size) const {
    assert(_size[0]<=m_size[0] && _size[1]<=m_size[1]);

    mindex_t rv(0,0);

    if (_mi[0]<0) rv[0]=_mi[0];
    else if (_mi[0]+_size[0]>=m_size[0]) rv[0]=(_mi[0]+_size[0])-m_size[0];

    if (_mi[1]<0) rv[1]=_mi[1];
    else if (_mi[1]+_size[1]>=m_size[1]) rv[1]=(_mi[1]+_size[1])-m_size[1];

    return rv;
  }

  /** Same as cp_block() but provide dimensions for `_dst`.
      \param _dst destination buffer, partially filled if `_dims[1]>_size[1]`
      \param _dims for `_dst`, `_dims[1]>=_size[1]`, `_dims[0]` is ignored.
      \param _mi block origin
      \param _size block size
   */
  template <typename TT>
  void cp_block(TT* _dst,const msize_t& _dims,
                const mindex_t& _mi,const msize_t& _size) const {
    cp_block(_dst,_dims,index(_mi),_size);
  }
  /// Same as above.
  template <typename TT>
  void cp_block(TT* _dst,const msize_t& _dims,
                index_t _i,const msize_t& _size) const {
    assert(_dims[1]>=_size[1]);
    unsigned n=0;
    for (int j=0;j<_size[1];++j,n+=_dims[1])
      for (int i=0;i<_size[0];++i)
        _dst[n+i]=m_data[_i+j*m_size[0]+i]; // let compiler do optimizations!
  }

protected:
  T*      m_data; //!< pointer to data
  msize_t m_size; //!< size of data
};

//-----------------------------------------------------------------------------

/** \class grid
    \brief Describe 2d regular grid.
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
# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t;  //!< coordinates
# else
  typedef VC::math::VecN<R,DIMS> coord_t;
# endif
  typedef T value_type;      //!< value type
  typedef marray<T> array_t; //!< array type

  /// undefined grid
  grid() : marray<T>(nullptr,mindex_t(0,0)) {}

  /// define grid from `_data` and `_size` in `[0,1]`
  grid(T* _data,msize_t _size)
    : marray<T>(_data,_size) {
    assert(_size[0]>1 && _size[1]>1);
    set_origin(coord_t(0,0));
    set_extent(coord_t(1,1));
  }
  /// copy grid (does _not_ copy data!)
  grid(const grid& _other)
    : marray<T>(_other),
      m_origin(_other.m_origin), m_extent(_other.m_extent), m_h(_other.m_h) {
  }

  /// origin of domain
  const coord_t& origin() const { return this->m_origin; }
  ///  extent of rectangular domain
  const coord_t& extent() const { return this->m_extent; }
  /// cell size
  const coord_t& h() const { return this->m_h; }

  /// set origin
  void set_origin(const coord_t& _origin) { this->m_origin=_origin; }
  /// set extent of rectangular domain
  void set_extent(const coord_t& _extent) {
    assert(_extent[0]>real_t(0) && _extent[1]>real_t(0));
    this->m_extent=_extent;
    this->m_h=this->m_extent/coord_t(this->size()[0]-1,this->size()[1]-1);
  }
  /// set bounding box
  template <typename RR>
  void set_bbox(const RR* _bbox) {
    R bbox[2*DIMS];
    bbox[0]=std::min(_bbox[0],_bbox[1]);
    bbox[1]=std::max(_bbox[0],_bbox[1]);
    bbox[2]=std::min(_bbox[2],_bbox[3]);
    bbox[3]=std::max(_bbox[2],_bbox[3]);
    set_origin(coord_t(bbox[0],bbox[2]));
    set_extent(coord_t(bbox[1],bbox[3])-this->m_origin);
  }

  /// Is `_x` inside?
  bool inside(const coord_t& _x) const {
    return
      this->m_origin[0]<=_x[0] && _x[0]<=this->m_origin[0]+this->m_extent[0] &&
      this->m_origin[1]<=_x[1] && _x[1]<=this->m_origin[1]+this->m_extent[1];
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
    return coord_t(::floor(_x[0]),::floor(_x[1]));
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

    Note that derivatives evaluated fx_local() and fy_local() are generally
    _discontinuous_! (fx_local() does _not_ apply a finite difference
    scheme over multiple grid cells!)

    \sa d1::kernel::linear_interpolatory_c0
 */
template <typename GRID,
          typename T=typename GRID::value_type,
          typename R=typename GRID::real_t>
class linear_interpolatory_c0 {
public:
  typedef GRID grid_t;        //!< the grid
  typedef T value_type;       //!< value type
  typedef R real_t;           //!< real numbers
#if !defined(_MSC_VER)
  typedef vec_t<R> coord_t; //!< coordinates
# else
  typedef VC::math::VecN<R,DIMS> coord_t;
# endif
  typedef d1::kernel::linear_interpolatory_c0<grid_t,value_type,real_t> k1_t;

  enum { SUPPORT_BUFFER_SIZE=4 }; // _MSC_VER fix: should be support_buffer_size()

  /// size of support
  static msize_t support_size() { return msize_t(2,2); }
  /// dimensions of of support (same as support_size() for 1d)
  static constexpr size_t support_buffer_size() { return 4; }
  /// offset to index returned by get_cell()
  static constexpr size_t offset() { return 0; }


  /** Copy support_size() coefficients from `_grid` at `_idx`.
      Eventually apply boundary conditions.
      \param _grid the grid: data source
      \param _idx source index
      \param _c destination buffer of size support_buffer_size()
      \return `false` if `_idx` is (too far) outside grid
   */
  static bool get_coeffs(const grid_t& _grid,const mindex_t& _idx,value_type* _c) {

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
  static value_type f_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    c[0]=k1_t::f_local(_s[0],_c  );
    c[1]=k1_t::f_local(_s[0],_c+2);
    return k1_t::f_local(_s[1],c);
  }

  /// evaluate first derivative d/dx to f_local() (independent of grid size)
  static value_type fx_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    c[0]=k1_t::fx_local(_s[0],_c  );
    c[1]=k1_t::fx_local(_s[0],_c+2);
    return k1_t::f_local(_s[1],c);
  }

  /// evaluate first derivative d/dy to f_local() (independent of grid size)
  static value_type fy_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    c[0]=k1_t::f_local(_s[0],_c  );
    c[1]=k1_t::f_local(_s[0],_c+2);
    return k1_t::fx_local(_s[1],c);
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

    coord_t s=x-xk;
    mindex_t k=xk;
    k+=offset();

    for (int i=0;i<DIMS;++i) {
      if (k[i]>=1 && s[i]==real_t(0)) {
        // prevent invalid evaluation at right boundary
        --k[i];
        s[i]=real_t(1);
      }
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

    \sa d1::kernel::cubic_convolution_interpolatory_c1
 */
template <typename GRID,
          typename T=typename GRID::value_type,
          typename R=typename GRID::real_t>
class cubic_convolution_interpolatory_c1 {
public:
  typedef GRID grid_t;        //!< the grid
  typedef T value_type;       //!< value type
  typedef R real_t;           //!< real numbers
# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t; //!< coordinates
# else
  typedef VC::math::VecN<R,DIMS> coord_t;
# endif
  typedef d1::kernel::cubic_convolution_interpolatory_c1<grid_t,value_type,real_t> k1_t;


  enum { SUPPORT_BUFFER_SIZE=16 }; // _MSC_VER fix: should be support_buffer_size()

  static msize_t support_size() { return msize_t(4,4); }
  static constexpr size_t support_buffer_size() { return 16; }
  static constexpr size_t offset() { return -1; }

  static void set_left_bcond(value_type* _c,const msize_t& _inc) {
    for (int i=0;i<4;++i)
      k1_t::set_left_bcond(_c+i*_inc[0],_inc[1]);
  }
  static void set_right_bcond(value_type* _c,const msize_t& _inc) {
    for (int i=0;i<4;++i)
      k1_t::set_right_bcond(_c+i*_inc[0],_inc[1]);
  }

  /// copy support_size() coefficients from `grid` at `_idx`
  static bool get_coeffs(const grid_t& _grid,
                         const mindex_t& _mi,value_type* _c) {

    mindex_t overlap=_grid.overlap(_mi,support_size());

    if (is_zero(overlap))
      _grid.cp_block(_c,_mi,support_size());
    else {
      // copy partial block
      {
        value_type* c0=_c;
        mindex_t    mi(_mi); // new origin
        mindex_t    sz(4,4); // new size
        const int   INC[DIMS]={1,4};

        for (int k=0;k<DIMS;++k) {
           if ( overlap[k]<0)     {
             if (overlap[k]<-1) return false;
             --sz[k]; ++mi[k]; c0+=INC[k];
           }
           else if (overlap[k]>0) {
             if (overlap[k]>1)  return false;
             --sz[k];
           }
        }

        _grid.cp_block(c0,msize_t(1,4),mi,sz);
      }
      // compute missing coefficients

      if (overlap[0]<0)
        set_left_bcond(_c,msize_t(4,1)); // slice i
      else if (overlap[0]>0)
        set_right_bcond(_c,msize_t(4,1));

      if (overlap[1]<0)
        set_left_bcond(_c,msize_t(1,4)); // slice j
      else if (overlap[1]>0)
        set_right_bcond(_c,msize_t(1,4));

    }
    return true;
  }

  /// evaluate function in local coordinates `_s`
  static value_type f_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];
    c[0]=k1_t::f_local(_s[0],_c  );
    c[1]=k1_t::f_local(_s[0],_c+4);
    c[2]=k1_t::f_local(_s[0],_c+8);
    c[3]=k1_t::f_local(_s[0],_c+12);

    return k1_t::f_local(_s[1],c);
  }

  /// evaluate first derivative d/dx to f_local() (independent of grid size)
  static value_type fx_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];
    c[0]=k1_t::fx_local(_s[0],_c  );
    c[1]=k1_t::fx_local(_s[0],_c+4);
    c[2]=k1_t::fx_local(_s[0],_c+8);
    c[3]=k1_t::fx_local(_s[0],_c+12);

    return k1_t::f_local(_s[1],c);
  }

  /// evaluate first derivative d/dy to f_local()(independent of grid size)
  static value_type fy_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];
    c[0]=k1_t::f_local(_s[0],_c  );
    c[1]=k1_t::f_local(_s[0],_c+4);
    c[2]=k1_t::f_local(_s[0],_c+8);
    c[3]=k1_t::f_local(_s[0],_c+12);

    return k1_t::fx_local(_s[1],c);
  }

  /// get cell index and local coordinates
  static std::pair<mindex_t,coord_t>
  get_cell(const grid_t& _grid,const typename grid_t::coord_t& _x) {

    auto x =_grid.to_local(_x);
    auto xk=_grid.floor(x);

    coord_t s=x-xk;
    mindex_t k=xk;
    k+=offset();

    for (int i=0;i<DIMS;++i) {
      if (k[i]>=1 && s[i]==real_t(0)) {
        // prevent invalid evaluation at right boundary
        --k[i];
        s[i]=real_t(1);
      }
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
    VC::mvfields::d2::kernel::cubic_convolution_interpolatory_c1.

    Note that the evaluator has a state (in contrast to kernels),
    i.e., multiple instances are required for parallel evaluation!

    \sa d1::evaluator
 */
template <typename KERNEL>
class evaluator {
public:
  typedef KERNEL kernel_t; //!< the kernel
  typedef typename kernel_t::grid_t grid_t; //!< the grid
  typedef typename kernel_t::value_type value_type; //!< value type
  typedef typename kernel_t::real_t real_t; //!< real numbers
  typedef typename kernel_t::coord_t coord_t;


# if !defined(_MSC_VER)  /// tagging cell invalid forces loadc() to call `kernel_t::get_cell()`
  static const index_t INVALIDCELL=std::numeric_limits<index_t>::min()+1;
# else
  enum { INVALIDCELL=INT_MIN+1 };
# endif

  /// create evaluator for `_grid`
  evaluator(const grid_t& _grid) : grid(_grid), cell(INVALIDCELL) {}

  /// load coefficients that are required to evaluate at `_x` to c
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

  /// evaluate function value `_f` at `_x`
  bool value(value_type& _f,const typename grid_t::coord_t& _x) {
    auto p=loadc(_x);
    if (p.first)
      _f=kernel_t::f_local(p.second,c);
    return p.first;
  }
  /// evaluate function gradient `_fx` at `_x`
  bool gradient(value_type* _g,const typename grid_t::coord_t& _x) {
    auto p=loadc(_x);
    if (p.first) {
      auto h=grid.h();
      _g[0]=kernel_t::fx_local(p.second,c)/h[0];
      _g[1]=kernel_t::fy_local(p.second,c)/h[1];
    }
    return p.first;
  }

  value_type    c[kernel_t::SUPPORT_BUFFER_SIZE];  //!< "local" coefficients
  const grid_t& grid; //!< the grid
  index_t       cell; //!< current (linear) cell index
};

# ifndef _MSC_VER
#  undef constexpr // TODO: remove this ASAP!
# endif

//-----------------------------------------------------------------------------
} // namespace d2
} // namespace mvfields
} // namespace VC

#endif // VC_MVFIELDS_GRID_D2_HH
