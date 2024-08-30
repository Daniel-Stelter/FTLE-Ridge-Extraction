//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MVFIELDS_GRID_D4_HH
#define VC_MVFIELDS_GRID_D4_HH

#include <cassert>
#include <climits> // _MSC_VER

#include <algorithm>
#include <iostream>
#include <utility>

#include "VecN.hh"

#include "grid_d3.hh"

namespace VC {
namespace mvfields {

/// four-variate reconstruction [\ref vc_mvfields_grid] \ingroup vc_mvfields_grid
namespace d4 {
//-----------------------------------------------------------------------------

/// number of dimensions
const int DIMS=4;

/// vector type for indices and coordinates
# if !defined(_MSC_VER)
template <typename T> using vec_t = VC::math::VecN<T,DIMS>;
# else
// TODO: remove this ASAP !!!
#  define constexpr
# endif

typedef int index_t;        //!< linear index type

# if !defined(_MSC_VER)
typedef vec_t<int> mindex_t; //!< multi-dimensional index type
# else
typedef VC::math::VecN<int,DIMS> mindex_t;
# endif

typedef mindex_t msize_t;   //!< multi-dimensional size type

//-----------------------------------------------------------------------------

inline static index_t prod(const mindex_t& _mi) {
  return _mi[0]*_mi[1]*_mi[2]*_mi[3];
}
inline static bool is_zero(const mindex_t& _mi) {
  return _mi[0]==0 && _mi[1]==0 && _mi[2]==0 && _mi[3]==0;
}

//-----------------------------------------------------------------------------

/** \class marray
    \brief Describe 4d array (first index fastest).
    \ingroup vc_mvfields_grid

    \tparam T `value_type`

    This class provides
    \arg indexed access to data,
    \arg copying of blocks.

    Note that marray does *not* manage memory for data!

    \sa d2::marray, d1::marray
 */
template <typename T>
class marray {
public:
  typedef T value_type;

  marray(T* _data,msize_t _size)
    : m_data(_data), m_size(_size) {
  }
  marray(const marray& _other)
    : m_data(_other.m_data), m_size(_other.m_size) {
  }

  static constexpr unsigned ndims() { return DIMS; }
  const msize_t& size() const { return m_size; }
  unsigned numel() const { return prod(m_size); }

  T* data() const { return m_data; }
  T* begin() const { return m_data; }
  T* end() const { return m_data+numel(); }

  index_t index(const mindex_t& _mi) const {
    return ((_mi[3]*m_size[2]+_mi[2])*m_size[1]+_mi[1])*m_size[0]+_mi[0]; // x (0-th dim) fastest
  }
  mindex_t mindex(index_t _i) const {
    int i=_i%m_size[0];
    _i/=m_size[0];
    int j=_i%m_size[1];
    _i/=m_size[1];
    int k=_i%m_size[2];
    _i/=m_size[2];
    return mindex_t(i,j,k,_i);
  }

  template <typename TT>
  void cp_block(TT* _dst,const mindex_t& _mi,const msize_t& _size) const {
    cp_block(_dst,index(_mi),_size);
  }
  template <typename TT>
  void cp_block(TT* _dst,index_t _i,const msize_t& _size) const {
    unsigned n=0;
    for (int p=0;p<_size[3];++p)
      for (int k=0;k<_size[2];++k)
        for (int j=0;j<_size[1];++j)
          for (int i=0;i<_size[0];++i)
            _dst[n++]=TT(m_data[_i+((p*m_size[2]+k)*m_size[1]+j)*m_size[0]+i]);
    // let compiler do optimizations!
  }

  mindex_t overlap(const mindex_t& _mi,const msize_t& _size) const {
    assert(_size[0]<=m_size[0] && _size[1]<=m_size[1] &&
           _size[2]<=m_size[2] && _size[3]<=m_size[3]);

    mindex_t rv(0,0,0,0);

    for (int k=0;k<DIMS;++k) {
      if (_mi[k]<0)                        rv[k]=_mi[k];
      else if (_mi[k]+_size[k]>=m_size[k]) rv[k]=(_mi[k]+_size[k])-m_size[k];
    }
    return rv;
  }

  template <typename TT>
  void cp_block(TT* _dst,const msize_t& _dims,
                const mindex_t& _mi,const msize_t& _size) const {
    cp_block(_dst,_dims,index(_mi),_size);
  }
  template <typename TT>
  void cp_block(TT* _dst,const msize_t& _dims,
                index_t _i,const msize_t& _size) const {
    assert(_dims[1]>=_size[1] && _dims[2]>=_size[2]  && _dims[3]>=_size[3]);
    for (int p=0;p<_size[3];++p)
      for (int k=0;k<_size[2];++k)
        for (int j=0;j<_size[1];++j)
          for (int i=0;i<_size[0];++i)
            _dst[((p*_dims[3]+k)*_dims[2]+j)*_dims[1]+i]=
              m_data[_i+((p*m_size[2]+k)*m_size[1]+j)*m_size[0]+i];
    // let compiler do optimizations!
  }

protected:
  T*      m_data;
  msize_t m_size;
};

//-----------------------------------------------------------------------------

/** \class grid
    \brief Describe 4d regular grid.
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

    \sa d2::grid, d1::grid
 */
template <typename R,typename T>
class grid : public marray<T> {
public:
  typedef R real_t;
# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t;
# else
  typedef VC::math::VecN<R,DIMS> coord_t;
# endif
  typedef T value_type;
  typedef marray<T> array_t;

  /// undefined grid
  grid() : marray<T>(nullptr,mindex_t(0,0,0,0)) {}

  grid(T* _data,msize_t _size)
    : marray<T>(_data,_size) {
    assert(_size[0]>1 && _size[1]>1 && _size[2]>1 && _size[3]>1);
    set_origin(coord_t(0,0,0,0));
    set_extent(coord_t(1,1,1,1));
  }
  grid(const grid& _other)
    : marray<T>(_other),
      m_origin(_other.m_origin), m_extent(_other.m_extent), m_h(_other.m_h) {
  }

  const coord_t& origin() const { return this->m_origin; }
  const coord_t& extent() const { return this->m_extent; }
  const coord_t& h() const { return this->m_h; }

  void set_origin(const coord_t& _origin) { this->m_origin=_origin; }
  void set_extent(const coord_t& _extent) {
    assert(_extent[0]>real_t(0) && _extent[1]>real_t(0) &&
           _extent[2]>real_t(0) && _extent[3]>real_t(0));
    this->m_extent=_extent;
    this->m_h=this->m_extent/
      coord_t(this->size()[0]-1,this->size()[1]-1,
              this->size()[2]-1,this->size()[3]-1);
  }
  template <typename RR>
  void set_bbox(const RR* _bbox) {
    R bbox[2*DIMS];
    for (int k=0;k<DIMS;++k) {
      bbox[k+k  ]=std::min(_bbox[k+k  ],_bbox[k+k+1]);
      bbox[k+k+1]=std::max(_bbox[k+k  ],_bbox[k+k+1]);
    }
    set_origin(coord_t(bbox[0],bbox[2],bbox[4],bbox[6]));
    set_extent(coord_t(bbox[1],bbox[3],bbox[5],bbox[7])-this->m_origin);
  }

  bool inside(const coord_t& _x) const {
    return
      this->m_origin[0]<=_x[0] && _x[0]<=this->m_origin[0]+this->m_extent[0] &&
      this->m_origin[1]<=_x[1] && _x[1]<=this->m_origin[1]+this->m_extent[1] &&
      this->m_origin[2]<=_x[2] && _x[2]<=this->m_origin[2]+this->m_extent[2] &&
      this->m_origin[3]<=_x[3] && _x[3]<=this->m_origin[3]+this->m_extent[3];
  }

  coord_t to_world(const coord_t& _x) const {
    return this->m_origin+_x*m_h;
  }

  coord_t to_local(const coord_t& _x) const {
    return (_x-this->m_origin)/m_h;
  }

  static coord_t floor(const coord_t& _x) {
    return coord_t(::floor(_x[0]),::floor(_x[1]),
                   ::floor(_x[2]),::floor(_x[3]));
  }

protected:
  coord_t m_origin;
  coord_t m_extent;
  coord_t m_h;
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
# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t; //!< coordinates
# else
  typedef VC::math::VecN<R,DIMS> coord_t;
# endif
  typedef d1::kernel::linear_interpolatory_c0<grid_t,value_type,real_t> k1_t;
  typedef d2::kernel::linear_interpolatory_c0<grid_t,value_type,real_t> k2_t;
  typedef d3::kernel::linear_interpolatory_c0<grid_t,value_type,real_t> k3_t;

  enum { SUPPORT_BUFFER_SIZE=16 }; // _MSC_VER fix: should be support_buffer_size()

  /// size of support
  static msize_t support_size() { return msize_t(2,2,2,2); }
  /// dimensions of of support (same as support_size() for 1d)
  static constexpr size_t support_buffer_size() { return 16; }
  /// offset to index returned by get_cell()
  static constexpr size_t offset() { return 0; }

  /** Copy support_size() coefficients from `grid` at `_idx`.
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
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::f_local(s3,_c  );
    c[1]=k3_t::f_local(s3,_c+8);
    return k1_t::f_local(_s[3],c);
  }

  /// evaluate first derivative d/dx to f_local() (independent of grid size)
  static value_type fx_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::fx_local(s3,_c  );
    c[1]=k3_t::fx_local(s3,_c+8);
    return k1_t::f_local(_s[3],c);
  }

  /// evaluate first derivative d/dy to f_local() (independent of grid size)
  static value_type fy_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::fy_local(s3,_c  );
    c[1]=k3_t::fy_local(s3,_c+8);
    return k1_t::f_local(_s[3],c);
  }

  /// evaluate first derivative d/dz to f_local() (independent of grid size)
  static value_type fz_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::fz_local(s3,_c  );
    c[1]=k3_t::fz_local(s3,_c+8);
    return k1_t::f_local(_s[3],c);
  }

  /// evaluate first derivative d/dt to f_local() (independent of grid size)
  static value_type ft_local(const coord_t& _s,const value_type* _c) {
    value_type c[2];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::f_local(s3,_c  );
    c[1]=k3_t::f_local(s3,_c+8);
    return k1_t::fx_local(_s[3],c);
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

    \sa d1::cubic_convolution_interpolatory_c1
 */
template <typename GRID,
          typename T=typename GRID::value_type,
          typename R=typename GRID::real_t>
class cubic_convolution_interpolatory_c1 {
public:
  typedef GRID grid_t;
  typedef T value_type;
  typedef R real_t;

# if !defined(_MSC_VER)
  typedef vec_t<R> coord_t;
# else
  typedef VC::math::VecN<R,DIMS> coord_t;
# endif
  typedef d1::kernel::cubic_convolution_interpolatory_c1<grid_t,value_type,real_t> k1_t;
  typedef d2::kernel::cubic_convolution_interpolatory_c1<grid_t,value_type,real_t> k2_t;
  typedef d3::kernel::cubic_convolution_interpolatory_c1<grid_t,value_type,real_t> k3_t;

  enum { SUPPORT_BUFFER_SIZE=256 }; // _MSC_VER fix: should be support_buffer_size()

  static msize_t support_size() { return msize_t(4,4,4,4); }
  static constexpr size_t support_buffer_size() { return 256; }
  static constexpr size_t offset() { return -1; }

  static void set_left_bcond(value_type* _c,const msize_t& _inc) {
    for (int i=0;i<4;++i)
      k3_t::set_left_bcond(_c+i*_inc[0],d3::msize_t(_inc[1],_inc[2],_inc[3]));
  }
  static void set_right_bcond(value_type* _c,const msize_t& _inc) {
    for (int i=0;i<4;++i)
      k3_t::set_right_bcond(_c+i*_inc[0],d3::msize_t(_inc[1],_inc[2],_inc[3]));
  }

  static bool get_coeffs(const grid_t& _grid,
                         const mindex_t& _mi,value_type* _c) {

    mindex_t overlap=_grid.overlap(_mi,support_size());

    if (is_zero(overlap))
      _grid.cp_block(_c,_mi,support_size());
    else {
      // copy partial block
      {
        value_type* c0=_c;
        mindex_t    mi(_mi);     // new origin
        mindex_t    sz(4,4,4,4); // new size
        const int   INC[DIMS]={1,4,16,64};

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

        _grid.cp_block(c0,msize_t(1,4,4,4),mi,sz);
      }
      // compute missing coefficients
      if (overlap[0]<0)
        set_left_bcond(_c,msize_t(64,16,4,1));
      else if (overlap[0]>0)
        set_right_bcond(_c,msize_t(64,16,4,1));

      if (overlap[1]<0)
        set_left_bcond(_c,msize_t(64,16,1,4));
      else if (overlap[1]>0)
        set_right_bcond(_c,msize_t(64,16,1,4));

      if (overlap[2]<0)
        set_left_bcond(_c,msize_t(1,4,64,16));
      else if (overlap[2]>0)
        set_right_bcond(_c,msize_t(1,4,64,16));

      if (overlap[3]<0)
        set_left_bcond(_c,msize_t(1,4,16,64));
      else if (overlap[3]>0)
        set_right_bcond(_c,msize_t(1,4,16,64));
    }
    return true;
  }

  static value_type f_local(const coord_t& _s,const value_type* _c) {
# if 1
    // faster!
    value_type c[4];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::f_local(s3,_c    );
    c[1]=k3_t::f_local(s3,_c+ 64);
    c[2]=k3_t::f_local(s3,_c+128);
    c[3]=k3_t::f_local(s3,_c+192);

    return k1_t::f_local(_s[3],c);
# else
    // slower!
    value_type c[16];
    typename k2_t::coord_t s2(_s[0],_s[1]);
    for (int i=0;i<16;++i)
      c[i]=k2_t::f_local(s2,_c+i*16);
    return k2_t::f_local(typename k2_t::coord_t(_s[2],_s[3]),c);
# endif
  }

  static value_type fx_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::fx_local(s3,_c    );
    c[1]=k3_t::fx_local(s3,_c+ 64);
    c[2]=k3_t::fx_local(s3,_c+128);
    c[3]=k3_t::fx_local(s3,_c+192);

    return k1_t::f_local(_s[3],c);
  }

  static value_type fy_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::fy_local(s3,_c    );
    c[1]=k3_t::fy_local(s3,_c+ 64);
    c[2]=k3_t::fy_local(s3,_c+128);
    c[3]=k3_t::fy_local(s3,_c+192);

    return k1_t::f_local(_s[3],c);
  }

  static value_type fz_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];

    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::fz_local(s3,_c    );
    c[1]=k3_t::fz_local(s3,_c+ 64);
    c[2]=k3_t::fz_local(s3,_c+128);
    c[3]=k3_t::fz_local(s3,_c+192);

    return k1_t::f_local(_s[3],c);
  }

  /// evaluate first derivative d/dt to f_local()
  static value_type ft_local(const coord_t& _s,const value_type* _c) {
    value_type c[4];
    typename k3_t::coord_t s3(_s[0],_s[1],_s[2]);
    c[0]=k3_t::f_local(s3,_c    );
    c[1]=k3_t::f_local(s3,_c+ 64);
    c[2]=k3_t::f_local(s3,_c+128);
    c[3]=k3_t::f_local(s3,_c+192);

    return k1_t::fx_local(_s[3],c);
  }


  static std::pair<mindex_t,coord_t>
  get_cell(const grid_t& _grid,const typename grid_t::coord_t& _x) {

    auto x =_grid.to_local(_x);
    auto xk=_grid.floor(x);

    coord_t  s=x-xk;
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
    VC::mvfields::d4::cubic_convolution_interpolatory_c1.

    Note that the evaluator has a state (in contrast to kernels),
    i.e., multiple instances are required for parallel evaluation!

    \sa d2::evaluator
 */
template <typename KERNEL>
class evaluator {
public:
  typedef KERNEL kernel_t;
  typedef typename kernel_t::grid_t grid_t;
  typedef typename kernel_t::value_type value_type;
  typedef typename kernel_t::real_t real_t;
  typedef typename kernel_t::coord_t coord_t;

# if !defined(_MSC_VER)
  static const index_t INVALIDCELL=std::numeric_limits<index_t>::min()+1;
# else
  enum { INVALIDCELL=INT_MIN+1 };
# endif

  evaluator(const grid_t& _grid) : grid(_grid), cell(-1) {}

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

  bool value(value_type& _f,const typename grid_t::coord_t& _x) {
    auto p=loadc(_x);
    if (p.first)
      _f=kernel_t::f_local(p.second,c);
    return p.first;
  }
  bool gradient(value_type* _g,const typename grid_t::coord_t& _x) {
    auto p=loadc(_x);
    if (p.first) {
      auto h=grid.h();
      _g[0]=kernel_t::fx_local(p.second,c)/h[0];
      _g[1]=kernel_t::fy_local(p.second,c)/h[1];
      _g[2]=kernel_t::fz_local(p.second,c)/h[2];
      _g[3]=kernel_t::ft_local(p.second,c)/h[3];
    }
    return p.first;
  }

  value_type    c[kernel_t::SUPPORT_BUFFER_SIZE];
  const grid_t& grid;
  index_t       cell;
};

# ifndef _MSC_VER
#  undef constexpr // TODO: remove this ASAP!
# endif

//-----------------------------------------------------------------------------
} // namespace d4
} // namespace mvfields
} // namespace VC

#endif // VC_MVFIELDS_GRID_D4_HH
