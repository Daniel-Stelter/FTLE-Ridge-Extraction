//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MVFIELDS_GRID_HH
#define VC_MVFIELDS_GRID_HH

#include <cassert>

#include "grid_d1.hh"
#include "grid_d2.hh"
#include "grid_d3.hh"
#include "grid_d4.hh"

/** \defgroup vc_mvfields_grid Gridded data reconstruction.
    \ingroup vc_mvfields

    This module provides algorithms for gridded data reconstruction,
    which are define in the namespaces

    - VC::mvfields::d1 univariate reconstruction,
    - VC::mvfields::d2 bivariate reconstruction,
    - VC::mvfields::d3 trivariate reconstruction,
    - VC::mvfields::d4 four-variate reconstruction.

    Each namespace (e.g., VC::mvfields::d1) provides classes

    - VC::mvfields::d1::marray to describe a (multi-dimensional) array,
    - VC::mvfields::d1::grid to describe a regular grid over a
      rectangular domain,
    - VC::mvfields::d1::evaluator to evaluate reconstruction kernels
      such as
      VC::mvfields::d1::kernel::cubic_convolution_interpolatory_c1.

    Each namespace provides a namespace for reconstruction kernels, e.g.,
    VC::mvfields::d1::kernel.

    *Notes*

    - The module does *not* provide memory management for
      data. Classes such as VC::mvfields::d2::marray or
      VC::mvfields::d3::grid are just "wrappers" and take a pointer to
      data.

    - Array indexing is _first coordinate fastest_!

    - Data may be any _field_ type, in particular `float`, `double` or
      vectors of `float` or `double` that support the required
      operations.

    - Precision of data may be different from the precision used for
      the numerical reconstruction.

    - Vector values data such as multi-dimensional indices (e.g.,
      VC::mvfields::d2::mindex_t or VC::mvfields::d2::msize_t) or
      coordinates (e.g., VC::mvfields::d3::grid::coord_t) are
      currently based on VC::math::VecN. _This is "hard coded" but
      could be changed easily,_ see, e.g., definition of
      VC::mvfields::d2::vec_t.

    - The reconstruction is provided by **kernel** classes such as
      VC::mvfields::d2::kernel::cubic_convolution_interpolatory_c1 in
      a separate namespaces `kernel`, e.g.,
      VC::mvfields::d2::kernel.<br> The rationale is to first _load_
      the required coefficients (eventually applying boundary
      conditions), e.g., by
      VC::mvfields::d2::kernel::cubic_convolution_interpolatory_c1::get_coeffs(),
      and to then _evaluate_ the function or its derivatives, e.g, by
      VC::mvfields::d2::kernel::cubic_convolution_interpolatory_c1::f_local().<br>
      Note that kernels _don't_ store any persistent information and
      are therefore _stateless_.

    - Reconstruction of **derivatives** by the **kernel**, e.g., by
      VC::mvfields::d2::kernel::cubic_convolution_interpolatory_c1::fx_local(),
      uses **local coordinates** that are independent of the size of
      grid cells, e.g., VC::mvfields::d2::grid::h()!

    - The **evaluator** classes (e.g., VC::mvfields::d1::evaluator)
      simplify the use of kernels. They provide a local storage for
      coefficients (eventually with caching) and are therefore
      _stateful_. (Use kernels directly if this is an issue! The
      implementation of evaluators such as VC::mvfields::d4::evaluator
      shows how to use kernels.)

    - The **evaluator** computes derivatives in "world coordinates",
      i.e., it takes into account the partiuclar parametrization,
      i.e., the size of grid cells (e.g.,
      VC::mvfields::d2::grid::h()).

    - For _documentation_ see VC::mvfields::d1 (most documentation)
      and VC::mvfields::d2 (less documentation)! The structure is
      similar for all dimensions.
 */

# endif // DOXYGEN_SKIP

namespace VC {
namespace mvfields {

//-----------------------------------------------------------------------------

// TODO: kernel: bicubic_linear, tricubic_linear

// TODO: kernel: cubic spline (approximating), provide code for solving the
//                 interpolation problem: data->control points (trivial with
//                 constant bidiagonal factors)

// TODO: kernel: quadratic spline (with data->control points)

// TODO: marray: add copy transpose/reshape? or copy grid->grid

// TODO: try  "__builtin_prefetch (const void *addr, ...)"

// TODO: turn test_grid*.cc into unit tests (but keep framework)

//-----------------------------------------------------------------------------
} // namespace mvfields
} // namespace VC
