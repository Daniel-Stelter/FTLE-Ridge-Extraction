//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_DEBUG_HH
#define VC_BASE_DEBUG_HH

// NOTE: Removed all dependencies, which effectively disables
//       all debuggungs macros (same effect as defining NDEBUG)

#  define VC_DBG_LOCATION()
#  define VC_DBG_HERE(msg)
#  define VC_DBG_HERE_SHORT(msg)
#  define VC_DBG_TRACE(msg)
#  define VC_DBG_P(var)
#  define VC_DBG_PS(var)
#  define VC_DBG_PQS(var)
#  define VC_DBG_PTYPE(var)
#  define VC_DBG_PTYPEOF(var)
#  define VC_DBG_P_IF(condition,var)
#  define VC_DBG_COUNT(msg)
#  define VC_DBG_TRACE_ONCE(msg)
#  define VC_DBG_FIXME(msg)
#  define VC_DBG_BACKTRACE(msg)
#  define VC_DBG_BREAK(msg)
#  define VC_DBG_BREAK_IF(cond)
#  define VC_DBG_DBREAK(msg)
#  define VC_DBG_DBREAK_IF(cond)
#  define VC_DBG_SCOPE(msg)
#  define VC_DBG_TIC()
#  define VC_DBG_TOC()

//=============================================================================
#endif // VC_BASE_DEBUG_HH defined
