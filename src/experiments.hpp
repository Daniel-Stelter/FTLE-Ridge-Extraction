#pragma once
//--------------------------------------------------------------------------//
#include "flowsim/ftlefield.hpp"
#include "flowsim/marchingridges.hpp"
#include "flowsim/particlesystem.hpp"
#include "images/visualization.hpp"
//--------------------------------------------------------------------------//
namespace dst::experiments
{
  //--------------------------------------------------------------------------//
  // save directories (on creation only):
  struct SubDirectories
  {
    SubDirectories() = delete;

    // clang-format off
    constexpr static char const* DATA_SC_FIELD        = "/data-saves/sc-field-grid/";   // file ending: .data
    constexpr static char const* DATA_PARTICLES       = "/data-saves/ridge-particles/"; // file ending: .data

    constexpr static char const* IMGS_SC_FIELD        = "/out/sc-field/";               // file ending: .ppm
    constexpr static char const* IMGS_RIDGES_TEXTURE  = "/out/ridges-texture/";         // file ending: .ppm
    constexpr static char const* IMGS_RIDGES_VALS     = "/out/ridges-vals/";            // file ending: .ppm

    constexpr static char const* OUT_CLOUDCOMPARE     = "/out/cloudcompare/";           // file ending: .bin
    constexpr static char const* OUT_BLENDER          = "/out/blender/";                // file ending: .data
    // clang-format on
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  std::shared_ptr<flowsim::InterpolFTLEFieldSeries<n>>
  loadOrCreateItplFTLEFieldSeries(std::shared_ptr<flowsim::Flow<n>> p_flow,
                                  real t,
                                  const SteppedRange &stepped_tau,
                                  const VecI<n> &resolution,
                                  const Domain<n> &domain,
                                  const std::string &base_dir)
  {
    using namespace std;
    using namespace flowsim;
    // subdirectory for saving
    string save_dir = base_dir + SubDirectories::DATA_SC_FIELD;
    // if save exists: try to load
    if (InterpolFTLEFieldSeries<n>::existsSave(save_dir))
    {
      // load the data and check if it coincides to current user input
      shared_ptr<InterpolFTLEFieldSeries<n>> p_itpl_ftle_series =
          make_shared<InterpolFTLEFieldSeries<n>>(p_flow, t, save_dir);
      if (!p_itpl_ftle_series->hasSameParams(stepped_tau, resolution, domain))
        throw runtime_error("Loading InterpolFTLEFieldSeries: different saved and provided parameters!");
      else
        return p_itpl_ftle_series;
    }
    // if no save exists: initialize series
    else
      return make_shared<InterpolFTLEFieldSeries<n>>(p_flow, t, stepped_tau, resolution, domain, save_dir);
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  void
  executeParSysSeriesStelter(std::shared_ptr<flowsim::InterpolScalarFieldSeries<n>> &p_itpl_sc_series,
                             const params::ParamsStelter<n> &params,
                             const std::string &base_dir)
  {
    std::string save_dir = base_dir + SubDirectories::DATA_PARTICLES;
    flowsim::ParticleSystemSeriesStelter<n> par_sys_series(p_itpl_sc_series, params);
    par_sys_series.searchAndSaveRidges(save_dir);
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  void
  executeParSysSeriesKindlmann(std::shared_ptr<flowsim::InterpolScalarFieldSeries<n>> &p_itpl_sc_series,
                               const params::ParamsKindlmann<n> &params,
                               const std::string &base_dir)
  {
    std::string save_dir = base_dir + SubDirectories::DATA_PARTICLES;
    flowsim::ParticleSystemSeriesKindlmann<n> par_sys_series(p_itpl_sc_series, params);
    par_sys_series.searchAndSaveRidges(save_dir);
  }
  //--------------------------------------------------------------------------//
  void
  executeMarRidgesSeries(std::shared_ptr<flowsim::InterpolScalarFieldSeries<3>> &p_itpl_sc_series,
                         const params::ParamsMarRidges<3> &params,
                         const std::string &base_dir);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  createScFieldSeriesImages(std::shared_ptr<flowsim::ScalarFieldSeriesBase<2>> &sc_field_series,
                            const VecI<2> &img_resolution,
                            const std::string &base_dir);
  //--------------------------------------------------------------------------//
  void
  createScFieldSeriesImages(std::shared_ptr<flowsim::ScalarFieldSeriesBase<2>> &sc_field_series,
                            const VecI<2> &img_resolution,
                            const Domain<2> &rendered_domain,
                            const std::string &base_dir);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  createScFieldSeriesRidgeImagesOnTexture(std::shared_ptr<flowsim::ScalarFieldSeriesBase<2>> &p_sc_field_series,
                                          const std::string &base_dir);
  //--------------------------------------------------------------------------//
  void
  createScFieldSeriesRidgeImagesOnTexture(std::shared_ptr<flowsim::ScalarFieldSeriesBase<2>> &p_sc_field_series,
                                          const Domain<2> &rendered_domain,
                                          const std::string &base_dir);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  createScFieldSeriesRidgeImagesWithValue(std::shared_ptr<flowsim::ScalarFieldSeriesBase<2>> &p_sc_field_series,
                                          const VecI<2> &img_resolution,
                                          const std::string &base_dir);
  //--------------------------------------------------------------------------//
  void
  createScFieldSeriesRidgeImagesWithValue(std::shared_ptr<flowsim::ScalarFieldSeriesBase<2>> &p_sc_field_series,
                                          const VecI<2> &img_resolution,
                                          const Domain<2> &rendered_domain,
                                          const std::string &base_dir);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  enum class ParSysNormals
  {
    Stelter,
    Kindlmann
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  createParSysCloudCompareInput(std::shared_ptr<flowsim::ScalarFieldSeriesBase<3>> &p_sc_field_series,
                                ParSysNormals par_sys_normals,
                                const std::string &base_dir);
  //--------------------------------------------------------------------------//
  void
  createParSysCloudCompareInput(std::shared_ptr<flowsim::ScalarFieldSeriesBase<3>> &p_sc_field_series,
                                std::shared_ptr<flowsim::ScalarFieldSeriesBase<3>> &p_sc_field_series_for_values,
                                ParSysNormals par_sys_normals,
                                const std::string &base_dir);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  createParSysBlenderInput(std::shared_ptr<flowsim::ScalarFieldSeriesBase<3>> &p_sc_field_series,
                           ParSysNormals par_sys_normals,
                           const std::string &base_dir);
  //--------------------------------------------------------------------------//
  void
  createParSysBlenderInput(std::shared_ptr<flowsim::ScalarFieldSeriesBase<3>> &p_sc_field_series,
                           std::shared_ptr<flowsim::ScalarFieldSeriesBase<3>> &p_sc_field_series_for_values,
                           ParSysNormals par_sys_normals,
                           const std::string &base_dir);
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//