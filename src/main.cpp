#include "experiments.hpp"
#include "setup.hpp"
//--------------------------------------------------------------------------//
using namespace dst;
//--------------------------------------------------------------------------//
template <uint n>
using ItplScFieldSeries_ptr_t = std::shared_ptr<flowsim::InterpolScalarFieldSeries<n>>;
//--------------------------------------------------------------------------//
template <uint n>
using ScFieldSeriesBase_t = dst::flowsim::ScalarFieldSeriesBase<n>;
// -------------------------------------------------------------------------//
template <uint n>
using ScFieldSeriesBase_ptr_t = std::shared_ptr<ScFieldSeriesBase_t<n>>;
// -------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
template <uint n>
ItplScFieldSeries_ptr_t<n> prepareItplScFieldSeries(const nlohmann::json &data)
{
  params::ParamsGeneral<n> params = setup::getParamsGeneral<n>(data);
  std::string base_dir = setup::getSaveDirectory(data);

  ItplScFieldSeries_ptr_t<n> p_itpl_sc_series;
  if (setup::hasImplementedFlowDataset(data))
  {
    auto p_flow = setup::makeImplementedFlow<n>(data);
    Domain<n> domain = params.domain.value_or(p_flow->spaceDomainStandard());
    SteppedRange stepped_tau(params.tau, params.steps, false);
    p_itpl_sc_series = experiments::loadOrCreateItplFTLEFieldSeries<n>(p_flow,
                                                                       params.t,
                                                                       stepped_tau,
                                                                       params.space_res,
                                                                       domain,
                                                                       base_dir);
  }
  else if (setup::hasAmiraFlowDataset(data))
  {
    auto p_flow = setup::loadAmiraFlow<n>(data, true);
    Domain<n> domain = params.domain.value_or(p_flow->spaceDomainStandard());
    SteppedRange stepped_tau(params.tau, params.steps, false);
    p_itpl_sc_series = experiments::loadOrCreateItplFTLEFieldSeries<n>(p_flow,
                                                                       params.t,
                                                                       stepped_tau,
                                                                       params.space_res,
                                                                       domain,
                                                                       base_dir);
    p_flow->dropLoadedData(); // prevent unnecessary memory consumption of flow data
  }
  else
    throw std::runtime_error("Setup: no valid dataset provided");
  return p_itpl_sc_series;
}
//--------------------------------------------------------------------------//
void generateImages2D(const nlohmann::json &data,
                      ScFieldSeriesBase_ptr_t<2> &p_sc_series_base)
{
  if (setup::hasImageResolution(data))
  {
    // obtain image data
    VecI<2> img_res = setup::getImageResolution(data);
    Domain<2> rendered_domain = setup::hasRenderedDomain(data)
                                    ? setup::getRenderedDomain(data)
                                    : p_sc_series_base->domain();
    // create images
    std::string base_dir = setup::getSaveDirectory(data);
    experiments::createScFieldSeriesImages(p_sc_series_base, img_res, rendered_domain, base_dir);
    experiments::createScFieldSeriesRidgeImagesWithValue(p_sc_series_base, img_res, rendered_domain, base_dir);
    // experiments::createScFieldSeriesRidgeImagesOnTexture(p_sc_series_base, base_dir);
  }
  else
    std::cout << "WARNING: image resolution missing - skip creating images" << std::endl;
}
//--------------------------------------------------------------------------//
void generateParSysOutput3D(const nlohmann::json &data,
                            ScFieldSeriesBase_ptr_t<3> &p_sc_series_base)
{
  std::string base_dir = setup::getSaveDirectory(data);
  setup::Method method = setup::getMethod(data);
  auto par_sys_normals = setup::Method::Stelter == method
                             ? experiments::ParSysNormals::Stelter
                             : experiments::ParSysNormals::Kindlmann;
  experiments::createParSysBlenderInput(p_sc_series_base, par_sys_normals, base_dir);
  // experiments::createParSysCloudCompareInput(p_sc_series_base, par_sys_normals, base_dir);
}

//--------------------------------------------------------------------------//
void executeMarRidges3D(const nlohmann::json &data)
{
  // obtain and print setup
  auto params_general = setup::getParamsGeneral<3>(data);
  auto params_mar_ridges = setup::getParamsMarRidges<3>(data);
  std::string base_dir = setup::getSaveDirectory(data);
  std::cout << params_general << '\n'
            << params_mar_ridges << std::endl;
  utils::printSeparator();
  // start execution
  auto p_itpl_sc_field = prepareItplScFieldSeries<3>(data);
  experiments::executeMarRidgesSeries(p_itpl_sc_field, params_mar_ridges, base_dir);
}
//--------------------------------------------------------------------------//
void executeStelter2D(const nlohmann::json &data)
{
  // obtain and print setup
  auto params_general = setup::getParamsGeneral<2>(data);
  auto params_stelter = setup::getParamsStelter<2>(data);
  std::string base_dir = setup::getSaveDirectory(data);
  std::cout << params_general << '\n'
            << params_stelter << std::endl;
  utils::printSeparator();
  // start execution
  ItplScFieldSeries_ptr_t<2> p_itpl_sc_series = prepareItplScFieldSeries<2>(data);
  experiments::executeParSysSeriesStelter(p_itpl_sc_series, params_stelter, base_dir);
  // cast pointer and generate images
  ScFieldSeriesBase_ptr_t<2> p_sc_series_base = std::static_pointer_cast<ScFieldSeriesBase_t<2>>(p_itpl_sc_series);
  generateImages2D(data, p_sc_series_base);
}
//--------------------------------------------------------------------------//
void executeKindlmann2D(const nlohmann::json &data)
{
  auto params_general = setup::getParamsGeneral<2>(data);
  auto params_kindlmann = setup::getParamsKindlmann<2>(data);
  std::string base_dir = setup::getSaveDirectory(data);
  std::cout << params_general << '\n'
            << params_kindlmann << std::endl;
  utils::printSeparator();
  // start execution
  ItplScFieldSeries_ptr_t<2> p_itpl_sc_series = prepareItplScFieldSeries<2>(data);
  experiments::executeParSysSeriesKindlmann(p_itpl_sc_series, params_kindlmann, base_dir);
  // cast pointer and generate images
  ScFieldSeriesBase_ptr_t<2> p_sc_series_base = std::static_pointer_cast<ScFieldSeriesBase_t<2>>(p_itpl_sc_series);
  generateImages2D(data, p_sc_series_base);
}
//--------------------------------------------------------------------------//
void executeStelter3D(const nlohmann::json &data)
{
  // obtain and print setup
  auto params_general = setup::getParamsGeneral<3>(data);
  auto params_stelter = setup::getParamsStelter<3>(data);
  std::string base_dir = setup::getSaveDirectory(data);
  std::cout << params_general << '\n'
            << params_stelter << std::endl;
  utils::printSeparator();
  // start execution
  ItplScFieldSeries_ptr_t<3> p_itpl_sc_series = prepareItplScFieldSeries<3>(data);
  experiments::executeParSysSeriesStelter(p_itpl_sc_series, params_stelter, base_dir);
  // cast pointer and generate output for external rendering software
  ScFieldSeriesBase_ptr_t<3> p_sc_series_base = std::static_pointer_cast<ScFieldSeriesBase_t<3>>(p_itpl_sc_series);
  generateParSysOutput3D(data, p_sc_series_base);
}
//--------------------------------------------------------------------------//
void executeKindlmann3D(const nlohmann::json &data)
{
  auto params_general = setup::getParamsGeneral<3>(data);
  auto params_kindlmann = setup::getParamsKindlmann<3>(data);
  std::string base_dir = setup::getSaveDirectory(data);
  std::cout << params_general << '\n'
            << params_kindlmann << std::endl;
  utils::printSeparator();
  // start execution
  ItplScFieldSeries_ptr_t<3> p_itpl_sc_series = prepareItplScFieldSeries<3>(data);
  experiments::executeParSysSeriesKindlmann(p_itpl_sc_series, params_kindlmann, base_dir);
  // cast pointer and generate output for external rendering software
  ScFieldSeriesBase_ptr_t<3> p_sc_series_base = std::static_pointer_cast<ScFieldSeriesBase_t<3>>(p_itpl_sc_series);
  generateParSysOutput3D(data, p_sc_series_base);
}
//--------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
int main(int argc, char **argv)
{
  utils::Timer::logCurrentTime();
  // check if provided file is valid
  if (argc < 2)
    throw std::runtime_error("No setup file provided");
  std::string json_file_path = argv[1];

  auto data = setup::getJSON(json_file_path);
  std::string base_dir = setup::getSaveDirectory(data);
  std::cout << "Dataset: " << setup::getDatasetName(data)
            << "\nSave directory: " << base_dir << std::endl;

  setup::Method method = setup::getMethod(data);
  uint dim = setup::getDatasetDimension(data);

  // handle different approaches
  if (setup::Method::MarRidges == method) // Marching Ridges
  {
    // checking and printing
    if (3 != dim)
      throw std::runtime_error("Invalid dataset: Marching Ridges is only implemented for 3D datasets");
    executeMarRidges3D(data);
  }
  else // particle system (Stelter or Kindlmann)
  {
    if (2 == dim)
    {
      if (setup::Method::Stelter == method)
        executeStelter2D(data);
      else
        executeKindlmann2D(data);
    }
    else if (3 == dim)
    {
      if (setup::Method::Stelter == method)
        executeStelter3D(data);
      else
        executeKindlmann3D(data);
    }
  }
}
//--------------------------------------------------------------------------//
