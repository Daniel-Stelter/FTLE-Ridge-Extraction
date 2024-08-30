#pragma once
//--------------------------------------------------------------------------//
#include "params.hpp"
//--------------------------------------------------------------------------//
#include "flowsim/datasets/amiraflow.hpp"
#include "flowsim/flow.hpp"
//--------------------------------------------------------------------------//
#include <nlohmann/json.hpp>
//--------------------------------------------------------------------------//
#include <fstream>
#include <memory>
//--------------------------------------------------------------------------//
namespace dst::setup
{
  //--------------------------------------------------------------------------//
  nlohmann::json getJSON(const std::string &json_file_path);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  bool isFieldEmpty(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  bool isFieldString(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  std::string getFieldString(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  bool isFieldStringList(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  std::vector<std::string> getFieldStringList(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  bool isFieldUInt(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  uint getFieldUInt(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  bool isFieldReal(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  real getFieldReal(const nlohmann::json &data, const std::string &field);
  //--------------------------------------------------------------------------//
  template <uint n>
  bool isFieldVecI(const nlohmann::json &data, const std::string &field)
  {
    return data[field].is_array() &&
           data[field].size() == n &&
           std::all_of(data[field].begin(), data[field].end(), [](const nlohmann::json &elem)
                       { return elem.is_number_unsigned(); });
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  VecI<n> getFieldVecI(const nlohmann::json &data, const std::string &field)
  {
    if (!isFieldVecI<n>(data, field))
      throw std::runtime_error("Setup: '" + field + "' should be a list of " + std::to_string(n) + " floats");
    return VecI<n>(data[field].get<std::array<uint, n>>());
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  bool isFieldDomain(const nlohmann::json &data, const std::string &field)
  {
    return data[field].is_array() &&
           data[field].size() == n &&
           std::all_of(data[field].begin(), data[field].end(), [](const nlohmann::json &elem1)
                       { return elem1.is_array() &&
                                elem1.size() == 2 &&
                                std::all_of(elem1.begin(), elem1.end(), [](const nlohmann::json &elem2)
                                            { return elem2.is_number_float(); }); });
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  Domain<n> getFieldDomain(const nlohmann::json &data, const std::string &field)
  {
    if (!isFieldDomain<n>(data, field))
      throw std::runtime_error("Setup: '" + field + "' should be a list of " + std::to_string(n) + " lists of 2 floats");
    Domain<n> domain;
    for (uint i = 0; i < n; ++i)
      domain[i] = Range(data[field][i][0], data[field][i][1]);
    return domain;
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  std::string getDatasetName(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  std::string getSaveDirectory(const nlohmann::json &data);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  enum class Method
  {
    Stelter,
    Kindlmann,
    MarRidges
  };
  //--------------------------------------------------------------------------//
  Method getMethod(const nlohmann::json &data);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  bool hasAmiraFlowDataset(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  bool hasImplementedFlowDataset(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  bool hasImplementedTimeDepScFieldDataset(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  uint getDatasetDimension(const nlohmann::json &data);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  std::shared_ptr<flowsim::Flow<n>> makeImplementedFlow(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  template <uint n>
  std::shared_ptr<flowsim::AmiraFlow<n>> loadAmiraFlow(const nlohmann::json &data, bool load_data);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  bool hasImageResolution(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  VecI<2> getImageResolution(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  bool hasRenderedDomain(const nlohmann::json &data);
  //--------------------------------------------------------------------------//
  Domain<2> getRenderedDomain(const nlohmann::json &data);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  params::ParamsGeneral<n> getParamsGeneral(const nlohmann::json &data)
  {
    params::ParamsGeneral<n> params;
    // required
    params.t = getFieldReal(data, "t");
    params.tau = getFieldReal(data, "tau");
    params.steps = getFieldUInt(data, "steps");
    params.space_res = getFieldVecI<n>(data, "space_res");
    // optional
    params.domain = isFieldEmpty(data, "domain")
                        ? std::optional<Domain<n>>()
                        : getFieldDomain<n>(data, "domain");
    return params;
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  params::ParamsStelter<n> getParamsStelter(const nlohmann::json &data)
  {
    params::ParamsStelter<n> params;
    params.par_init_res = getFieldVecI<n>(data, "par_init_res");
    params.voxel_res = getFieldVecI<n>(data, "voxel_res");
    params.min_ridge_strength = getFieldReal(data, "min_ridge_strength");
    params.elliptic_dis_fac = getFieldReal(data, "elliptic_dis_fac");
    params.full_init_period = getFieldUInt(data, "full_init_period");
    return params;
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  params::ParamsKindlmann<n> getParamsKindlmann(const nlohmann::json &data)
  {
    params::ParamsKindlmann<n> params;
    params.par_init_res = getFieldVecI<n>(data, "par_init_res");
    params.voxel_res = getFieldVecI<n>(data, "voxel_res");
    params.min_ridge_strength = getFieldReal(data, "min_ridge_strength");
    params.sigma = getFieldReal(data, "sigma");
    params.numeric_diff_delta = getFieldReal(data, "numeric_diff_delta");
    return params;
  }
  //--------------------------------------------------------------------------//
  template <uint n>
  params::ParamsMarRidges<n> getParamsMarRidges(const nlohmann::json &data)
  {
    params::ParamsMarRidges<n> params;
    params.min_ridge_strength = getFieldReal(data, "min_ridge_strength");
    return params;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//