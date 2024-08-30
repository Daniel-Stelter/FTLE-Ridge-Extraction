#include "setup.hpp"
//--------------------------------------------------------------------------//
#include "flowsim/datasets/abc.hpp"
#include "flowsim/datasets/doublegyre.hpp"
#include "flowsim/datasets/forcedduffing.hpp"
#include "flowsim/datasets/tornado.hpp"
//--------------------------------------------------------------------------//
namespace dst::setup
{
  //--------------------------------------------------------------------------//
  nlohmann::json getJSON(const std::string &json_file_path)
  {
    std::ifstream in(json_file_path);
    if (!in)
      throw std::runtime_error("Setup: could not read from file '" + json_file_path + "'");
    nlohmann::json data;
    try
    {
      data = nlohmann::json::parse(in);
    }
    catch (const std::exception &)
    {
      throw std::runtime_error("Setup: file '" + json_file_path + "' seems to have invalid JSON format");
    }
    return data;
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  bool isFieldEmpty(const nlohmann::json &data, const std::string &field)
  {
    return !data.contains(field) || data[field].is_null();
  }
  //--------------------------------------------------------------------------//
  bool isFieldString(const nlohmann::json &data, const std::string &field)
  {
    return data[field].is_string();
  }
  //--------------------------------------------------------------------------//
  std::string getFieldString(const nlohmann::json &data, const std::string &field)
  {
    if (!isFieldString(data, field))
      throw std::runtime_error("Setup: '" + field + "' should be a string");
    return data[field].get<std::string>();
  }
  //--------------------------------------------------------------------------//
  bool isFieldStringList(const nlohmann::json &data, const std::string &field)
  {
    return data[field].is_array() &&
           data[field].size() > 0 &&
           std::all_of(data[field].begin(), data[field].end(), [](const nlohmann::json &elem)
                       { return elem.is_string(); });
  }
  //--------------------------------------------------------------------------//
  std::vector<std::string> getStringList(const nlohmann::json &data, const std::string &field)
  {
    if (!isFieldStringList(data, field))
      throw std::runtime_error("Setup: '" + field + "' should be a list of at least one string");
    std::vector<std::string> list;
    list.reserve(data[field].size());
    for (const nlohmann::json &file : data[field])
      list.push_back(file.get<std::string>());
    return list;
  }
  //--------------------------------------------------------------------------//
  bool isFieldUInt(const nlohmann::json &data, const std::string &field)
  {
    return data[field].is_number_unsigned();
  }
  //--------------------------------------------------------------------------//
  uint getFieldUInt(const nlohmann::json &data, const std::string &field)
  {
    if (!isFieldUInt(data, field))
      throw std::runtime_error("Setup: '" + field + "' should be an unsigned int");
    return data[field].get<uint>();
  }
  //--------------------------------------------------------------------------//
  bool isFieldReal(const nlohmann::json &data, const std::string &field)
  {
    return data[field].is_number_float();
  }
  //--------------------------------------------------------------------------//
  real getFieldReal(const nlohmann::json &data, const std::string &field)
  {
    if (!isFieldReal(data, field))
      throw std::runtime_error("Setup: '" + field + "' should be a float");
    return data[field].get<real>();
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  std::string getDatasetName(const nlohmann::json &data)
  {
    return getFieldString(data, "dataset");
  }
  //--------------------------------------------------------------------------//
  std::string getSaveDirectory(const nlohmann::json &data)
  {
    std::string save_dir = getFieldString(data, "save_dir");
    if ("~/" == save_dir.substr(0, 2))
      return getenv("HOME") + save_dir.substr(1);
    return save_dir;
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  Method getMethod(const nlohmann::json &data)
  {
    std::string method = getFieldString(data, "method");
    if (std::string("Stelter") == method)
      return Method::Stelter;
    if (std::string("Kindlmann") == method)
      return Method::Kindlmann;
    if (std::string("Marching Ridges") == method)
      return Method::MarRidges;
    throw std::runtime_error("Setup: 'method' should be 'Stelter', 'Kindlmann' or 'Marching Ridges'");
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  bool hasAmiraFlowDataset(const nlohmann::json &data)
  {
    std::string name = getDatasetName(data);
    return std::string("Amira-") == name.substr(0, 6);
  }
  //--------------------------------------------------------------------------//
  bool hasImplementedFlowDataset(const nlohmann::json &data)
  {
    std::string name = getDatasetName(data);
    for (std::string dataset : {"DoubleGyre2D", "ForcedDuffing2D",
                                "ABC3D", "DoubleGyre3D", "Tornado3D"})
    {
      if (dataset == name)
        return true;
    }
    return false;
  }
  //--------------------------------------------------------------------------//
  uint getDatasetDimension(const nlohmann::json &data)
  {
    std::string name = getDatasetName(data);
    if (std::string("2D") == name.substr(name.length() - 2, name.length()))
      return 2;
    if (std::string("3D") == name.substr(name.length() - 2, name.length()))
      return 3;
    throw std::runtime_error("Setup: 'dataset' must end with either 2D or 3D");
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <>
  std::shared_ptr<flowsim::Flow<2>> makeImplementedFlow<2>(const nlohmann::json &data)
  {
    using namespace std;
    using namespace flowsim;

    std::string name = getDatasetName(data);
    // try to find correct flow
    if ("DoubleGyre2D" == name)
      return make_shared<DoubleGyre2D>();
    if ("ForcedDuffing2D" == name)
      return make_shared<ForcedDuffing2D>();
    // if no flow was found, throw exception
    throw runtime_error("Setup: makeImplementedFlow<2> failed (2D dataset '" + name + "' is not implemented)");
  }
  //--------------------------------------------------------------------------//
  template <>
  std::shared_ptr<flowsim::Flow<3>> makeImplementedFlow<3>(const nlohmann::json &data)
  {
    using namespace std;
    using namespace flowsim;

    std::string name = getDatasetName(data);
    // try to find correct flow
    if ("ABC3D" == name)
      return make_shared<ABC3D>();
    if ("DoubleGyre3D" == name)
      return make_shared<DoubleGyre3D>();
    if ("Tornado3D" == name)
      return make_shared<Tornado3D>();
    // if no flow was found, throw exception
    throw runtime_error("Setup: makeImplementedFlow<3> failed (3D dataset '" + name + "' is not implemented)");
  }
  //--------------------------------------------------------------------------//
  template <>
  std::shared_ptr<flowsim::AmiraFlow<2>> loadAmiraFlow<2>(const nlohmann::json &data, bool load_data) // const Range &time_range
  {
    using namespace std;
    using namespace flowsim;

    std::string name = getDatasetName(data);
    // check dimension
    if (2 == getDatasetDimension(data))
      throw runtime_error("Setup: loadAmiraFlow<2> failed (dataset dimension was not 2)");
    // check Amira flow
    if (!hasAmiraFlowDataset(data))
      throw runtime_error("Setup: loadAmiraFlow<2> failed (dataset must begin with 'Amira-')");
    // check Amira files
    if (!isFieldStringList(data, "amira_files"))
      throw runtime_error("Setup: loadAmiraFlow<2> failed (no input files for Amira dataset provided)");
    vector<string> amira_files = getStringList(data, "amira_files");
    // create and return flow
    return make_shared<AmiraFlow2D>(amira_files[0], name, load_data);
  }
  //--------------------------------------------------------------------------//
  template <>
  std::shared_ptr<flowsim::AmiraFlow<3>> loadAmiraFlow<3>(const nlohmann::json &data, bool load_data)
  {
    using namespace std;
    using namespace flowsim;

    real t = getFieldReal(data, "t");
    real tau = getFieldReal(data, "tau");
    Range time_range(t, t + tau);

    std::string name = getDatasetName(data);
    // check dimension
    if (3 == getDatasetDimension(data))
      throw runtime_error("Setup: loadAmiraFlow<3> failed (dataset dimension was not 3)");
    // check Amira flow
    if (!hasAmiraFlowDataset(data))
      throw runtime_error("Setup: loadAmiraFlow<3> failed (dataset must begin with 'Amira-')");
    // check Amira files
    if (!isFieldStringList(data, "amira_files"))
      throw runtime_error("Setup: loadAmiraFlow<3> failed (no input files for Amira dataset provided)");
    vector<string> amira_files = getStringList(data, "amira_files");
    // create and return flow
    return make_shared<AmiraFlow3D>(amira_files, time_range, name, load_data);
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  bool hasImageResolution(const nlohmann::json &data) { return isFieldVecI<2>(data, "img_res"); }
  //--------------------------------------------------------------------------//
  VecI<2> getImageResolution(const nlohmann::json &data) { return getFieldVecI<2>(data, "img_res"); }
  //--------------------------------------------------------------------------//
  bool hasRenderedDomain(const nlohmann::json &data) { return isFieldDomain<2>(data, "rendered_domain"); }
  //--------------------------------------------------------------------------//
  Domain<2> getRenderedDomain(const nlohmann::json &data) { return getFieldDomain<2>(data, "rendered_domain"); }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//