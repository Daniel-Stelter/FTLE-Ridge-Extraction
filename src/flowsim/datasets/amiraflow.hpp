#pragma once
//--------------------------------------------------------------------------//
#include "../flow.hpp"
//--------------------------------------------------------------------------//
#include <vc/legacy/grid_d4.hh>
//--------------------------------------------------------------------------//
#include <memory>
#include <optional>
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  template <uint n>
  class AmiraFlow;
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <>
  class AmiraFlow<2> : public Flow<2>
  {
  public:
    //--------------------------------------------------------------------------//
    AmiraFlow(const std::string &filename,
              const std::string &name = "Amira Dataset",
              bool load_data = true);
    //--------------------------------------------------------------------------//
    virtual ~AmiraFlow();
    //--------------------------------------------------------------------------//
    virtual MaybeVec<2> v(const VecR<2> &pos, real t) const override;
    //--------------------------------------------------------------------------//
    void loadData();
    //--------------------------------------------------------------------------//
    void dropLoadedData();
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    void loadFile(float **pData);
    //--------------------------------------------------------------------------//
    void loadHeaderInfo(bool verbose = false,
                        bool print_header = false);
    //--------------------------------------------------------------------------//
    const std::string m_filepath;
    bool m_is_loaded;
    bool m_is_steady;
    //--------------------------------------------------------------------------//
    int m_dimX = 0;
    int m_dimY = 0;
    int m_dimT = 0;
    VecR<2> *mp_gridData = nullptr;
    //--------------------------------------------------------------------------//
    typedef VC::mvfields::d2::grid<real, VecR<2>> VecR2Grid;
    typedef VC::mvfields::d2::kernel::cubic_convolution_interpolatory_c1<VecR2Grid> VecR2Kernel;
    typedef VC::mvfields::d2::evaluator<VecR2Kernel> VecR2Field;
    //--------------------------------------------------------------------------//
    typedef VC::mvfields::d3::grid<real, VecR<2>> VecR3Grid;
    typedef VC::mvfields::d3::kernel::cubic_convolution_interpolatory_c1<VecR3Grid> VecR3Kernel;
    typedef VC::mvfields::d3::evaluator<VecR3Kernel> VecR3Field;
    //--------------------------------------------------------------------------//
    // Contains sampled vector grids, sorted to by tau value
    std::unique_ptr<VecR2Grid> mp_vec2Grid;
    std::unique_ptr<VecR3Grid> mp_vec3Grid;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
  using AmiraFlow2D = AmiraFlow<2>;
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <>
  class AmiraFlow<3> : public Flow<3>
  {
  public:
    //--------------------------------------------------------------------------//
    AmiraFlow(const std::string &filename,
              const std::string &name = "Amira Dataset",
              bool load_data = true);
    //--------------------------------------------------------------------------//
    AmiraFlow(const std::vector<std::string> &files,
              const Range &time_range,
              const std::string &name = "Amira Dataset",
              bool load_data = true);
    //--------------------------------------------------------------------------//
    virtual ~AmiraFlow();
    //--------------------------------------------------------------------------//
    /// compute flow
    virtual MaybeVec<3> v(const VecR<3> &pos, real t) const override;
    //--------------------------------------------------------------------------//
    void loadData();
    //--------------------------------------------------------------------------//
    void dropLoadedData();
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    void loadFile(const std::string &filename,
                  float **pData);
    //--------------------------------------------------------------------------//
    void loadHeaderInfo(const std::string &filename,
                        bool verbose = false,
                        bool print_header = false);
    //--------------------------------------------------------------------------//
    bool compareHeaderInfo(const std::string &filename) const;
    //--------------------------------------------------------------------------//
    const std::vector<std::string> m_files;
    bool m_is_loaded;
    bool m_is_steady;
    //--------------------------------------------------------------------------//
    int m_dimX = 0;
    int m_dimY = 0;
    int m_dimZ = 0;
    int m_dimT = 0;
    VecR<3> *mp_gridData = nullptr;
    //--------------------------------------------------------------------------//
    typedef VC::mvfields::d3::grid<real, VecR<3>> VecR3Grid;
    typedef VC::mvfields::d3::kernel::cubic_convolution_interpolatory_c1<VecR3Grid> VecR3Kernel;
    typedef VC::mvfields::d3::evaluator<VecR3Kernel> VecR3Field;
    //--------------------------------------------------------------------------//
    typedef VC::mvfields::d4::grid<real, VecR<3>> VecR4Grid;
    typedef VC::mvfields::d4::kernel::cubic_convolution_interpolatory_c1<VecR4Grid> VecR4Kernel;
    typedef VC::mvfields::d4::evaluator<VecR4Kernel> VecR4Field;
    //--------------------------------------------------------------------------//
    // Contains sampled vector grids, sorted to by tau value
    std::unique_ptr<VecR3Grid> mp_vec3Grid;
    std::unique_ptr<VecR4Grid> mp_vec4Grid;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
  using AmiraFlow3D = AmiraFlow<3>;
}
//--------------------------------------------------------------------------//