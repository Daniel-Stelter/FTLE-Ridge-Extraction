#include "amiraflow.hpp"
//--------------------------------------------------------------------------//
#include "../../utils.hpp"
//--------------------------------------------------------------------------//
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdio.h>
//--------------------------------------------------------------------------//
using namespace std;
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  AmiraFlow<2>::AmiraFlow(const string &filepath,
                          const string &name,
                          bool load_data)
      : Flow<2>({Range(), Range()}, Range(), name),
        m_filepath(filepath),
        m_is_loaded(false),
        m_is_steady(false)
  {
    loadHeaderInfo(true);
    if (load_data)
      loadData();
  }
  //--------------------------------------------------------------------------//
  AmiraFlow<2>::~AmiraFlow(void) { delete[] mp_gridData; }
  //--------------------------------------------------------------------------//
  MaybeVec<2> AmiraFlow<2>::v(const VecR<2> &pos, real t) const
  {
    assert(m_is_loaded);
    if (!this->m_space_domain_defined.isInside(pos) || !this->m_time_domain.isInside(t))
      return VC::odeint::OutOfDomain;
    VecR<2> v(0.);
    if (m_is_steady)
    {
      VecR2Field tempField(*mp_vec2Grid);
      tempField.value(v, {pos[0], pos[1]}); // ASDF here we get a warning
    }
    else
    {
      VecR3Field tempField(*mp_vec3Grid);
      tempField.value(v, {pos[0], pos[1], t}); // ASDF here we get a warning
    }
    return v;
  }
  //--------------------------------------------------------------------------//
  void AmiraFlow<2>::loadData()
  {
    if (m_is_loaded)
      return;
    mp_gridData = new VecR<2>[m_dimX * m_dimY * m_dimT];

    // load the data from file
    cout << "AmiraFlow: load data" << endl;
    float *data;
    loadFile(&data);

    // convert the float data into a grid
    int data_id = 0;
    for (int k = 0; k < m_dimT; k++)
      for (int j = 0; j < m_dimY; j++)
        for (int i = 0; i < m_dimX; i++)
        {
          float x = data[data_id * 2];
          float y = data[data_id * 2 + 1];
          VecR<2> gridVec(x, y);
          mp_gridData[data_id++] = gridVec;
        }
    delete[] data;

    real x_min = this->m_space_domain_defined[0].min;
    real x_max = this->m_space_domain_defined[0].max;
    real y_min = this->m_space_domain_defined[1].min;
    real y_max = this->m_space_domain_defined[1].max;
    real t_min = this->m_time_domain.min;
    real t_max = this->m_time_domain.max;
    if (m_is_steady)
    {
      VecR2Grid tempGrid(mp_gridData,
                         VC::mvfields::d2::msize_t(m_dimX, m_dimY));
      tempGrid.set_origin(VecR2Grid::coord_t(x_min, y_min));
      tempGrid.set_extent(VecR2Grid::coord_t(x_max - x_min, y_max - y_min));
      mp_vec2Grid = make_unique<VecR2Grid>(tempGrid);
    }
    else
    {
      VecR3Grid tempGrid(mp_gridData,
                         VC::mvfields::d3::msize_t(m_dimX, m_dimY, m_dimT));
      tempGrid.set_origin(VecR3Grid::coord_t(x_min, y_min, t_min));
      tempGrid.set_extent(VecR3Grid::coord_t(x_max - x_min, y_max - y_min, t_max - t_min));
      mp_vec3Grid = make_unique<VecR3Grid>(tempGrid);
    }
    m_is_loaded = true;
  }
  // ------------------------------------------------------------------------ //
  void AmiraFlow<2>::dropLoadedData()
  {
    cout << "AmiraFlow: drop data" << endl;
    m_is_loaded = false;
    mp_vec2Grid.release();
    mp_vec3Grid.release();
    delete[] mp_gridData;
    mp_gridData = nullptr;
  }
  //--------------------------------------------------------------------------//
  void AmiraFlow<2>::loadFile(float **p_data)
  {
    ifstream in(m_filepath, ios::in | ios::binary);
    if (!in)
      throw runtime_error("AmiraFlow: could not read from file '" + m_filepath + "'");

    char *raw_header_data = new char[2048];
    in.read(raw_header_data, sizeof(char) * 2048);
    string header = string(raw_header_data, 2048);
    delete[] raw_header_data;

    // find the beginning of the data section
    string substring = "# Data section follows";
    size_t string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find data section in file '" + m_filepath + "'");
    }
    in.seekg(string_pos, ios::beg);
    // consume the next two lines which area "# Data section follows" and "@1"
    string token;
    getline(in, token);
    getline(in, token);
    // read the data
    // - how much to read
    const size_t num_to_read = m_dimX * m_dimY * m_dimT * 2;
    // - prepare memory
    *p_data = new float[num_to_read];
    // - do it
    in.read((char *)*p_data, sizeof(float) * num_to_read);
    in.close();
  }
  //--------------------------------------------------------------------------//
  void AmiraFlow<2>::loadHeaderInfo(bool verbose,
                                    bool print_header)
  {
    ifstream in(m_filepath, ios::in | ios::binary);
    if (!in)
      throw runtime_error("AmiraFlow: could not read from file '" + m_filepath + "'");

    char *raw_header_data = new char[2048];
    in.read(raw_header_data, sizeof(char) * 2048);
    string header = string(raw_header_data, 2048);
    delete[] raw_header_data;

    if (print_header)
    {
      utils::printSeparator();
      cout << "HEADER:\n"
           << header << "\n";
      utils::printSeparator();
    }

    string substring = "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1";
    if (header.find(substring) == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: '" + m_filepath + "' is not a proper Amira file");
    }

    // Find the Lattice definition, i.e., the dimensions of the uniform grid
    substring = "define Lattice ";
    size_t string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find lattice information in file '" + m_filepath + "'");
    }
    string_pos += substring.length();
    in.seekg(string_pos, ios::beg);
    string token;
    // clang-format off
    getline(in, token, ' '); m_dimX = atoi(&token[0]);
    getline(in, token, ' '); m_dimY = atoi(&token[0]);
    getline(in, token, ' '); m_dimT = atoi(&token[0]);
    m_is_steady = m_dimT == 0;
    if (m_is_steady)
      m_dimT = 1;
    // clang-format on
    if (verbose)
    {
      cout << "Extracted data:\nLattice is: "
           << m_dimX << " " << m_dimY;
      if (m_is_steady)
        cout << "\nDataset is steady!";
      else
        cout << " " << m_dimT;
      cout << endl;
    }

    // Find the bounding box information
    substring = "BoundingBox ";
    string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find bounding box in file '" + m_filepath + "'");
    }
    string_pos += substring.length();
    in.seekg(string_pos, ios::beg);
    float x_min{1}, y_min{1}, t_min{1}, x_max{-1}, y_max{-1}, t_max{-1};
    // clang-format off
    getline(in, token, ' ');   x_min = atof(&token[0]);
    getline(in, token, ' ');   x_max = atof(&token[0]);
    getline(in, token, ' ');   y_min = atof(&token[0]);
    if (!m_is_steady)
    {
      getline(in, token, ' '); y_max = atof(&token[0]);
      getline(in, token, ' '); t_min = atof(&token[0]);
      getline(in, token);      t_max = atof(&token[0]);

    }
    else
    {
      getline(in, token);      y_max = atof(&token[0]);
    }
    // clang-format on
    this->m_space_domain_defined = this->m_space_domain_standard = {Range(x_min, x_max), Range(y_min, y_max)};
    if (!m_is_steady)
      this->m_time_domain = Range(t_min, t_max);
    if (verbose)
    {
      cout << "Bounding box in x-direction: " << this->m_space_domain_defined[0]
           << "\nBounding box in y-direction: " << this->m_space_domain_defined[1]
           << "\nBounding box in t-direction: " << this->m_time_domain
           << endl;
    }

    // check if we have uniform data
    substring = "CoordType \"uniform\"";
    string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find 'CoordType \"uniform\"' in file '" + m_filepath + "'");
    }

    // type of the field: scalar, vector
    int num_components = {0};
    substring = "Lattice { float Data }";
    string_pos = header.find(substring);
    if (string_pos != string::npos)
      num_components = 1;
    substring = "Lattice { float[";
    string_pos = header.find(substring);
    if (string_pos != string::npos)
    {
      string_pos += substring.length();
      in.seekg(string_pos, ios::beg);
      getline(in, token, ']');
      num_components = atoi(&token[0]);
    }

    // sanity check
    if (m_dimX <= 0 || m_dimY <= 0 || m_dimT <= 0 ||
        x_min > x_max || y_min > y_max ||
        (!m_is_steady && t_min > t_max) ||
        num_components != 2)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not load reasonable domain in file '" + m_filepath + "'");
    }
    if (verbose)
      cout << "AmiraFlow: successfully loaded info" << endl;
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  AmiraFlow<3>::AmiraFlow(const string &filepath,
                          const string &name,
                          bool load_data)
      : AmiraFlow({filepath}, Range(), name, load_data) {}
  //--------------------------------------------------------------------------//
  AmiraFlow<3>::AmiraFlow(const vector<string> &files,
                          const Range &time_range,
                          const string &name,
                          bool load_data)
      : Flow<3>({Range(), Range(), Range()}, time_range, name),
        m_files(files),
        m_is_loaded(false),
        m_is_steady(files.size() < 2),
        m_dimT(files.size())
  {
    // load one header info and compare all others
    loadHeaderInfo(m_files.front(), true);
    for (uint fid = 1; fid < m_files.size(); ++fid)
      if (compareHeaderInfo(m_files[fid]))
        throw runtime_error("AmiraFlow: headers of input files are inconsistent");
    if (load_data)
      loadData();
  }
  //--------------------------------------------------------------------------//
  AmiraFlow<3>::~AmiraFlow(void) { delete[] mp_gridData; }
  //--------------------------------------------------------------------------//
  MaybeVec<3> AmiraFlow<3>::v(const VecR<3> &pos, real t) const
  {
    assert(m_is_loaded);
    if (!this->m_space_domain_defined.isInside(pos) || !this->m_time_domain.isInside(t))
      return VC::odeint::OutOfDomain;
    VecR<3> v;
    if (m_is_steady)
    {
      VecR3Field tempField(*mp_vec3Grid);
      tempField.value(v, {pos[0], pos[1], pos[2]}); // ASDF here we get a warning
    }
    else
    {
      VecR4Field tempField(*mp_vec4Grid);
      tempField.value(v, {pos[0], pos[1], pos[2], t}); // ASDF here we get a warning
    }
    return v;
  }
  //--------------------------------------------------------------------------//
  void AmiraFlow<3>::loadData()
  {
    if (m_is_loaded)
      return;
    mp_gridData = new VecR<3>[m_dimX * m_dimY * m_dimZ * m_dimT];

    int grid_id = 0;
    utils::LoopPrinter loop_printer("AmiraFlow: load data", m_files.size());
    for (uint fid = 0; fid < m_files.size(); ++fid)
    {
      float *data;
      loadFile(m_files[fid], &data);

      // convert the float data into a grid
      int data_id = 0;
      for (int k = 0; k < m_dimZ; k++)
        for (int j = 0; j < m_dimY; j++)
          for (int i = 0; i < m_dimX; i++)
          {
            float x = data[data_id++];
            float y = data[data_id++];
            float z = data[data_id++];
            VecR<3> gridVec(x, y, z);
            mp_gridData[grid_id++] = gridVec;
          }
      delete[] data;
      loop_printer.endIteration();
    }

    real x_min = this->m_space_domain_defined[0].min;
    real x_max = this->m_space_domain_defined[0].max;
    real y_min = this->m_space_domain_defined[1].min;
    real y_max = this->m_space_domain_defined[1].max;
    real z_min = this->m_space_domain_defined[2].min;
    real z_max = this->m_space_domain_defined[2].max;
    real t_min = this->m_time_domain.min;
    real t_max = this->m_time_domain.max;
    if (m_is_steady)
    {
      VecR3Grid tempGrid(mp_gridData,
                         VC::mvfields::d3::msize_t(m_dimX, m_dimY, m_dimZ));
      tempGrid.set_origin(VecR3Grid::coord_t(x_min, y_min, z_min));
      tempGrid.set_extent(VecR3Grid::coord_t(x_max - x_min, y_max - y_min, z_max - z_min));
      mp_vec3Grid = make_unique<VecR3Grid>(tempGrid);
    }
    else
    {
      VecR4Grid tempGrid(mp_gridData,
                         VC::mvfields::d4::msize_t(m_dimX, m_dimY, m_dimZ, m_dimT));
      tempGrid.set_origin(VecR4Grid::coord_t(x_min, y_min, z_min, t_min));
      tempGrid.set_extent(VecR4Grid::coord_t(x_max - x_min, y_max - y_min, z_max - z_min, t_max - t_min));
      mp_vec4Grid = make_unique<VecR4Grid>(tempGrid);
    }
    m_is_loaded = true;
  }
  // ------------------------------------------------------------------------ //
  void AmiraFlow<3>::dropLoadedData()
  {
    cout << "AmiraFlow: drop data" << endl;
    m_is_loaded = false;
    mp_vec3Grid.release();
    mp_vec4Grid.release();
    delete[] mp_gridData;
    mp_gridData = nullptr;
  }
  //--------------------------------------------------------------------------//
  void AmiraFlow<3>::loadFile(const string &filename,
                              float **p_data)
  {
    ifstream in(filename, ios::in | ios::binary);
    if (!in)
      throw runtime_error("AmiraFlow: could not read from file '" + filename + "'");

    char *raw_header_data = new char[2048];
    in.read(raw_header_data, sizeof(char) * 2048);
    string header = string(raw_header_data, 2048);
    delete[] raw_header_data;

    // find the beginning of the data section
    string substring = "# Data section follows";
    size_t string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find data section in file '" + filename + "'");
    }
    in.seekg(string_pos, ios::beg);
    // consume the next two lines which area "# Data section follows" and "@1"
    string token;
    getline(in, token);
    getline(in, token);
    // read the data
    // - how much to read
    const size_t num_to_read = m_dimX * m_dimY * m_dimZ * m_dimT * 3;
    // - prepare memory
    *p_data = new float[num_to_read];
    // - do it
    in.read((char *)*p_data, sizeof(float) * num_to_read);
    in.close();
  }
  //--------------------------------------------------------------------------//
  void AmiraFlow<3>::loadHeaderInfo(const string &filename,
                                    bool verbose,
                                    bool print_header)
  {
    cout << "AmiraFlow: Load infos from header" << endl;
    ifstream in(filename, ios::in | ios::binary);
    if (!in)
      throw runtime_error("AmiraFlow: could not read from file '" + filename + "'");

    char *raw_header_data = new char[2048];
    in.read(raw_header_data, sizeof(char) * 2048);
    string header = string(raw_header_data, 2048);
    delete[] raw_header_data;

    if (print_header)
    {
      utils::printSeparator();
      cout << "HEADER:\n"
           << header << "\n";
      utils::printSeparator();
    }

    string substring = "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1";
    if (header.find(substring) == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: '" + filename + "' is not a proper Amira file");
    }

    // Find the Lattice definition, i.e., the dimensions of the uniform grid
    substring = "define Lattice ";
    size_t string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find lattice information in file '" + filename + "'");
    }
    string_pos += substring.length();
    in.seekg(string_pos, ios::beg);
    string token;
    // clang-format off
    getline(in, token, ' '); m_dimX = atoi(&token[0]);
    getline(in, token, ' '); m_dimY = atoi(&token[0]);
    getline(in, token);      m_dimZ = atoi(&token[0]);
    // clang-format on
    if (verbose)
    {
      cout << "Extracted data:\nLattice is: "
           << m_dimX << " " << m_dimY << " " << m_dimZ;
      if (m_is_steady)
        cout << "\nDataset is steady!";
      else
        cout << " " << m_dimT;
      cout << endl;
    }

    // Find the bounding box information
    substring = "BoundingBox ";
    string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find bounding box in file '" + filename + "'");
    }
    string_pos += substring.length();
    in.seekg(string_pos, ios::beg);
    float x_min{1}, y_min{1}, z_min{1}, x_max{-1}, y_max{-1}, z_max{-1};
    // clang-format off
    getline(in, token, ' '); x_min = atof(&token[0]);
    getline(in, token, ' '); x_max = atof(&token[0]);
    getline(in, token, ' '); y_min = atof(&token[0]);
    getline(in, token, ' '); y_max = atof(&token[0]);
    getline(in, token, ' '); z_min = atof(&token[0]);
    getline(in, token);      z_max = atof(&token[0]);
    // clang-format on
    this->m_space_domain_defined = this->m_space_domain_standard = {Range(x_min, x_max), Range(y_min, y_max), Range(z_min, z_max)};
    if (verbose)
    {
      cout << "Bounding box in x-direction: " << this->m_space_domain_defined[0]
           << "\nBounding box in y-direction: " << this->m_space_domain_defined[1]
           << "\nBounding box in z-direction: " << this->m_space_domain_defined[2]
           << "\nBounding box in t-direction: " << this->m_time_domain
           << endl;
    }

    // check if we have uniform data
    substring = "CoordType \"uniform\"";
    string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find 'CoordType \"uniform\"' in file '" + filename + "'");
    }

    // type of the field: scalar, vector
    int num_components = {0};
    substring = "Lattice { float Data }";
    string_pos = header.find(substring);
    if (string_pos != string::npos)
      num_components = 1;
    substring = "Lattice { float[";
    string_pos = header.find(substring);
    if (string_pos != string::npos)
    {
      string_pos += substring.length();
      in.seekg(string_pos, ios::beg);
      getline(in, token, ']');
      num_components = atoi(&token[0]);
    }

    // sanity check
    if (m_dimX <= 0 || m_dimY <= 0 || m_dimZ <= 0 || m_dimT <= 0 ||
        x_min > x_max || y_min > y_max || z_min > z_max ||
        num_components != 3)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not load reasonable domain in file '" + filename + "'");
    }
    if (verbose)
      cout << "AmiraFlow: successfully loaded info" << endl;
  }
  //--------------------------------------------------------------------------//
  bool AmiraFlow<3>::compareHeaderInfo(const string &filename) const
  {
    ifstream in(filename, ios::in | ios::binary);
    if (!in)
      throw runtime_error("AmiraFlow: could not read from file '" + filename + "'");

    char *raw_header_data = new char[2048];
    in.read(raw_header_data, sizeof(char) * 2048);
    string header = string(raw_header_data, 2048);
    delete[] raw_header_data;

    // search for file begin
    string substring = "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1";
    if (header.find(substring) == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: '" + filename + "' is not a proper Amira file");
    }

    // Find the Lattice definition, i.e., the dimensions of the uniform grid
    substring = "define Lattice ";
    size_t string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find lattice information in file '" + filename + "'");
    }
    string_pos += substring.length();
    in.seekg(string_pos, ios::beg);
    string token;
    getline(in, token, ' ');
    int dimX = atoi(&token[0]);
    getline(in, token, ' ');
    int dimY = atoi(&token[0]);
    getline(in, token);
    int dimZ = atoi(&token[0]);
    if (dimX != m_dimX || dimY != m_dimY || dimZ != m_dimZ)
    {
      in.close();
      throw runtime_error("AmiraFlow: inconsistent lattice information w.r.t. other previously loaded files in file '" + filename + "'");
    }

    // Find the bounding box information
    substring = "BoundingBox ";
    string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find bounding box in file '" + filename + "'");
    }
    string_pos += substring.length();
    in.seekg(string_pos, ios::beg);
    float x_min{1}, y_min{1}, z_min{1}, x_max{-1}, y_max{-1}, z_max{-1};
    // clang-format off
    getline(in, token, ' '); x_min = atof(&token[0]);
    getline(in, token, ' '); x_max = atof(&token[0]);
    getline(in, token, ' '); y_min = atof(&token[0]);
    getline(in, token, ' '); y_max = atof(&token[0]);
    getline(in, token, ' '); z_min = atof(&token[0]);
    getline(in, token);      z_max = atof(&token[0]);
    // clang-format on
    if (this->m_space_domain_defined[0].min != x_min || this->m_space_domain_defined[0].max != x_max ||
        this->m_space_domain_defined[1].min != y_min || this->m_space_domain_defined[1].max != y_max ||
        this->m_space_domain_defined[2].min != z_min || this->m_space_domain_defined[2].max != z_max)
    {
      in.close();
      throw runtime_error("AmiraFlow: inconsistent bounding box w.r.t. other previously loaded files in file '" + filename + "'");
    }

    // check if we have uniform data
    substring = "CoordType \"uniform\"";
    string_pos = header.find(substring);
    if (string_pos == string::npos)
    {
      in.close();
      throw runtime_error("AmiraFlow: could not find 'CoordType \"uniform\"' in file '" + filename + "'");
    }

    // type of the field: scalar, vector
    int num_components = {0};
    substring = "Lattice { float Data }";
    string_pos = header.find(substring);
    if (string_pos != string::npos)
      num_components = 1;
    substring = "Lattice { float[";
    string_pos = header.find(substring);
    if (string_pos != string::npos)
    {
      string_pos += substring.length();
      in.seekg(string_pos, ios::beg);
      getline(in, token, ']');
      num_components = atoi(&token[0]);
    }
    if (num_components != 3)
    {
      in.close();
      throw runtime_error("AmiraFlow: inconsistent number of components w.r.t. other previously loaded files in file '" + filename + "'");
    }

    return true;
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//