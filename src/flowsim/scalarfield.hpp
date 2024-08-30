#pragma once
//--------------------------------------------------------------------------//
#include "grid.hpp"
#include "../utils.hpp"
//--------------------------------------------------------------------------//
#include <filesystem>
#include <fstream>
#include <iostream>
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  template <uint n>
  class ScalarFieldBase
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ScalarFieldBase(const Domain<n> &domain) : m_domain(domain) {}
    //--------------------------------------------------------------------------//
    virtual ~ScalarFieldBase() {}
    //--------------------------------------------------------------------------//
    virtual const Domain<n> &domain() const { return m_domain; }
    virtual const Range &domain(uint i) const { return m_domain[i]; }
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const = 0;
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScalarFieldBase() : ScalarFieldBase(Domain<n>()) {}
    //--------------------------------------------------------------------------//
    Domain<n> m_domain;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  std::function<real(const VecR<n> &)> asFunction(std::shared_ptr<ScalarFieldBase<n>> &p_sc_field)
  {
    return [p_sc_field](const VecR<n> &pos)
    { return p_sc_field->value(pos); };
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ScalarField : public ScalarFieldBase<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    using Func_t = std::function<real(const VecR<n> &)>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    template <typename F> ScalarField(F&& function, const Domain<n> &domain)
        : ScFieldBase_t(domain),
          m_function(function) {}
    //--------------------------------------------------------------------------//
    virtual ~ScalarField() {}
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const
    {
      assert(m_domain.isInside(pos));
      return m_function(pos);
    }
    //--------------------------------------------------------------------------//
    const Func_t &function() const { return m_function; }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScalarField() : ScalarField([](const VecR<n> &) -> real
                                { return 0; },
                                Domain<n>()) {}
    //--------------------------------------------------------------------------//
    Func_t m_function;
    Domain<n> m_domain;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n, uint m>
  class ReducedScalarField : public ScalarFieldBase<n>
  {
    //--------------------------------------------------------------------------//
    using HigherScalarField_ptr_t = std::shared_ptr<ScalarFieldBase<n + m>>;
    //--------------------------------------------------------------------------//
    HigherScalarField_ptr_t mp_higher_field;
    const VecI<n> m_free_dims;
    const VecI<m> m_fixed_dims;
    VecR<m> m_fixed_vals;
    //--------------------------------------------------------------------------//
    VecR<n + m> higherPos(const VecR<n> &pos) const
    {
      VecR<n + m> higher;
      for (uint j = 0; j < m; ++j)
        higher[m_fixed_dims[j]] = m_fixed_vals[j];
      for (uint i = 0; i < n; ++i)
        higher[m_free_dims[i]] = pos[i];
      return higher;
    }
    //--------------------------------------------------------------------------//
    VecI<n> findFreeDims(const VecI<m> &fixed_dims) const
    {
      VecI<n> free_dims;
      for (uint i = 0, j = 0; i < n; ++i)
      {
        while (j < m && fixed_dims[j] == i + j)
          ++j;
        uint free_i = i + j;
        free_dims[i] = free_i;
      }
      return free_dims;
    }
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ReducedScalarField(HigherScalarField_ptr_t p_higher_field,
                       const VecI<m> &fixed_dims,
                       const VecR<m> &fixed_vals)
        : ScalarFieldBase<n>(Domain<n>()),
          mp_higher_field(p_higher_field),
          m_free_dims(findFreeDims(fixed_dims)),
          m_fixed_dims(fixed_dims),
          m_fixed_vals(fixed_vals)
    {
      for (uint i = 0; i < n; ++i)
        this->m_domain[i] = mp_higher_field->domain(m_free_dims[i]);
    }
    //--------------------------------------------------------------------------//
    virtual ~ReducedScalarField() {}
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const override { return mp_higher_field->value(higherPos(pos)); }
    //--------------------------------------------------------------------------//
    void setFixedVal(uint fixed_dim, real new_fixed_val)
    {
      for (uint i = 0; i < m; ++i)
        if (m_fixed_dims[i] == fixed_dim)
        {
          m_fixed_vals[i] = new_fixed_val;
          return;
        }
      throw std::runtime_error("ReducedScalarField: unable to set fixed value for free dimension");
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ScalarFieldSeriesBase : virtual public ScalarFieldBase<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ScalarFieldSeriesBase(const SteppedRange &stepped_time,
                          const Domain<n> &domain)
        : ScFieldBase_t(domain),
          m_stepped_time(stepped_time),
          m_cur_step(0) {}
    //--------------------------------------------------------------------------//
    virtual ~ScalarFieldSeriesBase() {}
    //--------------------------------------------------------------------------//
    const SteppedRange &steppedTime() const { return m_stepped_time; }
    //--------------------------------------------------------------------------//
    real curTime() const { return m_stepped_time.value(m_cur_step); }
    //--------------------------------------------------------------------------//
    virtual void loadStep(uint step) = 0;
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScalarFieldSeriesBase() : ScFieldBase_t() {}
    //--------------------------------------------------------------------------//
    SteppedRange m_stepped_time;
    uint m_cur_step;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ScalarFieldSeries : public ScalarFieldSeriesBase<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    using ScFieldSeriesBase_t = ScalarFieldSeriesBase<n>;
    //--------------------------------------------------------------------------//
    using TimeDepFunc_t = std::function<real(const VecR<n + 1> &)>;
    //--------------------------------------------------------------------------//
    using TimeDepScField_t = ScalarField<n + 1>;
    using TimeDepScField_ptr_t = std::shared_ptr<TimeDepScField_t>;
    //--------------------------------------------------------------------------//
    using ReducedScField_t = ReducedScalarField<n, 1>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ScalarFieldSeries(const TimeDepFunc_t &time_dep_function,
                      const SteppedRange &stepped_time,
                      const Domain<n> &domain)
        : ScFieldBase_t(domain),
          ScFieldSeriesBase_t(stepped_time, domain),
          m_reduced_sc_field(makeSharedTimeDepScFieldPtr(time_dep_function, stepped_time, domain),
                             {n},
                             {stepped_time.min}) {}
    //--------------------------------------------------------------------------//
    virtual ~ScalarFieldSeries() {}
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const override
    {
      return m_reduced_sc_field.value(pos);
    }
    //--------------------------------------------------------------------------//
    virtual void loadStep(uint step) override
    {
      if (step > this->m_stepped_time.steps)
        throw std::runtime_error("ScalarFieldSeries: try to load step out of range");
      this->m_cur_step = step;
      m_reduced_sc_field.setFixedVal(n, this->curTime());
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScalarFieldSeries(const SteppedRange &stepped_time,
                      const Domain<n> &domain)
        : ScFieldBase_t(domain),
          ScFieldSeriesBase_t(stepped_time, domain) {}
    //--------------------------------------------------------------------------//
    TimeDepScField_ptr_t
    makeSharedTimeDepScFieldPtr(const TimeDepFunc_t &time_dep_function,
                                const SteppedRange &stepped_time,
                                const Domain<n> &spatial_domain)
    {
      Domain<n + 1> full_domain;
      for (uint i = 0; i < n; ++i)
        full_domain[i] = spatial_domain[i];
      full_domain[n].min = stepped_time.min;
      full_domain[n].max = stepped_time.max;
      return std::make_shared<TimeDepScField_t>(time_dep_function, full_domain);
    }
    //--------------------------------------------------------------------------//
    ReducedScField_t m_reduced_sc_field;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class InterpolScalarField : virtual public ScalarFieldBase<n>, public CornerDataGrid<n, real>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using Func_t = std::function<real(const VecR<n> &)>;
    using ScFieldBase_t = ScalarFieldBase<n>;
    using CornerDataGrid_t = CornerDataGrid<n, real>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    InterpolScalarField(const std::string &file_path)
        : ScFieldBase_t(),
          CornerDataGrid_t()
    {
      readFile(file_path);
    }
    //--------------------------------------------------------------------------//
    InterpolScalarField(const Func_t &function,
                        const VecI<n> &cell_resolution,
                        const Domain<n> &domain)
        : InterpolScalarField(cell_resolution, domain)
    {
      init(function);
    }
    //--------------------------------------------------------------------------//
    virtual ~InterpolScalarField() {}
    //--------------------------------------------------------------------------//
    // ScalarFieldBase and CornerDataGrid both implement domain, but as they are
    // consistent, this is no problem.
    using ScFieldBase_t::domain;
    //--------------------------------------------------------------------------//
    virtual void readFile(const std::string &file_path) override
    {
      CornerDataGrid_t::readFile(file_path);
      ScFieldBase_t::m_domain = CornerDataGrid_t::m_domain;
    }
    //--------------------------------------------------------------------------//
    virtual void redefine(const VecI<n> &cell_resolution, const Domain<n> &domain) override
    {
      CornerDataGrid_t::redefine(cell_resolution, domain);
      ScFieldBase_t::m_domain = domain;
    }
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const override
    {
      const VecI<n> cell_coord = this->posToCellCoord(pos);
      const VecR<n> min_corner = this->cellMinPos(cell_coord);

      // find relative lengths from min corner to provided position inside the cell
      VecR<n> fractions;
      for (uint i = 0; i < n; ++i)
        fractions[i] = (pos[i] - min_corner[i]) * this->m_rec_edge_len[i];

      real result = 0.0;
      // calc for each corner (e.g. 4 corners in 2D, 8 corners in 3D ...)
      for (uint bin_id = 0; bin_id < std::pow(2, n); ++bin_id)
      {
        // 'unpack' bin_id to corresponding cell corner
        VecI<n> corner_coord(cell_coord);
        for (uint i = 0, div = 1; i < n; ++i, div *= 2)
          corner_coord[i] += (bin_id / div % 2); // depending on binId: add one or not
        // calc influence: small fraction -> near to corner -> great influence
        real influence = 1;
        for (uint i = 0; i < n; ++i)
          influence *= (corner_coord[i] == cell_coord[i]) ? 1 - fractions[i] : fractions[i];
        result += influence * this->cornerData(corner_coord);
      }
      return result;
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    InterpolScalarField() : ScFieldBase_t(), CornerDataGrid_t() {}
    //--------------------------------------------------------------------------//
    InterpolScalarField(const VecI<n> &cell_resolution,
                        const Domain<n> &domain)
        : ScFieldBase_t(domain),
          CornerDataGrid_t(cell_resolution, domain) {}
    //--------------------------------------------------------------------------//
    void init(const Func_t &function)
    {
      std::cout << "Init InterpolScalarField" << std::endl;
#pragma omp parallel for schedule(dynamic)
      for (uint corner_idx = 0; corner_idx < this->numCorners(); ++corner_idx)
      {
        const VecR<n> pos = this->cornerIndexToPos(corner_idx);
        this->m_data[corner_idx] = function(pos);
      }
      utils::printSeparator();
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class InterpolScalarFieldSeries : virtual public InterpolScalarField<n>, public ScalarFieldSeriesBase<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using TimeDepFunc_t = std::function<real(const VecR<n + 1> &)>;
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    using ItplScField_t = InterpolScalarField<n>;
    using ScFieldSeriesBase_t = ScalarFieldSeriesBase<n>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    static bool existsSave(const std::string &save_dir)
    {
      return std::filesystem::exists(save_dir + "./series.info");
    }
    //--------------------------------------------------------------------------//
    InterpolScalarFieldSeries(const std::string &save_dir)
        : ScFieldBase_t(),
          ItplScField_t(),
          ScFieldSeriesBase_t(),
          m_save_dir(save_dir)
    {
      readSeriesInfoFile();
      loadStep(0);
    }
    //--------------------------------------------------------------------------//
    InterpolScalarFieldSeries(const TimeDepFunc_t &time_dep_function,
                              const SteppedRange &stepped_time,
                              const VecI<n> &cell_resolution,
                              const Domain<n> &domain,
                              const std::string &save_dir)
        : ScFieldBase_t(domain),
          ItplScField_t(cell_resolution, domain),
          ScFieldSeriesBase_t(stepped_time, domain),
          m_save_dir(save_dir)
    {
      initAll(time_dep_function);
    }
    //--------------------------------------------------------------------------//
    virtual ~InterpolScalarFieldSeries() {}
    //--------------------------------------------------------------------------//
    virtual void loadStep(uint step) override
    {
      if (step >= this->m_stepped_time.steps)
        throw std::runtime_error("InterpolScalarFieldSeries: try to load step out of range");
      this->m_cur_step = step;
      Domain<n> old_domain = this->domain();
      VecI<n> old_res = this->cellResolution();
      this->readFile(curFilePath());
      if (old_domain != this->domain() || old_res != this->cellResolution())
        throw std::runtime_error("InterpolScalarFieldSeries: loaded file has inconsistent domain or resolution");
    }
    //--------------------------------------------------------------------------//
    bool hasSameParams(const SteppedRange &stepped_time,
                       const VecI<n> &cell_resolution,
                       const Domain<n> &domain) const
    {
      return this->steppedTime() == stepped_time &&
             this->cellResolution() == cell_resolution &&
             this->domain() == domain;
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    InterpolScalarFieldSeries(const SteppedRange &stepped_time,
                              const VecI<n> &cell_resolution,
                              const Domain<n> &domain,
                              const std::string &save_dir)
        : ScFieldBase_t(domain),
          ItplScField_t(cell_resolution, domain),
          ScFieldSeriesBase_t(stepped_time, domain),
          m_save_dir(save_dir) {}
    //--------------------------------------------------------------------------//
    void initAll(const TimeDepFunc_t &time_dep_function)
    {
      // create reduced scalar field from time_dep_function
      Domain<n + 1> time_dep_domain;
      for (uint i = 0; i < n; ++i)
        time_dep_domain[i] = this->domain(i);
      time_dep_domain[n] = {this->m_stepped_time.min, this->m_stepped_time.max};
      std::shared_ptr<ScalarField<n + 1>> p_time_dep_sc_field =
          std::make_shared<ScalarField<n + 1>>(time_dep_function, time_dep_domain);
      ReducedScalarField<n, 1> reduced_sc_field(p_time_dep_sc_field, {n}, {this->curTime()});

      // iterate over all time steps
      uint &cur_step = this->m_cur_step;
      const uint num_steps = this->m_stepped_time.steps;
      utils::LoopPrinter printer("Init InterpolScalarFieldSeries", num_steps);
      for (cur_step = 0; cur_step < num_steps; ++cur_step)
      {
        reduced_sc_field.setFixedVal(n, this->curTime());
        this->init([&reduced_sc_field](const VecR<n> &pos)
                   { return reduced_sc_field.value(pos); });
        this->writeFile(curFilePath());

        printer.endIteration();
      }
      writeSeriesInfoFile();
      loadStep(0);
    }
    //--------------------------------------------------------------------------//
    std::string curFilePath() const { return m_save_dir + "/" + std::to_string(this->m_cur_step) + ".data"; }
    //--------------------------------------------------------------------------//
    void writeSeriesInfoFile() const
    {
      using namespace std;

      string file_path = m_save_dir + "/series.info";
      filesystem::create_directories(filesystem::absolute(file_path).parent_path());
      ofstream out(file_path, ios::binary);
      if (!out)
        throw runtime_error("InterpolScalarFieldSeries: could not write to file '" + file_path + "'");

      // write dimension, cell resolution, domain and stepped time
      uint dim = n;
      out.write((char *)&dim, sizeof(dim));
      out.write((char *)&this->cellResolution(), sizeof(this->cellResolution()));
      out.write((char *)&this->domain(), sizeof(this->domain()));
      out.write((char *)&this->steppedTime(), sizeof(this->steppedTime()));

      out.close();
    }
    //--------------------------------------------------------------------------//
    void readSeriesInfoFile()
    {
      using namespace std;

      string file_path = m_save_dir + "/series.info";
      ifstream in(file_path, ios::binary);
      if (!in)
        throw runtime_error("InterpolScalarFieldSeries: could not read from file '" + file_path + "'");

      // check dimension
      uint dim;
      in.read((char *)&dim, sizeof(dim));
      if (dim != n)
        throw runtime_error("InterpolScalarFieldSeries: expected dimension " + to_string(n) + ", but was " + to_string(dim));
      // get cell resolution and domain, then redefine
      VecI<n> cell_res;
      Domain<n> domain;
      in.read((char *)&cell_res, sizeof(cell_res));
      in.read((char *)&domain, sizeof(domain));
      this->redefine(cell_res, domain);
      // read stepped time
      in.read((char *)&this->steppedTime(), sizeof(this->steppedTime()));
    }
    //--------------------------------------------------------------------------//
    const std::string m_save_dir;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//
