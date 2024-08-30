#pragma once
//--------------------------------------------------------------------------//
#include "flow.hpp"
#include "numericdiff.hpp"
#include "scalarfield.hpp"
#include "../utils.hpp"
//--------------------------------------------------------------------------//
#include <chrono>
#include <Eigen/Eigenvalues>
#include <filesystem>
#include <fstream>
#include <iostream>
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  template <uint n>
  class FTLEField : virtual public ScalarFieldBase<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    FTLEField(std::shared_ptr<Flow<n>> p_flow,
              real t,
              real tau,
              const Domain<n> &domain,
              real delta = 0.0001)
        : FTLEField(p_flow, t, tau, domain, VecR<n>(delta)) {}
    //--------------------------------------------------------------------------//
    FTLEField(std::shared_ptr<Flow<n>> p_flow,
              real t,
              real tau,
              const Domain<n> &domain,
              const VecR<n> &deltas)
        : ScFieldBase_t(domain),
          mp_flow(p_flow),
          m_t(t),
          m_tau(tau),
          m_deltas(deltas) {}
    //--------------------------------------------------------------------------//
    virtual ~FTLEField() {}
    //--------------------------------------------------------------------------//
    const Flow<n> &flow() const { return *mp_flow; }
    real t() const { return m_t; }
    real tau() const { return m_tau; }
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const override
    {
      if (m_tau == 0.0)
        return 0.0;
      // calculate centered differences for all n dimensions and save in gradient matrix
      EMat<n, n> grad;
      for (uint i = 0; i < n; ++i)
      {
        VecR<n> pos_m1(pos), pos_p1(pos);
        // try to move sideways, but prevent to step outside of the rim
        pos_m1[i] = std::max(pos_m1[i] - m_deltas[i], this->domain(i).min);
        pos_p1[i] = std::min(pos_p1[i] + m_deltas[i], this->domain(i).max);
        real real_dis_start = pos_p1[i] - pos_m1[i];
        // calculate gradient by start and end differences
        VecR<n> vec_dif_end = mp_flow->integrate(pos_p1, m_t, m_tau) -
                              mp_flow->integrate(pos_m1, m_t, m_tau);
        if (!VC::vecn::isfinitenorm(vec_dif_end))
          return NAN;
        for (uint j = 0; j < n; ++j)
          grad(j, i) = vec_dif_end[j] / real_dis_start;
      }
      // find max eigenvalue and calculate FTLE value
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(grad.transpose() * grad);
      auto max_eig = solver.eigenvalues().maxCoeff();
      return log(sqrt(max_eig)) / abs(m_tau);
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    std::shared_ptr<Flow<n>> mp_flow;
    real m_t, m_tau;
    // deltas determine distance of flow map samples from the actual position
    VecR<n> m_deltas;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class InterpolFTLEField : virtual public InterpolScalarField<n>, public FTLEField<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using Flow_ptr_t = std::shared_ptr<Flow<n>>;
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    using ItplScField_t = InterpolScalarField<n>;
    using FTLEField_t = FTLEField<n>;
    //--------------------------------------------------------------------------//
    using VecGrid_t = CornerDataGrid<n, VecR<n>>;
    using RealGrid_t = CornerDataGrid<n, real>;
    using EigenSolver_t = Eigen::SelfAdjointEigenSolver<EMat<n, n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    InterpolFTLEField(Flow_ptr_t p_flow,
                      real t,
                      real tau,
                      const std::string &file_path)
        : InterpolFTLEField(p_flow, t, tau) // protected constructor
    {
      ItplScField_t::readFile(file_path); // automatically sets domain of ScFieldBase
    }
    //--------------------------------------------------------------------------//
    InterpolFTLEField(Flow_ptr_t p_flow,
                      real t,
                      real tau,
                      const VecI<n> &resolution,
                      const Domain<n> &domain)
        : InterpolFTLEField(p_flow, t, tau, resolution, domain, true) {}
    //--------------------------------------------------------------------------//
    virtual ~InterpolFTLEField() {}
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const override { return ItplScField_t::value(pos); }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    InterpolFTLEField(Flow_ptr_t p_flow,
                      real t,
                      real tau)
        : ScFieldBase_t(),
          ItplScField_t(),
          FTLEField_t(p_flow, t, tau, Domain<n>()) {}
    //--------------------------------------------------------------------------//
    InterpolFTLEField(Flow_ptr_t p_flow,
                      real t,
                      real tau,
                      const VecI<n> &resolution,
                      const Domain<n> &domain,
                      bool execute_init)
        : ScFieldBase_t(domain),
          ItplScField_t(resolution, domain),
          FTLEField_t(p_flow, t, tau, domain)
    {
      if (execute_init)
        init(p_flow);
    }
    //--------------------------------------------------------------------------//
    void init(Flow_ptr_t p_flow)
    {
      std::cout << "Init InterpolFTLEField" << std::endl;

      this->mp_flow = p_flow;

      utils::RuntimeAnalyzer analyzer("Integration", "FTLE");

      // start analytics for flow integration
      analyzer.startTimer("Integration");
      VecGrid_t par_grid = setupParticleGrid();
      integrateParticleGrid(par_grid, this->t(), this->tau());
      analyzer.stopTimer("Integration");

      // use positions to calculate FTLE values
      // use one resource per thread for better efficiency
      EMat<n, n> gradient;
      EigenSolver_t solver;
      analyzer.startTimer("FTLE");
#pragma omp parallel for schedule(dynamic) private(gradient, solver)
      for (uint corner_idx = 0; corner_idx < par_grid.numCorners(); ++corner_idx)
      {
        this->m_data[corner_idx] = ftleFromParticleGrid(par_grid,
                                                        corner_idx,
                                                        gradient,
                                                        solver);
      }
      analyzer.stopTimer("FTLE");
      analyzer.logAnalysis();
    }
    //--------------------------------------------------------------------------//
    /** Creates a grid which contains the positions of particles. Each corner
     *  contains a particle with exactly the same position of the corner itself.
     *  Use the integrateParticleGrid() method to advect the contained particles. */
    VecGrid_t setupParticleGrid() const
    {
      // setup grid with same layout and nD vectors at its corners
      VecGrid_t par_grid(this->cellResolution(), this->domain());
#pragma omp parallel for schedule(dynamic)
      for (uint corner_idx = 0; corner_idx < par_grid.numCorners(); ++corner_idx)
        par_grid.cornerData(corner_idx) = par_grid.cornerIndexToPos(corner_idx);
      return par_grid;
    }
    //--------------------------------------------------------------------------//
    /** Integrates all particles with start time t and integration time.
     *  Positions inside the container are overwritten. Optionally, the initial
     *  step sizes can be provided for each particle individually. In this case,
     *  the values are overwritten to the last step size of the integration. */
    void integrateParticleGrid(VecGrid_t &par_grid,
                               real t,
                               real tau,
                               RealGrid_t *p_init_step_grid = nullptr) const
    {
      // test whether p_init_step_grid has valid layout
      assert(!p_init_step_grid || p_init_step_grid->cornerResolution() == par_grid.cornerResolution());

#pragma omp parallel for schedule(dynamic)
      for (uint corner_idx = 0; corner_idx < par_grid.numCorners(); ++corner_idx)
      {
        VecR<n> &pos = par_grid.cornerData(corner_idx);
        if (!VC::vecn::isfinitenorm(pos))
          continue;
        // determine inital step size and setup solver
        real init_step = p_init_step_grid ? p_init_step_grid->cornerData(corner_idx) : 0.01;
        auto solver = ODE<n>::solver(VC::odeint::RK43, ODE<n>::make_options(VC::odeint::AbsTol = 1e-4,
                                                                            VC::odeint::RelTol = 1e-4,
                                                                            VC::odeint::InitialStep = init_step,
                                                                            VC::odeint::MaxStep = 1e-1));
        // start integration at position of the vector
        solver.integrate([this](real t, const VecR<n> &pos) -> MaybeVec<n>
                         { return this->mp_flow->v(pos, t); },
                         t, t + tau, pos);
        // overwrite position to result of the integration
        pos = solver.state == VC::odeint::AcceptStep ? solver.y : VecR<n>(NAN);
        // if provided, set init step
        if (p_init_step_grid)
          p_init_step_grid->cornerData(corner_idx) = solver.absh;
      }
    }
    //--------------------------------------------------------------------------//
    real ftleFromParticleGrid(const VecGrid_t &par_grid,
                              uint corner_idx,
                              EMat<n, n> &grad,
                              EigenSolver_t &solver) const
    {
      VecI<n> corner_coord = this->cornerIndexToCoord(corner_idx);
      bool is_finite = true;
      // go through all dimensions
      for (uint i = 0; i < n; ++i)
      {
        // calculate gradient with four neighbors
        VecR<n> dif(0.0);
        VecI<n> vec_i(0u);
        vec_i[i] = 1;
        // FORWARD DIFFERENCES
        if (corner_coord[i] < 2)
        {
          dif -= par_grid.cornerData(corner_coord) * (3.0 / 2.0);
          dif += par_grid.cornerData(corner_coord + vec_i) * 2.0;
          dif -= par_grid.cornerData(corner_coord + vec_i * 2) * (1.0 / 2.0);
        }
        // BACKWARD DIFFERENCES
        else if (corner_coord[i] > this->cornerResolution(i) - 3)
        {
          dif += par_grid.cornerData(corner_coord) * (3.0 / 2.0);
          dif -= par_grid.cornerData(corner_coord - vec_i) * 2.0;
          dif += par_grid.cornerData(corner_coord - vec_i * 2) * (1.0 / 2.0);
        }
        // CENTRAL DIFFERENCES
        else
        {
          dif += par_grid.cornerData(corner_coord - vec_i * 2) * (1.0 / 12.0);
          dif -= par_grid.cornerData(corner_coord - vec_i) * (2.0 / 3.0);
          dif += par_grid.cornerData(corner_coord + vec_i) * (2.0 / 3.0);
          dif -= par_grid.cornerData(corner_coord + vec_i * 2) * (1.0 / 12.0);
        }
        if (!VC::vecn::isfinitenorm(dif))
        {
          is_finite = false;
          break;
        }
        for (uint j = 0; j < n; ++j)
          grad(j, i) = dif[j] * this->m_rec_edge_len[i];
      }
      if (is_finite)
      {
        // find max eigenvalue and calculate FTLE value
        solver.computeDirect(grad.transpose() * grad);
        auto max_eig = solver.eigenvalues().maxCoeff();
        return log(sqrt(max_eig)) / this->m_tau;
      }
      else
        return NAN;
    }
    //--------------------------------------------------------------------------//
    real ftleFromFlowMap(const VecGrid_t &par_grid,
                         uint corner_idx) const
    {
      EMat<n, n> grad;
      EigenSolver_t solver;
      return ftleFromParticleGrid(par_grid, corner_idx, grad, solver);
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class InterpolFTLEFieldSeries : public InterpolScalarFieldSeries<n>, public InterpolFTLEField<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using Flow_ptr_t = std::shared_ptr<Flow<n>>;
    using FTLEField_t = FTLEField<n>;
    //--------------------------------------------------------------------------//
    using ScFieldBase_t = ScalarFieldBase<n>;
    using ItplScField_t = InterpolScalarField<n>;
    using ItplScFieldSeries_t = InterpolScalarFieldSeries<n>;
    using ItplFTLEField_t = InterpolFTLEField<n>;
    //--------------------------------------------------------------------------//
    using VecGrid_t = CornerDataGrid<n, VecR<n>>;
    using RealGrid_t = CornerDataGrid<n, real>;
    using EigenSolver_t = Eigen::SelfAdjointEigenSolver<EMat<n, n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    InterpolFTLEFieldSeries(Flow_ptr_t p_flow,
                            real t,
                            const std::string &save_dir)
        : ScFieldBase_t(),
          ItplScField_t(),
          ItplScFieldSeries_t(save_dir),  // loads series.info and step 0 data
          ItplFTLEField_t(p_flow, t, 0.0) // tau is port if ItplScFieldSeries
    {
      // time of ItplScFieldSeries corresponds to tau
      FTLEField_t::m_tau = ItplScFieldSeries_t::curTime();
    }
    //--------------------------------------------------------------------------//
    InterpolFTLEFieldSeries(Flow_ptr_t p_flow,
                            real t,
                            const SteppedRange &stepped_tau,
                            const VecI<n> &resolution,
                            const Domain<n> &domain,
                            const std::string &save_dir)
        : ScFieldBase_t(domain),
          ItplScField_t(resolution, domain),
          ItplScFieldSeries_t(stepped_tau, resolution, domain, save_dir),
          ItplFTLEField_t(p_flow, t, stepped_tau.value(0), resolution, domain, false)
    {
      initAll(p_flow);
    }
    //--------------------------------------------------------------------------//
    virtual ~InterpolFTLEFieldSeries() {}
    //--------------------------------------------------------------------------//
    virtual real value(const VecR<n> &pos) const override { return ItplScField_t::value(pos); }
    //--------------------------------------------------------------------------//
    virtual void loadStep(uint step) override
    {
      ItplScFieldSeries_t::loadStep(step);
      FTLEField_t::m_tau = ItplScFieldSeries_t::curTime();
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    void initAll(Flow_ptr_t p_flow)
    {
      this->mp_flow = p_flow;

      // setup grids with particle positions and initial step sizes for flow integration
      VecGrid_t par_grid = ItplFTLEField_t::setupParticleGrid();
      RealGrid_t step_size_grid(this->cellResolution(), this->domain(), 0.01);

      // iterate over all integration time steps
      uint &cur_step = ItplScFieldSeries_t::m_cur_step;
      const uint num_steps = this->m_stepped_time.steps;
      real integ_start_t = this->t();
      const real delta_time = this->m_stepped_time.delta();

      utils::LoopPrinter printer("Init InterpolFTLEFieldSeries", num_steps, false);
      utils::RuntimeAnalyzer analyzer("Integration", "FTLE");
      for (cur_step = 0; cur_step < num_steps; ++cur_step)
      {
        // m_tau MUST BE SET for correct FTLE value computation in parent class
        FTLEField_t::m_tau = ItplScFieldSeries_t::curTime();

        // start analytics for flow integration, up to next
        analyzer.startTimer("Integration");
        if (cur_step == 0)
          ItplFTLEField_t::integrateParticleGrid(par_grid, integ_start_t, this->tau(), &step_size_grid);
        else
          ItplFTLEField_t::integrateParticleGrid(par_grid, integ_start_t, delta_time, &step_size_grid);
        analyzer.stopTimer("Integration");

        // use positions to calculate FTLE values
        // predefine one resource per thread for better efficiency
        EMat<n, n> gradient;
        EigenSolver_t solver;
        analyzer.startTimer("FTLE");
#pragma omp parallel for schedule(dynamic) private(gradient, solver)
        for (uint corner_idx = 0; corner_idx < par_grid.numCorners(); ++corner_idx)
        {
          this->m_data[corner_idx] = ItplFTLEField_t::ftleFromParticleGrid(par_grid,
                                                                           corner_idx,
                                                                           gradient,
                                                                           solver);
        }

        analyzer.stopTimer("FTLE");
        ItplScField_t::writeFile(this->curFilePath());

        integ_start_t += delta_time;

        printer.endIteration();
      }
      analyzer.logAnalysis();

      // after initializing all steps, write info file and load step 0
      this->writeSeriesInfoFile();
      this->loadStep(0);
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//