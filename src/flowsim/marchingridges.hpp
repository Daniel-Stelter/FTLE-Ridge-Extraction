#pragma once
//--------------------------------------------------------------------------//
#include "numericdiff.hpp"
#include "../params.hpp"
//--------------------------------------------------------------------------//
#include <optional>
#include <map>
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  class MarchingRidges3D
  {
    //--------------------------------------------------------------------------//
    using ItplScField_ptr_t = std::shared_ptr<InterpolScalarField<3>>;
    using ScFieldBase_ptr_t = std::shared_ptr<ScalarFieldBase<3>>;
    using Params_t = params::ParamsMarRidges<3>;
    //--------------------------------------------------------------------------//
    using VertIndex = uint;
    using EdgeIndex = uint;
    using TriangleIndex = uint;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    MarchingRidges3D(ItplScField_ptr_t p_itpl_sc_field,
                     real min_ridge_strength);
    //--------------------------------------------------------------------------//
    MarchingRidges3D(ItplScField_ptr_t p_itpl_sc_field,
                     const Params_t &params);
    //--------------------------------------------------------------------------//
    MarchingRidges3D(ScFieldBase_ptr_t p_sc_field_base,
                     real min_ridge_strength,
                     const VecI<3> &cell_resolution);
    //--------------------------------------------------------------------------//
    virtual ~MarchingRidges3D() {}
    //--------------------------------------------------------------------------//
    /**
     * @brief Uses the provided scalar field to execute an extensive search for ridge
     * points as intersections of the ridge surface with edges of a regular grid.
     * The points are triangulated and saved internally.
     */
    void calcRidgeSurface();
    //--------------------------------------------------------------------------//
    /**
     * @brief Saves calculated triangulation to the provided file. Also saves scalar
     * field and ridge strength values per vertex.
     * @param file_path
     */
    void saveRidgeSurface(const std::string &file_path) const;
    //--------------------------------------------------------------------------//
  private:
    //--------------------------------------------------------------------------//
    /**
     * @brief Tries to find a ridge intersection in the provided edge between pos1
     * and pos2 using pre-calculated local information.
     * @param pos1
     * @param grad1
     * @param pos2
     * @param grad2
     * @param avg_transverse
     * @return optional<VecR<3>>
     */
    std::optional<VecR<3>> calcLineIntersection(const VecR<3> &pos1,
                                                real grad1,
                                                const VecR<3> &pos2,
                                                real grad2,
                                                const VecR<3> &avg_transverse) const;
    //--------------------------------------------------------------------------//
    /**
     * @brief Searches cell for intersections with a ridge surface. If possible,
     * creates and saves the corresponding triangles.
     * @param cell_coord
     */
    void findCellTriangles(const VecI<3> &cell_coord);
    //--------------------------------------------------------------------------//
    /// @brief Scalar field (classically interpolated) which is used for the
    /// ridge search.
    ScFieldBase_ptr_t mp_sc_field_base;
    /// @brief A grid which defines the cells for the ridge search. In case of an
    /// InterpolScalarField, the parameters of both grids coincide.
    Grid<3> m_grid;
    /// @brief Numeric differentiation of the provided scalar field.
    SmoothNumericDiff<3> m_num_diff;
    /// @brief Minimum required ridge strength for an intersection to be accepted.
    real m_min_ridge_strength;
    /// @brief Container for all vertices of the triangle mesh.
    std::vector<VecR<3>> m_vertices;
    /// @brief Mapping from edge indices of the grid to vertex indices of m_vertices.
    /// It helps for efficiently finding already calculated intersections.
    std::map<EdgeIndex, VertIndex> m_edge_to_vertex_map;
    /// @brief List of triples of vertex indices which together represent a single
    /// triangle.
    std::vector<std::array<VertIndex, 3>> m_triangles;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  class MarchingRidgesSeries3D
  {
    //--------------------------------------------------------------------------//
    using ItplScFieldSeries_ptr_t = std::shared_ptr<InterpolScalarFieldSeries<3>>;
    using ScFieldSeriesBase_ptr_t = std::shared_ptr<ScalarFieldSeriesBase<3>>;
    using Params_t = params::ParamsMarRidges<3>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    MarchingRidgesSeries3D(ItplScFieldSeries_ptr_t p_itpl_sc_field_series,
                           real min_ridge_strength)
        : mp_sc_field_series_base(p_itpl_sc_field_series),
          m_mar_ridges(p_itpl_sc_field_series, min_ridge_strength) {}
    //--------------------------------------------------------------------------//
    MarchingRidgesSeries3D(ItplScFieldSeries_ptr_t p_itpl_sc_field_series,
                           const Params_t &params)
        : MarchingRidgesSeries3D(p_itpl_sc_field_series, params.min_ridge_strength) {}
    //--------------------------------------------------------------------------//
    MarchingRidgesSeries3D(ScFieldSeriesBase_ptr_t p_sc_field_series_base,
                           real min_ridge_strength,
                           const VecI<3> &cell_resolution)
        : mp_sc_field_series_base(p_sc_field_series_base),
          m_mar_ridges(p_sc_field_series_base,
                       min_ridge_strength,
                       cell_resolution) {}
    //--------------------------------------------------------------------------//
    virtual ~MarchingRidgesSeries3D() {}
    //--------------------------------------------------------------------------//
    void searchAndSaveRidges(const std::string &save_dir)
    {
      // create directory if necessary
      std::filesystem::create_directories(save_dir);

      // prepare statistics
      utils::RuntimeAnalyzer analyzer("Execution");

      // go through all time steps
      uint num_steps = mp_sc_field_series_base->steppedTime().steps;
      utils::LoopPrinter printer("Search and save ridges", num_steps);
      for (uint time_step = 0; time_step < num_steps; ++time_step)
      {
        // load current step and calculate ridges
        mp_sc_field_series_base->loadStep(time_step);
        analyzer.startTimer("Execution");
        m_mar_ridges.calcRidgeSurface();
        analyzer.stopTimer("Execution");
        // save particles to file
        std::string file_path = save_dir + "/" + std::to_string(time_step) + ".data";
        m_mar_ridges.saveRidgeSurface(file_path);

        printer.endIteration();
      }
      analyzer.logAnalysis();
    }
    //--------------------------------------------------------------------------//
  private:
    //--------------------------------------------------------------------------//
    ScFieldSeriesBase_ptr_t mp_sc_field_series_base;
    MarchingRidges3D m_mar_ridges;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void transformMarRidgesScValues(ScalarFieldSeriesBase<3> &sc_series_base,
                                  const std::string &input_dir,
                                  const std::string &output_dir);
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//