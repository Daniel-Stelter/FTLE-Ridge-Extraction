#include "flowsim/ftlefield.hpp"
#include "flowsim/datasets/movingcircleridges.hpp"
#include "flowsim/datasets/doublegyre.hpp"
#include "images/visualization.hpp"
#include "experiments.hpp"
//--------------------------------------------------------------------------//
namespace dst
{
  //--------------------------------------------------------------------------//
  using namespace std;
  using namespace flowsim;
  using namespace images;
  //--------------------------------------------------------------------------//
  const string RESULTS_DIR = "./results/test-ridges/";
  //--------------------------------------------------------------------------//
  /**
   * Helper function for drawing the trajectories of particles to a texture.
   * Draws the whole trajectory as well as the discrete positions the particles
   * moved to. Additionally highlights the start and end positions.
   *
   * This function distinguishes between "accepted" and "rejected paths" which
   * is used to show whether a particle successfully found a ridge. `domain`
   * defines the visible domain pictured by the `image`. `paths` contains all
   * paths. A single path consists of a pair: The first component is a container
   * of 2D vectors which is an ordered list describing the trajectory. The
   * second component is a bool: true -> accepted and false -> rejected.
   * Ultimately, this results in different color encodings of the final
   * positions.
   */
  inline void drawPaths(Texture &texture,
                        const Domain<2> &domain,
                        const vector<pair<vector<VecR<2>>, bool>> &paths)
  {
    color col_trajectory = color(0.7);
    color col_discrete_path_points = color(0.5);
    color col_start_pos = color(0.3);
    color col_end_pos_accepted = color{0.0, 1.0, 0.0};
    color col_end_pos_rejected = color{1.0, 0.5, 0.0};

    // distinguish between accepted and rejected paths -> save respective final positions
    vector<VecR<2>> accepted, rejected;
    for (auto &path : paths)
    {
      if (path.second)
        accepted.push_back(path.first.back());
      else
        rejected.push_back(path.first.back());
    }
    // draw path as lines between positions of the path
    for (auto &path : paths)
      drawPath(texture, path.first, domain, col_trajectory);
    // draw discrete positions on top in a different color (to see step sizes)
    for (auto &path : paths)
      drawParticles(texture, path.first, domain, col_discrete_path_points);
    // draw first position with a larger point
    for (auto &path : paths)
      drawParticles(texture, {path.first[0]}, domain, col_start_pos, col_start_pos * 0.5);
    // draw rejected and accepted particles with a larger point in different colors
    drawParticles(texture, rejected, domain, col_end_pos_rejected, col_end_pos_rejected * 0.5);
    drawParticles(texture, accepted, domain, col_end_pos_accepted, col_end_pos_accepted * 0.5);
  }
  //--------------------------------------------------------------------------//
  /**
   * Helper function for drawing the trajectories of particles to a texture.
   * Draws the whole trajectory as well as the discrete positions the particles
   * moved to. Additionally highlights the start and end positions.
   *
   * `domain` defines the visible domain pictured by the `image`. `paths`
   * contains all paths. A single path consists of a container of 2D vectors
   * which is an ordered list describing the trajectory.
   */
  inline void drawPaths(Texture &texture,
                        const Domain<2> &domain,
                        const vector<vector<VecR<2>>> &paths)
  {
    color col_trajectory = color(0.7);
    color col_discrete_path_points = color(0.5);
    color col_start_pos = color(0.3);
    color col_end_pos = color(0.0, 0.7, 0.0);

    // draw path as lines between positions of the path
    for (auto &path : paths)
      drawPath(texture, path, domain, col_trajectory);
    // draw discrete positions on top in a different color (to see step sizes)
    for (auto &path : paths)
      drawParticles(texture, path, domain, col_discrete_path_points);
    // draw first position with a larger point
    for (auto &path : paths)
      drawParticles(texture, {path[0]}, domain, col_start_pos, col_start_pos * 0.5);
    // draw last position with a larger point
    for (auto &path : paths)
      drawParticles(texture, {path.back()}, domain, col_end_pos, col_end_pos * 0.5);
  }
  //--------------------------------------------------------------------------//
  /**
   * Test of the Stelter particle system on a small region of the DoubleGyre2D
   * dataset. Searches ridges for two time steps:
   * 9  -> initialize and constrain particles
   * 10 -> only constrain
   * The respective particle paths are drawn and saved.
   */
  void testStelterSeriesPathsDoubleGyre2D()
  {
    cout << "testStelterSeriesPathsDoubleGyre2D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "Stelter-DoubleGyre2D-paths";

    // SETUP: define flow and all further parameters
    shared_ptr<Flow<2>> p_flow = make_shared<DoubleGyre2D>();
    Domain<2> domain{Range(0.9, 1.15), Range(0.25, 0.5)};
    real t = 0.0;
    SteppedRange stepped_time(1.1, 11, false);
    VecI<2> img_res(1000);
    VecI<2> space_res = {100, 100};
    // set params
    params::ParamsStelter<2> params;
    params.par_init_res = params.voxel_res = space_res / 8;
    params.min_ridge_strength = 1.0;
    params.elliptic_dis_fac = 4.0;
    params.full_init_period = 1;
    // create series
    shared_ptr<InterpolScalarFieldSeries<2>> p_series =
        experiments::loadOrCreateItplFTLEFieldSeries<2>(p_flow,
                                                        t,
                                                        stepped_time,
                                                        space_res,
                                                        domain,
                                                        base_dir);
    auto p_series_base = dynamic_pointer_cast<ScalarFieldSeriesBase<2>>(p_series);
    experiments::createScFieldSeriesImages(p_series_base, img_res, base_dir);

    // INIT AND FIRST CONSTRAIN
    uint step = 9;
    ParticleSystemStelter<2> par_sys(p_series, params);
    p_series->loadStep(step);
    par_sys.initNewParticles(false);
    vector<pair<vector<VecR<2>>, bool>> paths;
    par_sys.constrainParticles(2, 5, &paths);
    // create image of the scalar field with particle trajectories
    Texture texture(base_dir + experiments::SubDirectories::IMGS_SC_FIELD + to_string(step) + ".ppm");
    drawPaths(texture, domain, paths);
    texture.writePPM(base_dir + "/paths/" + to_string(step) + ".ppm");

    // SECOND CONSTRAIN
    ++step;
    p_series->loadStep(step);
    par_sys.constrainParticles(5, 5, &paths);
    // create image of the scalar field with particle trajectories
    texture.readPPM(base_dir + experiments::SubDirectories::IMGS_SC_FIELD + to_string(step) + ".ppm");
    drawPaths(texture, domain, paths);
    texture.writePPM(base_dir + "/paths/" + to_string(step) + ".ppm");
  }
  //--------------------------------------------------------------------------//
  /**
   * Comparison of the constraining methods of both particle systems (Stelter
   * and Kindlmann) on the MovingCircleRidges2D example. The example uses a time
   * step for which both, very sharp and very shallow ridges, exist. Images
   * showing the trajectories of particles are generated for both methods.
   */
  void testParSysPathsCircles2D()
  {
    cout << "testParSysPathsCircles2D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "ParSys-Circles2D-paths/";

    // SETUP: define flow and all further parameters
    // Domain<2> domain{Range(0.4, 0.55), Range(0.4, 0.55)};
    Domain<2> domain{Range(0, 1), Range(0, 1)};
    shared_ptr<ScalarFieldBase<3>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges2D>(domain);
    real fixed_t = 4.0;
    VecI<2> img_res(1000);
    // set params
    VecI<2> par_init_res(25);
    VecI<2> voxel_res = par_init_res;
    real min_ridge_strength = 25.0;
    real sigma = 0.01;
    real elliptic_dis_fac = 4.0;
    VecR<2> num_diff_deltas(0.001);
    // prepare scalar field and texture
    shared_ptr<ScalarFieldBase<2>> p_reduced_sc_field =
        make_shared<ReducedScalarField<2, 1>>(p_time_dep_sc_field, VecI<1>{2}, VecR<1>{fixed_t});
    Texture texture_sc = createScalarFieldImage(*p_reduced_sc_field, img_res, Range(0, 1));
    texture_sc.writePPM(base_dir + "field.ppm");

    // execute Stelter par sys and save paths of constraining
    {
      vector<pair<vector<VecR<2>>, bool>> paths;
      ParticleSystemStelter<2> par_sys_stelter(p_reduced_sc_field,
                                               par_init_res,
                                               voxel_res,
                                               min_ridge_strength,
                                               sigma,
                                               elliptic_dis_fac,
                                               num_diff_deltas);
      par_sys_stelter.initNewParticles(false);
      par_sys_stelter.constrainParticles(2, 0, &paths);
      // create image of the scalar field with particle trajectories
      Texture texture_stelter(texture_sc);
      drawPaths(texture_stelter, domain, paths);
      texture_stelter.writePPM(base_dir + "stelter.ppm");
    }
    // execute Kindlmann par sys and save paths of constraining
    {
      vector<pair<vector<VecR<2>>, bool>> paths;
      ParticleSystemKindlmann<2> par_sys_kindlmann(p_reduced_sc_field,
                                                   par_init_res,
                                                   voxel_res,
                                                   min_ridge_strength,
                                                   sigma,
                                                   num_diff_deltas);
      par_sys_kindlmann.initNewParticles(false);
      par_sys_kindlmann.constrainParticles(2, 0, &paths);
      // create image of the scalar field with particle trajectories
      Texture texture_kindlmann(texture_sc);
      drawPaths(texture_kindlmann, domain, paths);
      texture_kindlmann.writePPM(base_dir + "kindlmann.ppm");
    }
  }
  //--------------------------------------------------------------------------//
  /**
   * Generates images which show information about the spectral decomposition
   * for the MovingCircleRidges2D example, including directions of eigenvectors
   * and magnitude of eigenvalues.
   */
  void testSpectralDecompositionCircles2D()
  {
    cout << "testSpectralDecompositionCircles2D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "SpecDecomp-Circles2D/";

    // SETUP: define flow and all further parameters
    shared_ptr<ScalarFieldBase<3>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges2D>();
    real fixed_t = 4.0;
    VecI<2> img_res(1000);
    real numeric_diff_delta = 0.001;
    // prepare scalar field
    shared_ptr<ScalarFieldBase<2>> p_reduced_sc_field =
        make_shared<ReducedScalarField<2, 1>>(p_time_dep_sc_field, VecI<1>{2}, VecR<1>{fixed_t});
    Domain<2> domain = p_reduced_sc_field->domain();

    // create scalar fields for spectral decomposition:
    // - min eigenvalue
    // - max eigenvalue
    // - product of eigenvalues
    // - direction of min eigenvector
    // - direction of max eigenvector
    NumericDiff<2> num_diff(p_reduced_sc_field, numeric_diff_delta);
    ScalarField<2> eig_val_min_sc([&](const VecR<2> &pos) -> real
                                  {
                                    auto H = num_diff.hessian(pos);
                                    Eigen::SelfAdjointEigenSolver<EMat<2, 2>> solver;
                                    solver.computeDirect(H);
                                    return solver.eigenvalues()[0]; },
                                  domain);
    ScalarField<2> eig_val_max_sc([&](const VecR<2> &pos) -> real
                                  {
                                    auto H = num_diff.hessian(pos);
                                    Eigen::SelfAdjointEigenSolver<EMat<2, 2>> solver;
                                    solver.computeDirect(H);
                                    return solver.eigenvalues()[1]; },
                                  domain);
    ScalarField<2> eig_val_mult_sc([&](const VecR<2> &pos) -> real
                                   {
                                    auto H = num_diff.hessian(pos);
                                    Eigen::SelfAdjointEigenSolver<EMat<2, 2>> solver;
                                    solver.computeDirect(H);
                                    return solver.eigenvalues()[0] * solver.eigenvalues()[1]; },
                                   domain);
    ScalarField<2> eig_dir_min_sc([&](const VecR<2> &pos) -> real
                                  {
                                    auto H = num_diff.hessian(pos);
                                    Eigen::SelfAdjointEigenSolver<EMat<2, 2>> solver;
                                    solver.computeDirect(H);
                                    EVec<2> eig_vec = solver.eigenvectors().col(0);
                                    return eig_vec[1] > 0.0
                                      ? acos(eig_vec.dot(EVec<2>(1.0, 0.0)))
                                      : acos(eig_vec.dot(EVec<2>(1.0, 0.0))) + M_PI; },
                                  domain);
    ScalarField<2> eig_dir_max_sc([&](const VecR<2> &pos) -> real
                                  {
                                    auto H = num_diff.hessian(pos);
                                    Eigen::SelfAdjointEigenSolver<EMat<2, 2>> solver;
                                    solver.computeDirect(H);
                                    EVec<2> eig_vec = solver.eigenvectors().col(1);
                                    return eig_vec[1] > 0.0
                                      ? acos(eig_vec.dot(EVec<2>(1.0, 0.0)))
                                      : acos(eig_vec.dot(EVec<2>(1.0, 0.0))) + M_PI; },
                                  domain);
    createScalarFieldImage(*p_reduced_sc_field, img_res, Range(0, 1)).writePPM(base_dir + "1-sc-field.ppm");
    createScalarFieldImage(eig_val_min_sc, img_res, Range(-10, 10), colormap::RED_BLUE).writePPM(base_dir + "2-eig-val-min.ppm");
    createScalarFieldImage(eig_val_max_sc, img_res, Range(-10, 10), colormap::RED_BLUE).writePPM(base_dir + "3-eig-val-max.ppm");
    createScalarFieldImage(eig_val_mult_sc, img_res, Range(-10, 10), colormap::RED_BLUE).writePPM(base_dir + "4-eig-vals-prod.ppm");
    createScalarFieldImage(eig_dir_min_sc, img_res, Range(0, 2 * M_PI), colormap::RING).writePPM(base_dir + "5-eig-dir-min.ppm");
    createScalarFieldImage(eig_dir_max_sc, img_res, Range(0, 2 * M_PI), colormap::RING).writePPM(base_dir + "6-eig-dir-max.ppm");
    colormap::createTextureRing(1000).writePPM(base_dir + "/ring.ppm");
  }
  //--------------------------------------------------------------------------//
  /**
   * Application of the Stelter particle system to the MovingCircleRidges2D dataset.
   */
  void testParSysSeriesStelterCircles2D()
  {
    cout << "testParSysSeriesStelterCircles2D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "Stelter-Circles2D";

    // SETUP: define flow and all further parameters
    shared_ptr<ScalarField<3>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges2D>();
    Domain<2> domain{p_time_dep_sc_field->domain(0), p_time_dep_sc_field->domain(1)};
    SteppedRange stepped_time(4.0, 100, false);
    VecI<2> img_res(1000);
    // set params
    VecI<2> par_init_res(200 / 8);
    VecI<2> voxel_res = par_init_res;
    real min_ridge_strength = 25.0;
    real sigma = 1.0 / 200;
    real elliptic_dis_fac = 4.0;
    uint full_init_period = 5;
    VecR<2> num_diff_deltas(sigma);
    // prepare scalar field and texture
    shared_ptr<ScalarFieldSeriesBase<2>> p_sc_series_base =
        make_shared<ScalarFieldSeries<2>>(p_time_dep_sc_field->function(), stepped_time, domain);
    experiments::createScFieldSeriesImages(p_sc_series_base, img_res, base_dir);

    // EXECUTE
    ParticleSystemSeriesStelter<2> par_sys_series(p_sc_series_base,
                                                  par_init_res,
                                                  voxel_res,
                                                  min_ridge_strength,
                                                  sigma,
                                                  elliptic_dis_fac,
                                                  full_init_period,
                                                  num_diff_deltas);
    par_sys_series.searchAndSaveRidges(base_dir + experiments::SubDirectories::DATA_PARTICLES);
    experiments::createScFieldSeriesRidgeImagesOnTexture(p_sc_series_base, base_dir);
  }
  //--------------------------------------------------------------------------//
  /**
   * Application of the Kindlmann particle system to the MovingCircleRidges2D dataset.
   */
  void testParSysSeriesKindlmannCircles2D()
  {
    cout << "testParSysSeriesKindlmannCircles2D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "Kindlmann-Circles2D";

    // SETUP: define flow and all further parameters
    shared_ptr<ScalarField<3>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges2D>();
    Domain<2> domain{p_time_dep_sc_field->domain(0), p_time_dep_sc_field->domain(1)};
    SteppedRange stepped_time(4.0, 100, false);
    VecI<2> img_res(1000);
    // set params
    params::ParamsKindlmann<2> params;
    params.par_init_res = params.voxel_res = VecI<2>(200 / 8);
    params.min_ridge_strength = 25.0;
    params.sigma = params.numeric_diff_delta = 1.0 / 200;
    // prepare scalar field and texture
    shared_ptr<ScalarFieldSeriesBase<2>> p_sc_series_base =
        make_shared<ScalarFieldSeries<2>>(p_time_dep_sc_field->function(), stepped_time, domain);
    experiments::createScFieldSeriesImages(p_sc_series_base, img_res, base_dir);

    // EXECUTE
    ParticleSystemSeriesKindlmann<2> par_sys_series(p_sc_series_base, params);
    par_sys_series.searchAndSaveRidges(base_dir + experiments::SubDirectories::DATA_PARTICLES);
    experiments::createScFieldSeriesRidgeImagesOnTexture(p_sc_series_base, base_dir);
  }
  //--------------------------------------------------------------------------//
  /**
   * Application of the Stelter particle system to the MovingCircleRidges3D dataset.
   */
  void testParSysSeriesStelterCircles3D()
  {
    cout << "testParSysSeriesStelterCircles3D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "Stelter-Circles3D";

    // SETUP: define flow and all further parameters
    shared_ptr<ScalarField<4>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges3D>();
    Domain<3> domain{p_time_dep_sc_field->domain(0), p_time_dep_sc_field->domain(1), p_time_dep_sc_field->domain(2)};
    SteppedRange stepped_time(4.0, 100, false);
    // set params
    VecI<3> par_init_res(200 / 8);
    VecI<3> voxel_res = par_init_res;
    real min_ridge_strength = 25.0;
    real sigma = 1.0 / 200;
    real elliptic_dis_fac = 4.0;
    uint full_init_period = 5;
    VecR<3> num_diff_deltas(sigma);
    // prepare scalar field
    shared_ptr<ScalarFieldSeriesBase<3>> p_sc_series =
        make_shared<ScalarFieldSeries<3>>(p_time_dep_sc_field->function(), stepped_time, domain);

    // EXECUTE
    ParticleSystemSeriesStelter<3> par_sys_series(p_sc_series,
                                                  par_init_res,
                                                  voxel_res,
                                                  min_ridge_strength,
                                                  sigma,
                                                  elliptic_dis_fac,
                                                  full_init_period,
                                                  num_diff_deltas);
    par_sys_series.searchAndSaveRidges(base_dir + experiments::SubDirectories::DATA_PARTICLES);
    // visualize distance to closest ridge
    {
      shared_ptr<ScalarField<4>> p_time_dep_sc_field_distance =
          make_shared<MovingCircleRidgesMinDistance3D>();
      shared_ptr<ScalarFieldSeriesBase<3>> p_sc_series_distance =
          make_shared<ScalarFieldSeries<3>>(p_time_dep_sc_field_distance->function(), stepped_time, domain);
      experiments::createParSysBlenderInput(p_sc_series, p_sc_series_distance, experiments::ParSysNormals::Stelter, base_dir);
      // move created files
      string blender_out = base_dir + experiments::SubDirectories::OUT_BLENDER;
      filesystem::create_directories(blender_out + "dist-error");
      for (filesystem::path p : filesystem::directory_iterator(blender_out))
        if (filesystem::is_regular_file(p))
          filesystem::rename(p, (blender_out + "dist-error") / p.filename());
    }
    // standard visualization with coloring of the scalar field
    experiments::createParSysBlenderInput(p_sc_series, experiments::ParSysNormals::Stelter, base_dir);
  }
  //--------------------------------------------------------------------------//
  /**
   * Application of the Kindlmann particle system to the MovingCircleRidges3D dataset.
   */
  void testParSysSeriesKindlmannCircles3D()
  {
    cout << "testParSysSeriesKindlmannCircles3D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "Kindlmann-Circles3D";

    // SETUP: define flow and all further parameters
    shared_ptr<ScalarField<4>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges3D>();
    Domain<3> domain{p_time_dep_sc_field->domain(0), p_time_dep_sc_field->domain(1), p_time_dep_sc_field->domain(2)};
    SteppedRange stepped_time(4.0, 100, false);

    params::ParamsKindlmann<3> params;
    params.par_init_res = params.voxel_res = VecI<3>(50);
    params.min_ridge_strength = 25.0;
    params.sigma = params.numeric_diff_delta = 1.0 / 200;

    // prepare scalar field
    shared_ptr<ScalarFieldSeriesBase<3>> p_sc_series =
        make_shared<ScalarFieldSeries<3>>(p_time_dep_sc_field->function(), stepped_time, domain);

    // execute particle system steps
    ParticleSystemSeriesKindlmann<3> par_sys_series(p_sc_series, params);
    par_sys_series.searchAndSaveRidges(base_dir + experiments::SubDirectories::DATA_PARTICLES);

    // visualizing distance to closest ridge
    {
      shared_ptr<ScalarField<4>> p_time_dep_sc_field_distance =
          make_shared<MovingCircleRidgesMinDistance3D>();
      shared_ptr<ScalarFieldSeriesBase<3>> p_sc_series_distance =
          make_shared<ScalarFieldSeries<3>>(p_time_dep_sc_field_distance->function(), stepped_time, domain);
      experiments::createParSysBlenderInput(p_sc_series, p_sc_series_distance, experiments::ParSysNormals::Kindlmann, base_dir);
      // move created files
      string blender_out = base_dir + experiments::SubDirectories::OUT_BLENDER;
      filesystem::create_directories(blender_out + "dist-error");
      for (filesystem::path p : filesystem::directory_iterator(blender_out))
        if (filesystem::is_regular_file(p))
          filesystem::rename(p, (blender_out + "dist-error") / p.filename());
    }

    // standard visualization with coloring of the scalar field
    experiments::createParSysBlenderInput(p_sc_series, experiments::ParSysNormals::Kindlmann, base_dir);
  }
  //--------------------------------------------------------------------------//
  /**
   * Application of Marching Ridges to the MovingCircleRidges3D dataset.
   */
  void testMarRidgesSeriesCircles3D()
  {
    cout << "testMarRidgesSeriesCircles3D" << endl;
    utils::printSeparator();

    // define output directory
    string base_dir = RESULTS_DIR + "MarRidges-Circles3D";

    // SETUP: define flow and all further parameters
    shared_ptr<ScalarField<4>> p_time_dep_sc_field =
        make_shared<MovingCircleRidges3D>();
    Domain<3> domain{p_time_dep_sc_field->domain(0), p_time_dep_sc_field->domain(1), p_time_dep_sc_field->domain(2)};
    SteppedRange stepped_time(4.0, 100, false);

    real min_ridge_strength = 25.0;
    VecI<3> cell_resolution(200);

    // prepare scalar field and texture
    shared_ptr<ScalarFieldSeriesBase<3>> p_sc_series_base =
        make_shared<ScalarFieldSeries<3>>(p_time_dep_sc_field->function(), stepped_time, domain);

    MarchingRidgesSeries3D mar_ridges_series(p_sc_series_base, min_ridge_strength, cell_resolution);
    mar_ridges_series.searchAndSaveRidges(base_dir + experiments::SubDirectories::OUT_BLENDER);

    // // visualizing distance to closest ridge
    shared_ptr<ScalarField<4>> p_time_dep_sc_field_distance =
        make_shared<MovingCircleRidgesMinDistance3D>();
    shared_ptr<ScalarFieldSeriesBase<3>> p_sc_series_distance =
        make_shared<ScalarFieldSeries<3>>(p_time_dep_sc_field_distance->function(), stepped_time, domain);
    transformMarRidgesScValues(*p_sc_series_distance,
                               base_dir + experiments::SubDirectories::OUT_BLENDER,
                               base_dir + experiments::SubDirectories::OUT_BLENDER + "dist-error/");
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//
int main()
{
  using namespace dst;
  utils::Timer::logCurrentTime(true);

  // testStelterSeriesPathsDoubleGyre2D();

  testParSysPathsCircles2D();

  // testSpectralDecompositionCircles2D();

  testParSysSeriesStelterCircles2D();
  testParSysSeriesKindlmannCircles2D();

  testParSysSeriesStelterCircles3D();
  testParSysSeriesKindlmannCircles3D();
  testMarRidgesSeriesCircles3D();
}
//--------------------------------------------------------------------------//