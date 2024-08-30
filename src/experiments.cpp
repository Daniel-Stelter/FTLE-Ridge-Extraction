#include "experiments.hpp"
//--------------------------------------------------------------------------//
#include <numeric>
//--------------------------------------------------------------------------//
namespace dst::experiments
{
  //--------------------------------------------------------------------------//
  using namespace std;
  using ScFieldSeriesBase2D_ptr_t = shared_ptr<flowsim::ScalarFieldSeriesBase<2>>;
  using ScFieldSeriesBase3D_ptr_t = shared_ptr<flowsim::ScalarFieldSeriesBase<3>>;
  using ItplScFieldSeries3D_ptr_t = shared_ptr<flowsim::InterpolScalarFieldSeries<3>>;
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  executeMarRidgesSeries(ItplScFieldSeries3D_ptr_t &p_itpl_sc_series,
                         const params::ParamsMarRidges<3> &params,
                         const string &base_dir)
  {
    string save_dir = base_dir + SubDirectories::OUT_BLENDER;
    flowsim::MarchingRidgesSeries3D mar_ridges_series(p_itpl_sc_series, params);
    mar_ridges_series.searchAndSaveRidges(save_dir);
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void createScFieldSeriesImages(ScFieldSeriesBase2D_ptr_t &p_sc_field_series,
                                 const VecI<2> &img_resolution,
                                 const string &base_dir)
  {
    createScFieldSeriesImages(p_sc_field_series, img_resolution, p_sc_field_series->domain(), base_dir);
  }
  //--------------------------------------------------------------------------//
  void createScFieldSeriesImages(ScFieldSeriesBase2D_ptr_t &p_sc_field_series,
                                 const VecI<2> &img_resolution,
                                 const Domain<2> &rendered_domain,
                                 const string &base_dir)
  {
    // subdirectory for saving
    string save_dir = base_dir + SubDirectories::IMGS_SC_FIELD;
    // go over all time steps
    uint num_steps = p_sc_field_series->steppedTime().steps;
    utils::LoopPrinter printer("Create scalar field images", num_steps);
    for (uint step = 0; step < num_steps; ++step)
    {
      string file_path = save_dir + to_string(step) + ".ppm";
      p_sc_field_series->loadStep(step);
      images::createScalarFieldImage(*p_sc_field_series, img_resolution, rendered_domain).writePPM(file_path);
      printer.endIteration();
    }
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void createScFieldSeriesRidgeImagesOnTexture(ScFieldSeriesBase2D_ptr_t &p_sc_field_series,
                                               const string &base_dir)
  {
    createScFieldSeriesRidgeImagesOnTexture(p_sc_field_series,
                                            p_sc_field_series->domain(),
                                            base_dir);
  }
  //--------------------------------------------------------------------------//
  void createScFieldSeriesRidgeImagesOnTexture(ScFieldSeriesBase2D_ptr_t &p_sc_field_series,
                                               const Domain<2> &rendered_domain,
                                               const string &base_dir)
  {
    using namespace std;
    // subdirectories for input and output
    string texture_dir = base_dir + SubDirectories::IMGS_SC_FIELD;
    string particle_dir = base_dir + SubDirectories::DATA_PARTICLES;
    string save_dir = base_dir + SubDirectories::IMGS_RIDGES_TEXTURE;

    uint num_steps = p_sc_field_series->steppedTime().steps;
    utils::LoopPrinter printer("Create ridge images (on texture)", num_steps);
    for (uint step = 0; step < num_steps; ++step)
    {
      string texture_path = texture_dir + to_string(step) + ".ppm";
      string particle_path = particle_dir + to_string(step) + ".data";
      images::Texture texture(texture_path);
      images::drawParticles(texture,
                            flowsim::VoxelContainer<2>::fetchFile(particle_path),
                            rendered_domain);
      texture.writePPM(save_dir + to_string(step) + ".ppm");
      printer.endIteration();
    }
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void createScFieldSeriesRidgeImagesWithValue(ScFieldSeriesBase2D_ptr_t &p_sc_field_series,
                                               const VecI<2> &img_resolution,
                                               const string &base_dir)
  {
    createScFieldSeriesRidgeImagesWithValue(p_sc_field_series,
                                            img_resolution,
                                            p_sc_field_series->domain(),
                                            base_dir);
  }
  //--------------------------------------------------------------------------//
  void createScFieldSeriesRidgeImagesWithValue(ScFieldSeriesBase2D_ptr_t &p_sc_field_series,
                                               const VecI<2> &img_resolution,
                                               const Domain<2> &rendered_domain,
                                               const string &base_dir)
  {
    using namespace std;
    // subdirectories for input and output
    string particle_dir = base_dir + SubDirectories::DATA_PARTICLES;
    string save_dir = base_dir + SubDirectories::IMGS_RIDGES_VALS;

    // prepare function for obtaining scalar values
    function<real(const VecR<2> &)> func = [&p_sc_field_series](const VecR<2> &pos)
    { return p_sc_field_series->value(pos); };

    // go over all time steps
    uint num_steps = p_sc_field_series->steppedTime().steps;
    utils::LoopPrinter printer("Create ridge images (with scalar value)", num_steps);
    for (uint step = 0; step < num_steps; ++step)
    {
      // try to find the saved particles
      string particle_path = particle_dir + to_string(step) + ".data";
      p_sc_field_series->loadStep(step);
      // create empty texture and load particles from file
      images::Texture texture(img_resolution[0], img_resolution[1]);
      vector<VecR<2>> particles = flowsim::VoxelContainer<2>::fetchFile(particle_path);
      // find function value for each particle
      vector<real> values(particles.size());
      transform(particles.begin(), particles.end(), values.begin(), func);
      // draw to texture and save as PPM file
      images::drawParticlesWithValue(texture,
                                     particles,
                                     values,
                                     rendered_domain);
      texture.writePPM(save_dir + to_string(step) + ".ppm");

      printer.endIteration();
    }
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void
  createParSysCloudCompareInput(ScFieldSeriesBase3D_ptr_t &p_sc_field_series,
                                ScFieldSeriesBase3D_ptr_t &p_sc_field_series_for_values,
                                ParSysNormals par_sys_normals,
                                const string &base_dir)
  {
    /** File structure (see https://www.cloudcompare.org/doc/wiki/index.php/BIN):
     *  1. Number of clouds (uint, 32 bit)
     *  2. PER CLOUD
     *    2.1 Number of points in cloud (uint, 32 bit)
     *    2.2 Flags / options (uchar, 8 bit)
     *    2.3 Opt: Name (uchar[...] until char 0)
     *  3. PER POINT
     *    3.1 Position (3 * float, 3 * 32 bit)
     *    3.2 Opt: Color (3 * uchar, 3 * 8 bit)
     *    3.3 Opt: Normal (3 * float, 3 * 32 bit)
     *    3.4 Opt: Scalar Field (double, 64 bit)
     */
    using namespace std;
    using VecCC = VC::vecn::vecn<3, float>;
    struct CloudCompareInfo
    {
      VecCC pos;
      VecCC normal;
      double value;
    };
    constexpr uint BUFFER_SIZE = 50000;

    string particle_dir = base_dir + SubDirectories::DATA_PARTICLES;
    string save_dir = base_dir + SubDirectories::OUT_CLOUDCOMPARE;
    filesystem::create_directories(save_dir);

    VecR<3> deltas(0.001);
    auto p_itpl_sc_series = dynamic_pointer_cast<flowsim::InterpolScalarFieldSeries<3>>(p_sc_field_series);
    if (p_itpl_sc_series)
      deltas = p_itpl_sc_series->cellSideLengths();
    flowsim::SmoothNumericDiff<3> num_diff(p_sc_field_series, deltas);

    uint num_steps = p_sc_field_series->steppedTime().steps;
    utils::LoopPrinter printer("Create CloudCompare input files", num_steps);
    for (uint step = 0; step < num_steps; ++step)
    {
      // load step of series
      p_sc_field_series->loadStep(step);
      if (p_sc_field_series.get() != p_sc_field_series_for_values.get())
        p_sc_field_series_for_values->loadStep(step);
      // read created particles, and write with extra information to outfile
      string particle_path = particle_dir + to_string(step) + ".data";
      string save_path = save_dir + to_string(step) + ".bin";
      ifstream infile(particle_path, ios::binary);
      ofstream outfile(save_path, ios::binary);
      if (!infile)
        throw runtime_error("createParSysCloudCompareInput: could not read particles from file '" + particle_path + "'");
      if (!outfile)
        throw runtime_error("createParSysCloudCompareInput: could not write to file '" + save_path + "'");

      // first entry is number of extracted particles
      size_t num_particles;
      infile.read((char *)&num_particles, sizeof(size_t));

      // 1. Number of clouds (uint, 32 bit)
      uint32_t temp_uint = 1;
      outfile.write((char *)&temp_uint, sizeof(temp_uint));
      // 2.1 Number of points in cloud (uint, 32 bit)
      temp_uint = num_particles;
      outfile.write((char *)&temp_uint, sizeof(temp_uint));
      // 2.2 Flags (uchar, 8 bit): always ON | normals | scalar field
      u_char temp_uchar = 1 << 0 | 1 << 2 | 1 << 3;
      outfile.write((char *)&temp_uchar, sizeof(temp_uchar));

      vector<VecR<3>> input(BUFFER_SIZE);
      vector<CloudCompareInfo> output(BUFFER_SIZE);
      // iterate in multiple steps over input file
      for (uint buffer_iter = 0; buffer_iter * BUFFER_SIZE < num_particles; ++buffer_iter)
      {
        // use start id to find out how much the buffer will be filled in this iteration
        uint buffer_start = buffer_iter * BUFFER_SIZE;
        uint buffer_fill = buffer_start + BUFFER_SIZE < num_particles ? BUFFER_SIZE : num_particles - buffer_start;
        infile.read((char *)&input[0], buffer_fill * sizeof(VecR<3>));
        // calculate extended info in parallel
#pragma omp parallel for schedule(dynamic)
        for (uint buffer_id = 0; buffer_id < buffer_fill; ++buffer_id)
        {
          VecR<3> pos = input[buffer_id];
          VecR<3> normal = ParSysNormals::Stelter == par_sys_normals
                               ? flowsim::ridgeinfo::stelter::direction<3>(num_diff, pos)
                               : flowsim::ridgeinfo::kindlmann::direction<3>(num_diff, pos);
          real value = p_sc_field_series_for_values->value(pos);
          output[buffer_id] = {VecCC(pos[0], pos[1], pos[2]),          // 3.1 Position (3 * float, 3 * 32 bit)
                               VecCC(normal[0], normal[1], normal[2]), // 3.3 Opt: Normal (3 * float, 3 * 32 bit)
                               (double)value};                         // 3.4 Opt: Scalar Field (double, 64 bit)
        }
        // actually write after all entries in buffer were used
        outfile.write((char *)&output[0], buffer_fill * sizeof(CloudCompareInfo));
      }
      infile.close();
      outfile.close();

      printer.endIteration();
    }
  }
  //--------------------------------------------------------------------------//
  void
  createParSysCloudCompareInput(ScFieldSeriesBase3D_ptr_t &p_sc_field_series,
                                ParSysNormals par_sys_normals,
                                const string &base_dir)
  {
    createParSysCloudCompareInput(p_sc_field_series, p_sc_field_series, par_sys_normals, base_dir);
  }
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  void createParSysBlenderInput(ScFieldSeriesBase3D_ptr_t &p_sc_field_series,
                                ScFieldSeriesBase3D_ptr_t &p_sc_field_series_for_values,
                                ParSysNormals par_sys_normals,
                                const string &base_dir)
  {
    using namespace std;

    using NumParsBlender_t = u_int32_t;
    using FloatBlender_t = float;
    using VecBlender_t = VC::vecn::vecn<3, FloatBlender_t>;

    string particle_dir = base_dir + SubDirectories::DATA_PARTICLES;
    string save_dir = base_dir + SubDirectories::OUT_BLENDER;
    filesystem::create_directories(save_dir);

    VecR<3> deltas(0.001);
    auto p_itpl_sc_series = dynamic_pointer_cast<flowsim::InterpolScalarFieldSeries<3>>(p_sc_field_series);
    if (p_itpl_sc_series)
      deltas = p_itpl_sc_series->cellSideLengths();
    flowsim::SmoothNumericDiff<3> num_diff(p_sc_field_series, deltas);

    real global_min = Globals::INF;
    real global_max = -Globals::INF;
    real global_sum = 0.0;
    uint global_num_particles = 0;

    uint num_steps = p_sc_field_series->steppedTime().steps;
    utils::LoopPrinter printer("Create Blender input files", num_steps);
    for (uint step = 0; step < num_steps; ++step)
    {
      // load step of series
      p_sc_field_series->loadStep(step);
      if (p_sc_field_series.get() != p_sc_field_series_for_values.get())
        p_sc_field_series_for_values->loadStep(step);
      // read created particles, and write with extra information to outfile
      string particle_path = particle_dir + to_string(step) + ".data";
      string save_path = save_dir + to_string(step) + ".data";
      ofstream outfile(save_path, ios::binary);
      if (!outfile)
        throw runtime_error("createParSysCloudCompareInput: could not write to file '" + save_path + "'");

      // load all particles and write population size
      vector<VecR<3>> particles = flowsim::VoxelContainer<3>::fetchFile(particle_path);
      NumParsBlender_t num_particles = particles.size();
      outfile.write((char *)&num_particles, sizeof(NumParsBlender_t));

      // pre-calculate all scalar values
      vector<FloatBlender_t> vals(num_particles);
#pragma omp parallel for schedule(dynamic)
      for (uint par_idx = 0; par_idx < num_particles; ++par_idx)
        vals[par_idx] = p_sc_field_series_for_values->value(particles[par_idx]);

      // find and write min and max scalar values
      auto iters_minmax = minmax_element(vals.begin(), vals.end());
      FloatBlender_t val_min = iters_minmax.first == vals.end() ? -Globals::INF : *iters_minmax.first;
      FloatBlender_t val_max = iters_minmax.second == vals.end() ? Globals::INF : *iters_minmax.second;
      outfile.write((char *)&val_min, sizeof(FloatBlender_t));
      outfile.write((char *)&val_max, sizeof(FloatBlender_t));

      // print out all particle positions with corresponding normals and scalar values
#pragma omp parallel for schedule(dynamic)
      for (uint par_idx = 0; par_idx < num_particles; ++par_idx)
      {
        const VecR<3> &real_pos = particles[par_idx];
        // normal computation may take some time, thus parallel execution
        VecR<3> real_normal = ParSysNormals::Stelter == par_sys_normals
                                  ? flowsim::ridgeinfo::stelter::direction<3>(num_diff, real_pos)
                                  : flowsim::ridgeinfo::kindlmann::direction<3>(num_diff, real_pos);
        // convert data
        VecBlender_t pos(real_pos[0], real_pos[1], real_pos[2]);
        VecBlender_t normal(real_normal[0], real_normal[1], real_normal[2]);
#pragma omp critical
        {
          outfile.write((char *)&pos, sizeof(VecBlender_t));
          outfile.write((char *)&normal, sizeof(VecBlender_t));
          outfile.write((char *)&vals[par_idx], sizeof(FloatBlender_t));
        }
      }
      outfile.close();

      // handle statistics
      if (val_min < global_min && val_min != -Globals::INF)
        global_min = val_min;
      if (val_max > global_max && val_max != Globals::INF)
        global_max = val_max;
      global_sum += std::accumulate(vals.begin(), vals.end(), 0.0);
      global_num_particles += num_particles;

      printer.endIteration();
    }
    string stats_save_path = save_dir + "stats.txt";
    ofstream out_stats(stats_save_path);
    out_stats << "Min:  " << global_min
              << "\nMax:  " << global_max
              << "\nMean: " << (global_sum / global_num_particles);
    out_stats.close();
  }
  //--------------------------------------------------------------------------//
  void createParSysBlenderInput(ScFieldSeriesBase3D_ptr_t &p_sc_field_series,
                                ParSysNormals par_sys_normals,
                                const string &base_dir)
  {
    createParSysBlenderInput(p_sc_field_series, p_sc_field_series, par_sys_normals, base_dir);
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//