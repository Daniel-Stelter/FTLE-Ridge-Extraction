#pragma once
//--------------------------------------------------------------------------//
#include "ridgeinfo.hpp"
#include "numericdiff.hpp"
#include "scalarfield.hpp"
#include "../params.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  real energyProfile(real r);
  //--------------------------------------------------------------------------//
  real energyProfileDerivative(real r);
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  struct EntryState
  {
    using State_t = unsigned char;
    static constexpr State_t NONE = 0;
    static constexpr State_t LOCKED = 1 << 0;  // entry is currently locked / inaccessible
    static constexpr State_t DELETE = 1 << 1;  // entry will be deleted by calling applyDelete()
    static constexpr State_t UPDATED = 1 << 2; // entry was updated lately
  };
  //--------------------------------------------------------------------------//
  template <uint n>
  class VoxelContainer
  {
    //--------------------------------------------------------------------------//
    using VecList_t = std::vector<VecR<n>>;
    using StateList_t = std::vector<EntryState::State_t>;
    using IndexList_t = std::vector<uint>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    VoxelContainer() { omp_init_nest_lock(&m_lck); }
    //--------------------------------------------------------------------------//
    VoxelContainer(const VecI<n> &resolution, const Domain<n> &domain)
        : m_voxel_grid(resolution, domain)
    {
      omp_init_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    VoxelContainer(const Grid<n> &grid)
        : VoxelContainer(grid.cellResolution(), grid.domain()) {}
    //--------------------------------------------------------------------------//
    ~VoxelContainer() { omp_destroy_nest_lock(&m_lck); }
    //--------------------------------------------------------------------------//
    const VecI<n> &resolution() const { return m_voxel_grid.cellResolution(); }
    //--------------------------------------------------------------------------//
    const Domain<n> &domain() const { return m_voxel_grid.domain(); }
    //--------------------------------------------------------------------------//
    const Grid<n> &baseGrid() const { return m_voxel_grid; }
    //--------------------------------------------------------------------------//
    size_t size() const { return m_entry_list.size(); }
    //--------------------------------------------------------------------------//
    const VecList_t &entries() const { return m_entry_list; }
    VecR<n> entry(uint eid) const { return m_entry_list[eid]; }
    //--------------------------------------------------------------------------//
    // CAUTION! This is not thread-safe and must be used in isolation -> locking
    const IndexList_t &voxelCell(const VecR<n> &pos) const
    {
      return m_voxel_grid.cellData(m_voxel_grid.posToCellCoord(pos));
    }
    //--------------------------------------------------------------------------//
    // CAUTION! This is not thread-safe and must be used in isolation -> locking
    IndexList_t &voxelCell(const VecR<n> &pos)
    {
      return m_voxel_grid.cellData(m_voxel_grid.posToCellCoord(pos));
    }
    //--------------------------------------------------------------------------//
    void reset()
    {
      omp_set_nest_lock(&m_lck);
      m_entry_list.clear();
      m_entry_state_list.clear();
      for (auto &cell : m_voxel_grid.cellData())
        cell.clear();
      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    void insert(const VecR<n> &entry)
    {
      omp_set_nest_lock(&m_lck);
      uint list_id = m_entry_list.size();

      m_entry_list.push_back(entry);
      // no flags set after creation
      m_entry_state_list.push_back(EntryState::NONE);

      // add reference as index into voxel grid
      voxelCell(entry).push_back(list_id);
      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    void update(uint eid, const VecR<n> &upd_entry, bool set_updated_state = true)
    {
      omp_set_nest_lock(&m_lck);
      assert(eid < m_entry_list.size());
      // get old an new voxel cells
      VecI<n> old_cell_coord = m_voxel_grid.posToCellCoord(m_entry_list[eid]);
      VecI<n> new_cell_coord = m_voxel_grid.posToCellCoord(upd_entry);
      // if cells differ, delete id from old and add to new voxel
      if (old_cell_coord != new_cell_coord)
      {
        auto &r_old_cell = m_voxel_grid.cellData(old_cell_coord);
        auto &r_new_cell = m_voxel_grid.cellData(new_cell_coord);
        r_old_cell.erase(std::find(r_old_cell.begin(), r_old_cell.end(), eid));
        r_new_cell.push_back(eid);
      }
      // in any case we replace the old position with the new one
      m_entry_list[eid] = upd_entry;
      if (set_updated_state)
        setState(eid, EntryState::UPDATED);
      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    void applyDelete()
    {
      omp_set_nest_lock(&m_lck);
      // map the old indices to the new ones
      IndexList_t id_mapping(size());
      uint current_id = 0;
      VecList_t new_entries;
      StateList_t new_entry_states;
      // sort out entries and find new indices
      for (uint eid = 0; eid < size(); ++eid)
      {
        if (isDeleted(eid))
          id_mapping[eid] = size(); // deletion id is encoded as out-of-range id 'size()'
        else
        {
          id_mapping[eid] = current_id++;
          new_entries.push_back(m_entry_list[eid]);
          new_entry_states.push_back(m_entry_state_list[eid]);
        }
      }
      // update indices
      for (IndexList_t &cell : m_voxel_grid.cellData())
      {
        // go over all entries inside this cell in reverse order
        for (int x = cell.size() - 1; x >= 0; --x)
        {
          uint old_id = cell[x];
          uint new_id = id_mapping[old_id];
          // delete step: overwrite with last element, then pop back (CHANGES ORDER OF IDS!)
          if (new_id == size())
          {
            cell[x] = cell.back();
            cell.pop_back();
          }
          else
            cell[x] = new_id;
        }
      }
      m_entry_list = std::move(new_entries);
      m_entry_state_list = std::move(new_entry_states);
      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    void unsetUpdatedStates()
    {
      omp_set_nest_lock(&m_lck);
#pragma omp parallel for schedule(dynamic)
      for (uint eid = 0; eid < m_entry_state_list.size(); ++eid)
        unsetState(eid, EntryState::UPDATED);
      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    IndexList_t findSupportNeighborIndices(const VecR<n> &entry,
                                           const real radius,
                                           EntryState::State_t required = EntryState::NONE,
                                           EntryState::State_t ignored = EntryState::NONE)
    {
      omp_set_nest_lock(&m_lck);
      IndexList_t result;
      // go through all neighboring cells (own cell included)
      for (VecI<n> cell_coord : m_voxel_grid.cellNeighborCoords(entry, radius))
      {
        for (uint nbh_eid : m_voxel_grid.cellData(cell_coord))
        {
          assert(nbh_eid < size());
          // skip on exact match
          if (same(entry, m_entry_list[nbh_eid]))
            continue;
          // apply filter
          if (!hasAllStates(nbh_eid, required) || hasAtLeastOneState(nbh_eid, ignored))
            continue;
          VecR<n> dif = m_entry_list[nbh_eid] - entry;
          // begin with very cheap norm1 check and only if needed check for norm2
          if (norm1(dif) < radius && norm2(dif) < radius)
            result.push_back(nbh_eid);
        }
      }
      omp_unset_nest_lock(&m_lck);
      return result;
    }
    //--------------------------------------------------------------------------//
    void setIsolation(uint eid, real radius)
    {
      assert(eid < size());
      bool all_neighbors_free;
      do
      {
        omp_set_nest_lock(&m_lck);
        // neighbors may change -> recalculation needed
        IndexList_t neighbor_ids = findSupportNeighborIndices(m_entry_list[eid], radius);

        all_neighbors_free = true;
        for (uint nid : neighbor_ids)
        {
          // find out if there is an isolated neighbor
          if (hasAllStates(nid, EntryState::LOCKED))
          {
            all_neighbors_free = false;
            break;
          }
        }
        // if all neighbors are ready, mark this entry as locked
        if (all_neighbors_free)
          setState(eid, EntryState::LOCKED);

        omp_unset_nest_lock(&m_lck);
      } while (!all_neighbors_free);
    }
    //--------------------------------------------------------------------------//
    void unsetIsolation(uint eid) { unsetState(eid, EntryState::LOCKED); }
    //--------------------------------------------------------------------------//
    void setDelete(uint eid) { setState(eid, EntryState::DELETE); }
    //--------------------------------------------------------------------------//
    bool isDeleted(uint eid) { return hasAllStates(eid, EntryState::DELETE); }
    //--------------------------------------------------------------------------//
    void writeFile(const std::string &file_path)
    {
      omp_set_nest_lock(&m_lck);
      std::filesystem::create_directories(std::filesystem::absolute(file_path).parent_path());
      std::ofstream out(file_path, std::ios::binary);
      if (!out)
        throw std::runtime_error("VoxelContainer: could not write to file '" + file_path + "'");

      size_t size = m_entry_list.size();
      out.write((char *)&size, sizeof(size_t));
      out.write((char *)&m_entry_list[0], size * sizeof(VecR<n>));
      out.close();

      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    void readFile(const std::string &file_path)
    {
      omp_set_nest_lock(&m_lck);
      // delete contents of current stuff
      for (auto &cell : m_voxel_grid.cellData())
        cell.clear();
      m_entry_list.clear();
      m_entry_state_list.clear();

      auto data = fetchFile(file_path);
      for (const VecR<n> &entry : data)
        insert(entry);
      omp_unset_nest_lock(&m_lck);
    }
    //--------------------------------------------------------------------------//
    static VecList_t fetchFile(const std::string &file_path)
    {
      VecList_t data;
      std::ifstream in(file_path, std::ios::binary);
      if (!in)
        throw std::runtime_error("VoxelContainer: could not read file '" + file_path + "'");

      size_t size;
      in.read((char *)&size, sizeof(size_t));
      data.resize(size);
      in.read((char *)&data[0], size * sizeof(VecR<n>));
      in.close();
      return data;
    }
    //--------------------------------------------------------------------------//
  private:
    //--------------------------------------------------------------------------//
    bool hasAllStates(uint eid, EntryState::State_t state)
    {
      assert(eid < size());
      return (m_entry_state_list[eid] & state) == state;
    }
    bool hasAtLeastOneState(uint eid, EntryState::State_t state)
    {
      assert(eid < size());
      return (m_entry_state_list[eid] & state) != EntryState::NONE;
    }
    void setState(uint eid, EntryState::State_t state)
    {
      assert(eid < size());
      m_entry_state_list[eid] |= state;
    }
    void unsetState(uint eid, EntryState::State_t state)
    {
      assert(eid < size());
      m_entry_state_list[eid] &= ~state;
    }
    //--------------------------------------------------------------------------//
    omp_nest_lock_t m_lck;
    VecList_t m_entry_list;
    StateList_t m_entry_state_list;
    CellDataGrid<n, IndexList_t> m_voxel_grid;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleBase
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ParticleBase(VoxelContainer<n> &particle_container,
                 const NumericDiff<n> &num_diff,
                 const Grid<n> &par_init_grid,
                 real sigma,
                 real min_ridge_strength,
                 real elliptic_factor)
        : m_particle_container(particle_container),
          m_num_diff(num_diff),
          m_par_init_grid(par_init_grid),
          m_sigma(sigma),
          m_inv_sigma(1.0 / sigma),
          m_min_ridge_strength(min_ridge_strength),
          m_elliptic_factor(elliptic_factor),
          m_sq_elliptic_factor(elliptic_factor * elliptic_factor) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleBase() {}
    //--------------------------------------------------------------------------//
    virtual void init(const VecR<n> &pos) = 0;
    //--------------------------------------------------------------------------//
    virtual bool constrain(uint max_voxel_distance = 2, std::vector<VecR<n>> *p_path = nullptr) = 0;
    //--------------------------------------------------------------------------//
    const VecR<n> &pos() const { return m_pos; }
    //--------------------------------------------------------------------------//
    real ridgeStrength() { return m_ridge_strength; }
    //--------------------------------------------------------------------------//
    bool hasMinRidgeStrength() { return m_ridge_strength >= m_min_ridge_strength; }
    //--------------------------------------------------------------------------//
    const EMat<n, n> &ridgeTangent() { return m_ridge_tangent; }
    //--------------------------------------------------------------------------//
    const VecR<n> &ridgeDir() { return m_ridge_dir; }
    //--------------------------------------------------------------------------//
    void advect(const VecR<n> &dir)
    {
      const Domain<n> &domain = this->m_num_diff.domain();
      assert(domain.isInside(m_pos));
      // assure that new position does not go out of scope
      VecR<n> new_pos(m_pos + dir);
      for (uint i = 0; i < n; ++i)
      {
        if (new_pos[i] < domain[i].min)
        {
          const real fac = (domain[i].min - new_pos[i]) / (m_pos[i] - new_pos[i]);
          new_pos = fac * m_pos + (1 - fac) * new_pos;
        }
        else if (new_pos[i] > domain[i].max)
        {
          const real fac = (new_pos[i] - domain[i].max) / (new_pos[i] - m_pos[i]);
          new_pos = fac * m_pos + (1 - fac) * new_pos;
        }
      }
      init(new_pos);
    }
    //--------------------------------------------------------------------------//
    /// Modified distance to another position using an elliptical distance which is halved
    /// in ridge tangent direction(s). Distance in tangent normal direction stays the same.
    real ellipticDistance(const VecR<n> &other) const
    {
      const real ss = m_sq_elliptic_factor;
      const VecR<n> r = m_pos - other;
      const real norm_r = VC::vecn::norm(r);
      const real cos_pow_2 = pow(m_ridge_dir | r, 2) / pow(norm_r, 2);
      return norm_r * sqrt(cos_pow_2 + (1 - cos_pow_2) / ss);
    }
    //--------------------------------------------------------------------------//
    VecR<n> ellipticDistanceDeriv(const VecR<n> &other) const
    {
      const real s = m_elliptic_factor;
      const real ssm1 = m_sq_elliptic_factor - 1; // squared s minus 1
      const VecR<n> r = m_pos - other;
      const real norm_r = VC::vecn::norm(r);
      const real dot = m_ridge_dir | r;
      return (m_ridge_dir * (ssm1 * dot) + r) / (s * sqrt(pow(norm_r, 2) + ssm1 * pow(dot, 2)));
    }
    //--------------------------------------------------------------------------//
    std::vector<uint> neighborIds(real range_mult = 1.0,
                                  EntryState::State_t required = EntryState::NONE,
                                  EntryState::State_t ignored = EntryState::NONE)
    {
      const real s = m_elliptic_factor;
      const real max_dis = range_mult * m_sigma;
      // use sigma for neighborhood
      std::vector<uint> neighborhood =
          m_particle_container.findSupportNeighborIndices(m_pos, s * max_dis, required, ignored);
      // check all neighbors if they lie inside elliptic neighborhood
      for (int id = neighborhood.size() - 1; id >= 0; --id)
      {
        VecR<n> pos = m_particle_container.entry(neighborhood[id]);
        if (ellipticDistance(pos) > max_dis)
          neighborhood.erase(neighborhood.begin() + id);
      }
      return neighborhood;
    }
    //--------------------------------------------------------------------------//
    std::vector<VecR<n>> neighbors(real range_mult = 1.0,
                                   EntryState::State_t required = EntryState::NONE,
                                   EntryState::State_t ignored = EntryState::NONE)
    {
      std::vector<uint> nbh_ids = neighborIds(range_mult, required, ignored);
      std::vector<VecR<n>> neighborhood;
      neighborhood.reserve(nbh_ids.size());
      for (uint nid : nbh_ids)
        neighborhood.push_back(m_particle_container.entry(nid));
      return neighborhood;
    }
    //--------------------------------------------------------------------------//
    real localEnergyInfluence(const std::vector<VecR<n>> &neighbors)
    {
      real energy = 0.0;
      for (const VecR<n> &other : neighbors)
        energy += energyProfile(ellipticDistance(other) * m_inv_sigma);
      return energy;
    }
    //--------------------------------------------------------------------------//
    real localEnergyInfluence(real range_radius = 1.0,
                              EntryState::State_t required = EntryState::NONE,
                              EntryState::State_t ignored = EntryState::DELETE)
    {
      return localEnergyInfluence(neighbors(range_radius, required, ignored));
    }
    //--------------------------------------------------------------------------//
    uint traveledVoxelDistance(const VecR<n> &start_pos) const
    {
      return m_par_init_grid.cellCoordDistance(start_pos, m_pos);
    }
    //--------------------------------------------------------------------------//
    void updateByEnergy(uint iters = 5,
                        EntryState::State_t required = EntryState::NONE,
                        EntryState::State_t ignored = EntryState::DELETE)
    {
      for (uint iter = 0; iter < iters; ++iter)
      {
        // find step direction using all particles in the neighborhood and sum them up
        VecR<n> energy_dir(0.0);
        for (const VecR<n> other : neighbors(1.0, required, ignored))
        {
          const real rij = ellipticDistance(other);
          const VecR<n> d_rij = ellipticDistanceDeriv(other);
          const real d_Eij = energyProfileDerivative(rij * m_inv_sigma);
          energy_dir -= d_rij * (d_Eij * rij / pow(VC::vecn::norm(d_rij), 2));
        }
        // return if (nearly) no movement
        if (norminf(energy_dir) <= Globals::SMALL)
          return;

        // scale energy direction for reasonable movement
        energy_dir *= 2.0 * m_sigma;
        // project it onto the non-full rank tangent matrix
        VecR<n> proj_dir(EVec<n>(m_ridge_tangent * EVec<n>(energy_dir.data.data())).data());
        // VecR<n> start_pos(m_pos);
        advect(proj_dir);
        // // if not able to find better position: reset pos and break up
        // if (!constrain(false))
        // {
        //   init(start_pos);
        //   return;
        // }
      }
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    // Particle settings
    VoxelContainer<n> &m_particle_container;
    const NumericDiff<n> &m_num_diff;
    const Grid<n> &m_par_init_grid;
    const real m_sigma;
    const real m_inv_sigma;
    const real m_min_ridge_strength;
    const real m_elliptic_factor;
    const real m_sq_elliptic_factor;
    // Actual local particle infomation
    VecR<n> m_pos;
    EMat<n, n> m_ridge_tangent;
    real m_ridge_strength;
    VecR<n> m_ridge_dir; // direction to next ridge as a unit vector
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleStelter : public ParticleBase<n>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ParticleStelter(VoxelContainer<n> &particle_container,
                    const NumericDiff<n> &num_diff,
                    const Grid<n> &par_init_grid,
                    real sigma,
                    real min_ridge_strength,
                    real elliptic_factor)
        : ParticleBase<n>(particle_container, num_diff, par_init_grid, sigma, min_ridge_strength, elliptic_factor) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleStelter() {}
    //--------------------------------------------------------------------------//
    virtual void init(const VecR<n> &pos) override
    {
      if (VC::vecn::norminf(this->m_pos - pos) == 0)
        return;
      // set particle onto domain border if outside
      const Domain<n> &domain = this->m_num_diff.domain();
      for (uint i = 0; i < n; ++i)
        this->m_pos[i] = std::min(std::max(pos[i], domain[i].min), domain[i].max);
      // initialize local ridge information
      ridgeinfo::stelter::all<n>(this->m_num_diff,
                                 this->m_pos,
                                 &this->m_ridge_strength,
                                 &this->m_ridge_tangent,
                                 &this->m_ridge_dir);
      real norm_dir = VC::vecn::norm(this->m_ridge_dir);
      if (norm_dir > Globals::SMALL)
        this->m_ridge_dir /= norm_dir;
    }
    //--------------------------------------------------------------------------//
    virtual bool constrain(uint max_voxel_distance = 2, std::vector<VecR<n>> *p_path = nullptr) override
    {
      // use references to members for current position and direction
      VecR<n> last_pos, last_dir;
      const VecR<n> &cur_pos(this->m_pos), &cur_dir(this->m_ridge_dir);
      VecR<n> start_pos(cur_pos);

      // for debug / visualization purpose: save path positions
      auto appendCurPos = [&]()
      { if (p_path) p_path->push_back(cur_pos); };
      appendCurPos(); // insert first position

      // explicit Euler with abs step length until direction flips
      int step = 0;
      do // until angle between last and new dir near 180°
      {
        // check whether too many steps were executed
        if (step >= Globals::MAX_CONSTRAIN_STEPS)
          return false;
        ++step;
        // update position and all local information - step length is sigma
        last_pos = cur_pos;
        last_dir = cur_dir;
        this->advect(cur_dir * this->m_sigma);
        appendCurPos();

        // if voxel travel distance is limited and exceeded, return
        if (this->traveledVoxelDistance(start_pos) > max_voxel_distance)
          return false;
      } while ((cur_dir | last_dir) > -0.98); // indicates flip

      VecR<n> positions[2] = {cur_pos, last_pos};
      VecR<n> grad, middle;
      real min_dot_abs_val = Globals::INF;
      VecR<n> min_dot_pos(real(0));
      for (uint x = 0; x < 5; ++x) // bisection
      {
        middle = (positions[0] + positions[1]) * 0.5;
        grad = this->m_num_diff.grad(middle);
        real dot = cur_dir | grad;
        // update best position
        if (abs(dot) < min_dot_abs_val)
        {
          min_dot_abs_val = abs(dot);
          min_dot_pos = middle;
        }
        // update segment
        if (dot >= 0.0)
          positions[0] = middle;
        else
          positions[1] = middle;
      }
      init(min_dot_pos);

      appendCurPos();

      return this->hasMinRidgeStrength();
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleKindlmann : public ParticleBase<n>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ParticleKindlmann(VoxelContainer<n> &particle_container,
                      const NumericDiff<n> &num_diff,
                      const Grid<n> &par_init_grid,
                      real sigma,
                      real min_ridge_strength,
                      real elliptic_factor = 1)
        : ParticleBase<n>(particle_container, num_diff, par_init_grid, sigma, min_ridge_strength, elliptic_factor) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleKindlmann() {}
    //--------------------------------------------------------------------------//
    virtual void init(const VecR<n> &pos) override
    {
      if (VC::vecn::norminf(this->m_pos - pos) == 0)
        return;
      // set particle onto domain border if outside
      const Domain<n> &domain = this->m_num_diff.domain();
      for (uint i = 0; i < n; ++i)
        this->m_pos[i] = std::min(std::max(pos[i], domain[i].min), domain[i].max);
      // initialize local ridge information
      ridgeinfo::kindlmann::all<n>(this->m_num_diff,
                                   this->m_pos,
                                   &this->m_ridge_strength,
                                   &this->m_ridge_tangent,
                                   &this->m_ridge_dir);
      m_ridge_dir_len = VC::vecn::norm(this->m_ridge_dir);
      if (m_ridge_dir_len > Globals::SMALL)
        this->m_ridge_dir /= m_ridge_dir_len;
    }
    //--------------------------------------------------------------------------//
    virtual bool constrain(uint max_voxel_distance = 2, std::vector<VecR<n>> *p_path = nullptr) override
    {
      // idea: use repeated init() to update pos2 and dir2 which are references to member variables
      VecR<n> pos1(this->m_pos), dir1(this->m_ridge_dir);         // pos/dir at current step of integration
      real dir1_len(m_ridge_dir_len);                             // original length of dir1
      const VecR<n> &pos2(this->m_pos), &dir2(this->m_ridge_dir); // pos/dir at forward step for step size control
      const real &dir2_len(m_ridge_dir_len);                      // original length of dir2
      VecR<n> start_pos(this->m_pos);                             // very first position, required for checking traveledVoxelDistance()

      // for debug / visualization purpose: save path positions
      auto appendCurPos = [&]()
      { if (p_path) p_path->push_back(pos1); };
      appendCurPos(); // insert first position

      // explicit Euler with step size control
      int step = 0;
      const real REL_TOL = 0.1, MAX_STEP = 0.05, MIN_STEP = 0.0000001;
      real c = 0.1 * VC::vecn::minimum(this->m_particle_container.baseGrid().cellSideLengths()); // step size
      bool accepted_last_iter = true;
      while (!accepted_last_iter || !isConverged(c))
      {
        // check whether too many steps were executed
        if (step >= 10000)
          return false;
        ++step;
        // init calculates new pos2 and dir2 (which are references to member variables)
        init(pos1 + dir1 * (dir1_len * c)); // DO NOT USE advect()! (uses pos2)

        real rel_err = VC::vecn::norm(dir1 * dir1_len - dir2 * dir2_len) / m_ridge_dir_len;
        real fac = REL_TOL / rel_err;
        accepted_last_iter = fac > 1;
        // accept current step
        if (accepted_last_iter)
        {
          pos1 = pos2;
          dir1 = dir2;
          dir1_len = dir2_len;
          appendCurPos();
        }
        // set new step size depending on fac (becomes smaller if rejected)
        c = std::max(std::min(0.9 * fac * c, MAX_STEP), MIN_STEP);

        // if voxel travel distance is limited and exceeded, return false
        if (accepted_last_iter && // check only necessary if particle was advected
            this->traveledVoxelDistance(start_pos) > max_voxel_distance)
          return false;
      }
      return this->hasMinRidgeStrength();
    }
    //--------------------------------------------------------------------------//
    real percentOfDeletedNeighbors()
    {
      std::vector<uint> neighbor_ids = this->neighborIds();
      if (neighbor_ids.empty())
        return 0.0;
      int deleted = 0;
      for (uint nid : neighbor_ids)
        if (this->m_particle_container.isDeleted(nid))
          ++deleted;
      return (real)deleted / neighbor_ids.size();
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    bool isConverged(real step_size) const
    {
      const real abs_step_len = m_ridge_dir_len * step_size;
      const real max_rel_dif = 0.001 * VC::vecn::minimum(this->m_particle_container.baseGrid().cellSideLengths());
      return abs_step_len < max_rel_dif;
    }
    //--------------------------------------------------------------------------//
    real m_ridge_dir_len;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n, typename PAR>
  class ParticleSystemBase
  {
    //--------------------------------------------------------------------------//
    using Par_t = PAR;
    using Grid_t = Grid<n>;
    using NumDiff_ptr_t = std::shared_ptr<NumericDiff<n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ParticleSystemBase(NumDiff_ptr_t p_num_diff,
                       const VecI<n> &par_init_res,
                       const VecI<n> &voxel_res,
                       real min_ridge_strength,
                       real sigma,
                       real elliptic_dis_fac)
        : mp_num_diff(p_num_diff),
          m_particle_container(voxel_res, p_num_diff->domain()),
          m_par_init_grid(par_init_res, p_num_diff->domain()),
          m_min_ridge_strength(min_ridge_strength),
          m_sigma(sigma),
          m_elliptic_dis_fac(elliptic_dis_fac)
    {
      assert(2 == n || 3 == n);
    }
    //--------------------------------------------------------------------------//
    virtual ~ParticleSystemBase() {}
    //--------------------------------------------------------------------------//
    const Domain<n> &domain() const { return mp_num_diff->domain(); }
    //--------------------------------------------------------------------------//
    const Range &domain(uint i) const { return mp_num_diff->domain(i); }
    //--------------------------------------------------------------------------//
    const VoxelContainer<n> &particleContainer() const { return m_particle_container; }
    VoxelContainer<n> &particleContainer() { return m_particle_container; }
    //--------------------------------------------------------------------------//
    void initNewParticles(bool randomize_positions = true)
    {
      for (uint cid = 0; cid < m_par_init_grid.numCells(); ++cid)
      {
        VecI<n> cell_coord;
        for (uint i = 0, div = 1; i < n; div *= m_par_init_grid.cellResolution(i++))
          cell_coord[i] = cid / div % m_par_init_grid.cellResolution(i);
        // create particle in middle of the voxel
        VecR<n> pos(m_par_init_grid.cellCenterPos(cell_coord));
        // if true, give little random offset
        if (randomize_positions)
          for (uint i = 0; i < n; ++i)
            pos[i] += ((double)rand() / RAND_MAX - 0.5) * 0.33 * m_par_init_grid.cellSideLength(i);
        m_particle_container.insert(pos);
      }
    }
    //--------------------------------------------------------------------------//
    /**
     * @brief Constrain entire population to ridges.
     *
     * @param max_voxel_distance Restricts number of voxel cells which may be
     * passed by each particle.
     * @param update_iters Defines how many times each particle should be updated
     *    depending on the local energy.
     * @param p_paths Optional parameter for debugging and visualization purposes.
     *    Consists of a pair, with first part a vector which will save the path of the particle while
     *    constraining, and second part a bool representing if the constraining was successful.
     */
    void constrainParticles(uint max_voxel_distance = 2, uint update_iters = 5,
                            std::vector<std::pair<std::vector<VecR<n>>, bool>> *p_paths = nullptr)
    {
      // if exists, prepare p_paths
      if (p_paths)
        p_paths->resize(m_particle_container.size());

      // each thread uses own particle instance (prevents repeated allocation on construction)
      Par_t par = makeParticle();
#pragma omp parallel for schedule(dynamic) firstprivate(par)
      for (uint pid = 0; pid < m_particle_container.size(); ++pid)
      {
        const VecR<n> &pos = m_particle_container.entry(pid);

        // if exists, make sure that the container is empty
        std::pair<std::vector<VecR<n>>, bool> *p_path = p_paths ? &(p_paths->at(pid)) : nullptr;
        if (p_path)
          p_path->first.resize(0);

        // init position and constrain
        par.init(pos);
        bool success = par.constrain(max_voxel_distance, p_path ? &(p_path->first) : nullptr);
        if (p_path)
          p_path->second = success;
        // apply updates and check energy only if constrain was successful
        if (success)
        {
          // update particle and find out if it improves sampling quality
          par.updateByEnergy(update_iters, EntryState::UPDATED);
          success = par.localEnergyInfluence(1.0, EntryState::UPDATED) <= 0.0;
        }
        // either update position or set delete state
        if (success)
          m_particle_container.update(pid, par.pos());
        else
          m_particle_container.setDelete(pid);

        // if (p_path)
        //   p_path->second = success;
      }
      m_particle_container.applyDelete();
      m_particle_container.unsetUpdatedStates();
    }
    //--------------------------------------------------------------------------//
    void updateParticlesByEnergy(uint update_iters = 5)
    {
      // each thread uses own particle instance (prevents repeated allocation on construction)
      Par_t par = makeParticle();
#pragma omp parallel for schedule(dynamic) firstprivate(par)
      for (uint pid = 0; pid < m_particle_container.size(); ++pid)
      {
        par.init(m_particle_container.entry(pid));
        par.updateByEnergy(update_iters);
        m_particle_container.update(pid, par.pos(), false);
      }
    }
    //--------------------------------------------------------------------------//
    void createNewNeighborParticles(uint max_iters = 100, uint update_iters = 5)
    {
      assert(2 == n || n == 3);

      // create new particles for empty areas
      auto createNewNbh = [&](const VecR<n> &pos)
      {
        Par_t new_par = this->makeParticle();
        new_par.init(pos);
        if (new_par.constrain(true))
        {
          // update particle and find out if it improves sampling quality
          new_par.updateByEnergy(update_iters);
          if (new_par.localEnergyInfluence() <= 0.0)
            m_particle_container.insert(new_par.pos());
        }
      };
      // distance of new neighbor to the currently creating particle
      const real new_par_offset = Globals::ENERGY_W * m_elliptic_dis_fac * m_sigma;

      // keep track which particles already were handled (new particles are appended in container)
      uint creation_start_idx = 0, creation_stop_idx = m_particle_container.size();
      // also count iterations for cancel criterion
      uint cur_iter = 0;
      // finally, use one particle instance per thread
      Par_t par = makeParticle();
      while (creation_start_idx < creation_stop_idx && max_iters > cur_iter++)
      {
#pragma omp parallel for schedule(dynamic) firstprivate(par)
        for (uint pid = creation_start_idx; pid < creation_stop_idx; ++pid)
        {
          // handle particle in isolation to prevent multiple creations in same location
          // double sigma to prevent spawning new particles in between them
          // m_particle_container.setIsolation(pid, 2.0 * m_sigma);

          par.init(m_particle_container.entry(pid));
          if (2 == n)
          {
            // check neighbors in both directions
            bool has_nb_forward = false;
            bool has_nb_backward = false;
            // use eigenvector of greatest eigenvalue of Hessian
            Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
            solver.compute(mp_num_diff->hessian(par.pos()));
            VecR<n> tangent(solver.eigenvectors().col(1).data());
            // look at all neighbors (with smaller range for better gap filling)
            for (const VecR<n> &other : par.neighbors(0.95))
            {
              VecR<n> dif = other - par.pos();
              // check on which side the neighbor lies by using dot product
              if ((tangent | dif) > 0.0)
                has_nb_forward = true;
              else
                has_nb_backward = true;
              if (has_nb_forward && has_nb_backward)
                break;
            }
            if (!has_nb_forward)
              createNewNbh(par.pos() + tangent * new_par_offset);
            if (!has_nb_backward)
              createNewNbh(par.pos() - tangent * new_par_offset);
          }
          else if (3 == n)
          {
            // check neighbors in both directions
            bool has_nb[6] = {false, false, false, false, false, false};
            // create two orthogonal vectors inside the ridge surface using its normal
            // EVec<3> normal(par.ridgeDir().data.data());
            // EVec<3> temp(normal.cross(EVec<3>::Random()));
            // temp.normalize();
            // VecR<n> v1(temp.data());
            // VecR<n> v2(normal.cross(temp).normalized().data());
            EVec<3> normal(par.ridgeDir().data.data());
            EVec<3> rand_vec = EVec<3>::Random();
            VecR<n> o1(normal.cross(rand_vec).normalized().data());
            VecR<n> o2(normal.cross(EVec<3>(o1.data.data())).normalized().data());
            // span 6 directions which define desired areas of neighbors
            constexpr real SIN_ANGLE = 0.866025404, COS_ANGLE = 0.5;
            std::array<VecR<n>, 6> area_directions;
            area_directions[0] = o1;
            area_directions[1] = SIN_ANGLE * o1 + COS_ANGLE * o2;
            area_directions[2] = SIN_ANGLE * o1 - COS_ANGLE * o2;
            area_directions[3] = -area_directions[0];
            area_directions[4] = -area_directions[1];
            area_directions[5] = -area_directions[2];
            // look at all neighbors (with smaller range for better gap filling)
            for (VecR<n> other : par.neighbors(0.95))
            {
              real greatest_dot = -Globals::INF;
              uint greatest_id = 0;
              for (uint area_id = 0; area_id < area_directions.size(); ++area_id)
              {
                // find correct area by finding maximum scalar product (meaning smallest angle)
                VecR<n> dif = VC::vecn::normalized(other - par.pos());
                real dot_prod = VC::vecn::dot(area_directions[area_id], dif);
                if (dot_prod > greatest_dot)
                {
                  greatest_dot = dot_prod;
                  greatest_id = area_id;
                }
              }
              // set corresponding neighbor area to true
              has_nb[greatest_id] = true;
              if (has_nb[0] && has_nb[1] && has_nb[2] && has_nb[3] && has_nb[4] && has_nb[5])
                break;
            }
            for (uint area_id = 0; area_id < area_directions.size(); ++area_id)
              if (!has_nb[area_id])
                createNewNbh(par.pos() + area_directions[area_id] * new_par_offset);
          }
          // m_particle_container.unsetIsolation(pid);
        }
        creation_start_idx = creation_stop_idx;
        creation_stop_idx = m_particle_container.size();
      }
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    Par_t makeParticle()
    {
      return Par_t(m_particle_container, *mp_num_diff, m_par_init_grid, m_sigma, m_min_ridge_strength, m_elliptic_dis_fac);
    }
    //--------------------------------------------------------------------------//
    NumDiff_ptr_t mp_num_diff;
    VoxelContainer<n> m_particle_container;
    const Grid_t m_par_init_grid;
    const real m_min_ridge_strength;
    const real m_sigma;
    const real m_elliptic_dis_fac;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleSystemStelter : public ParticleSystemBase<n, ParticleStelter<n>>
  {
    //--------------------------------------------------------------------------//
    using Par_t = ParticleStelter<n>;
    using ParSysBase_t = ParticleSystemBase<n, Par_t>;
    using ItplScField_ptr_t = std::shared_ptr<InterpolScalarField<n>>;
    using ScFieldBase_ptr_t = std::shared_ptr<ScalarFieldBase<n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    /**
     * Constructs a particle system from an InterpolScalarField which leads to
     * automatically determining the particle neighborhood parameter sigma and
     * the differences for numerical differentiations.
     */
    ParticleSystemStelter(ItplScField_ptr_t p_itpl_sc_field,
                          const VecI<n> &par_init_res,
                          const VecI<n> &voxel_res,
                          real min_ridge_strength,
                          real elliptic_dis_fac)
        : ParSysBase_t(std::make_shared<SmoothNumericDiff<n>>(p_itpl_sc_field, p_itpl_sc_field->cellSideLengths()),
                       par_init_res,
                       voxel_res,
                       min_ridge_strength,
                       VC::vecn::minimum(p_itpl_sc_field->cellSideLengths()), // sigma as smallest cell side length
                       elliptic_dis_fac)
    {
    }
    //--------------------------------------------------------------------------//
    /**
     * Constructs a particle system from an InterpolScalarField which leads to
     * automatically determining the particle neighborhood parameter sigma and
     * the differences for numerical differentiations.
     */
    ParticleSystemStelter(ItplScField_ptr_t p_itpl_sc_field,
                          const params::ParamsStelter<n> &params)
        : ParticleSystemStelter(p_itpl_sc_field,
                                params.par_init_res,
                                params.voxel_res,
                                params.min_ridge_strength,
                                params.elliptic_dis_fac) {}
    //--------------------------------------------------------------------------//
    /**
     * Constructs a particle system from a ScalarFieldBase.
     */
    ParticleSystemStelter(ScFieldBase_ptr_t p_sc_field_base,
                          const VecI<n> &par_init_res,
                          const VecI<n> &voxel_res,
                          real min_ridge_strength,
                          real sigma,
                          real elliptic_dis_fac,
                          const VecR<n> &num_diff_deltas)
        : ParSysBase_t(std::make_shared<SmoothNumericDiff<n>>(p_sc_field_base, num_diff_deltas),
                       par_init_res,
                       voxel_res,
                       min_ridge_strength,
                       sigma,
                       elliptic_dis_fac) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleSystemStelter() {}
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleSystemKindlmann : public ParticleSystemBase<n, ParticleKindlmann<n>>
  {
    //--------------------------------------------------------------------------//
    using Par_t = ParticleKindlmann<n>;
    using ParSysBase_t = ParticleSystemBase<n, Par_t>;
    using ScFieldBase_ptr_t = std::shared_ptr<ScalarFieldBase<n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    /**
     * Constructs a particle system from a ScalarFieldBase.
     */
    ParticleSystemKindlmann(ScFieldBase_ptr_t p_sc_field,
                            const VecI<n> &par_init_res,
                            const VecI<n> &voxel_res,
                            real min_ridge_strength,
                            real sigma,
                            const VecR<n> &num_diff_deltas)
        : ParSysBase_t(std::make_shared<SmoothNumericDiff<n>>(p_sc_field, num_diff_deltas),
                       par_init_res,
                       voxel_res,
                       min_ridge_strength,
                       sigma,
                       1.0) {}
    //--------------------------------------------------------------------------//
    /**
     * Constructs a particle system from a ScalarFieldBase.
     */
    ParticleSystemKindlmann(ScFieldBase_ptr_t p_sc_field,
                            const VecI<n> &par_init_res,
                            const VecI<n> &voxel_res,
                            real min_ridge_strength,
                            real sigma,
                            real numeric_diff_delta)
        : ParticleSystemKindlmann(p_sc_field,
                                  par_init_res,
                                  voxel_res,
                                  min_ridge_strength,
                                  sigma,
                                  VecR<n>(numeric_diff_delta)) {}
    //--------------------------------------------------------------------------//
    /**
     * Constructs a particle system from a ScalarFieldBase.
     */
    ParticleSystemKindlmann(ScFieldBase_ptr_t p_sc_field,
                            const params::ParamsKindlmann<n> &params)
        : ParticleSystemKindlmann(p_sc_field,
                                  params.par_init_res,
                                  params.voxel_res,
                                  params.min_ridge_strength,
                                  params.sigma,
                                  params.numeric_diff_delta) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleSystemKindlmann() {}
    //--------------------------------------------------------------------------//
    void deleteParticles()
    {
      const real sigma = this->m_sigma;
      VoxelContainer<n> &particle_container = this->m_particle_container;

      ParticleKindlmann<n> par = this->makeParticle();
#pragma omp parallel for schedule(dynamic) firstprivate(par)
      for (uint pid = 0; pid < particle_container.size(); ++pid)
      {
        // particle must be in isolation for correct info about deleted neighbors
        particle_container.setIsolation(pid, 2.0 * sigma);

        par.init(particle_container.entry(pid));
        // definitely delete particle if its feature strength is too low
        if (!par.hasMinRidgeStrength())
          particle_container.setDelete(pid);
        // also delete particle if...
        // 1. less than 50% of all neighbors are already deleted OR
        // 2. energy is higher with this particle
        else if (par.percentOfDeletedNeighbors() < 0.5 &&
                 par.localEnergyInfluence() > 0.0)
          particle_container.setDelete(pid);

        particle_container.unsetIsolation(pid);
      }
      particle_container.applyDelete();
    }
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleSystemSeriesStelter
  {
    //--------------------------------------------------------------------------//
    using ParSysStelter_t = ParticleSystemStelter<n>;
    using ItplScFieldSeries_ptr_t = std::shared_ptr<InterpolScalarFieldSeries<n>>;
    using ScFieldSeriesBase_ptr_t = std::shared_ptr<ScalarFieldSeriesBase<n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ParticleSystemSeriesStelter(ItplScFieldSeries_ptr_t p_itpl_sc_series,
                                const params::ParamsStelter<n> &params)
        : mp_sc_series_base(p_itpl_sc_series),
          m_par_sys(p_itpl_sc_series, params), // MUST use p_itpl_sc_series since the member has other smart pointer type
          m_full_init_period(params.full_init_period)
    {
    }
    //--------------------------------------------------------------------------//
    ParticleSystemSeriesStelter(ScFieldSeriesBase_ptr_t p_sc_series_base,
                                const VecI<n> &par_init_res,
                                const VecI<n> &voxel_res,
                                real min_ridge_strength,
                                real sigma,
                                real elliptic_dis_fac,
                                uint full_init_period,
                                const VecR<n> &num_diff_deltas)
        : mp_sc_series_base(p_sc_series_base),
          m_par_sys(p_sc_series_base,
                    par_init_res,
                    voxel_res,
                    min_ridge_strength,
                    sigma,
                    elliptic_dis_fac,
                    num_diff_deltas),
          m_full_init_period(full_init_period) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleSystemSeriesStelter() {}
    //--------------------------------------------------------------------------//
    void searchAndSaveRidges(const std::string &save_dir)
    {
      // create directory if necessary
      std::filesystem::create_directories(save_dir);
      // easier access to particle container + reset for fresh start
      VoxelContainer<n> &par_container = m_par_sys.particleContainer();
      par_container.reset();

      // prepare statistics
      utils::RuntimeAnalyzer analyzer("Constrain", "Create", "Update");

      // go through all time steps
      const uint num_steps = mp_sc_series_base->steppedTime().steps;
      utils::LoopPrinter printer("Search and save ridges", num_steps, false);
      for (uint time_step = 0; time_step < num_steps; ++time_step)
      {
        mp_sc_series_base->loadStep(time_step);

        // START OF ACTUAL PARTICLE SYSTEM STEPS
        bool apply_init = time_step % m_full_init_period == 0 || par_container.size() == 0 || time_step == (num_steps - 1);
        if (apply_init)
          m_par_sys.initNewParticles();

        analyzer.startTimer("Constrain");
        // increase max voxel distance if no init (otherwise weak ridges tend to disintegrate)
        m_par_sys.constrainParticles(apply_init ? 2 : 5);
        analyzer.stopTimer("Constrain");

        analyzer.startTimer("Create");
        m_par_sys.createNewNeighborParticles();
        analyzer.stopTimer("Create");

        analyzer.startTimer("Update");
        m_par_sys.updateParticlesByEnergy();
        analyzer.stopTimer("Update");

        // END OF ACTUAL PARTICLE SYSTEM STEPS

        // save particles to file
        std::string path = save_dir + "/" + std::to_string(time_step) + ".data";
        par_container.writeFile(path);

        printer.endIteration();
      }
      std::cout << "Particles (last iteration): " << m_par_sys.particleContainer().size() << std::endl;
      analyzer.logAnalysis();
    }
    //--------------------------------------------------------------------------//
    const ParSysStelter_t &particleSystem() const { return m_par_sys; }
    ParSysStelter_t &particleSystem() { return m_par_sys; }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScFieldSeriesBase_ptr_t mp_sc_series_base;
    ParSysStelter_t m_par_sys;
    uint m_full_init_period;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class ParticleSystemSeriesKindlmann
  {
    //--------------------------------------------------------------------------//
    using ScFieldSeriesBase_ptr_t = std::shared_ptr<ScalarFieldSeriesBase<n>>;
    using ParSysKindlmann_t = ParticleSystemKindlmann<n>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    ParticleSystemSeriesKindlmann(ScFieldSeriesBase_ptr_t p_sc_field_series,
                                  const params::ParamsKindlmann<n> &params)
        : mp_sc_field_series(p_sc_field_series),
          m_par_sys(mp_sc_field_series, params) {}
    //--------------------------------------------------------------------------//
    virtual ~ParticleSystemSeriesKindlmann() {}
    //--------------------------------------------------------------------------//
    void searchAndSaveRidges(const std::string &save_dir)
    {
      // create directory if necessary
      std::filesystem::create_directories(save_dir);
      // easier access to particle container
      VoxelContainer<n> &par_container = m_par_sys.particleContainer();

      // prepare statistics
      utils::RuntimeAnalyzer analyzer("Constrain", "Update", "Create", "Delete");

      // go through all time steps
      uint num_steps = mp_sc_field_series->steppedTime().steps;
      utils::LoopPrinter printer("Search and save ridges", num_steps, false);
      for (uint time_step = 0; time_step < num_steps; ++time_step)
      {
        // setup fixed value for time dimension
        mp_sc_field_series->loadStep(time_step);

        // START OF ACTUAL PARTICLE SYSTEM STEPS

        // par_container.reset();
        m_par_sys.initNewParticles();

        analyzer.startTimer("Constrain");
        m_par_sys.constrainParticles(2, 0); // max_voxel_distance = 2 + no automatic updates
        analyzer.stopTimer("Constrain");

        int iter = 0;
        bool is_pop_growing = true;
        // iterate until either too many iterations passed or pop size is no longer growing
        while (101 > iter && is_pop_growing)
        {
          // apply population control every PC iterations
          if (iter % PC == 0)
          {
            int last_pop_size = par_container.size();
            analyzer.startTimer("Delete");
            m_par_sys.deleteParticles();
            analyzer.stopTimer("Delete");

            analyzer.startTimer("Create");
            m_par_sys.createNewNeighborParticles(1, PC); // update new particles PC times
            analyzer.stopTimer("Create");
            int new_pop_size = par_container.size();
            // track growing state
            is_pop_growing = last_pop_size < new_pop_size;
          }
          analyzer.startTimer("Update");
          m_par_sys.updateParticlesByEnergy(1); // one single update step
          analyzer.stopTimer("Update");

          analyzer.startTimer("Constrain");
          m_par_sys.constrainParticles(2, 0); // max_voxel_distance = 2 + no automatic updates
          analyzer.stopTimer("Constrain");

          ++iter;
        }
        // END OF ACTUAL PARTICLE SYSTEM STEPS

        // save particles to file
        std::string path = save_dir + "/" + std::to_string(time_step) + ".data";
        par_container.writeFile(path);

        printer.endIteration();
      }
      std::cout << "Particles (last iteration): " << m_par_sys.particleContainer().size() << std::endl;
      analyzer.logAnalysis();
    }
    //--------------------------------------------------------------------------//
    const ParSysKindlmann_t &particleSystem() const { return m_par_sys; }
    ParSysKindlmann_t &particleSystem() { return m_par_sys; }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScFieldSeriesBase_ptr_t mp_sc_field_series;
    ParSysKindlmann_t m_par_sys;
    //--------------------------------------------------------------------------//
    const int PC = 10; // periodicity for population control
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//
