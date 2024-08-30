#pragma once
//--------------------------------------------------------------------------//
#include "../globals.hpp"
//--------------------------------------------------------------------------//
#include <filesystem>
#include <fstream>
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  template <uint n>
  class Grid
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    Grid() : Grid(VecI<n>(0u), Domain<n>()) {}
    //--------------------------------------------------------------------------//
    Grid(const VecI<n> &cell_resolution, const Domain<n> &domain)
    {
      redefine(cell_resolution, domain);
    }
    //--------------------------------------------------------------------------//
    virtual ~Grid() {}
    //--------------------------------------------------------------------------//
    const VecI<n> &cellResolution() const { return m_cell_resolution; }
    uint cellResolution(uint i) const { return m_cell_resolution[i]; }
    //--------------------------------------------------------------------------//
    const VecI<n> &cornerResolution() const { return m_corner_resolution; }
    uint cornerResolution(uint i) const { return m_corner_resolution[i]; }
    //--------------------------------------------------------------------------//
    const Domain<n> &domain() const { return m_domain; }
    const Range &domain(uint i) const { return m_domain[i]; }
    //--------------------------------------------------------------------------//
    VecR<n> cellSideLengths() const { return m_edge_len; }
    real cellSideLength(uint i) const { return m_edge_len[i]; }
    //--------------------------------------------------------------------------//
    uint numCells() const { return VC::vecn::prod(m_cell_resolution); }
    //--------------------------------------------------------------------------//
    uint numCorners() const { return VC::vecn::prod(m_corner_resolution); }
    //--------------------------------------------------------------------------//
    bool isCellCoordInside(const VecI<n> &cell_coord) const
    {
      for (uint i = 0; i < n; ++i)
        if (cell_coord[i] >= m_cell_resolution[i])
          return false;
      return true;
    }
    //--------------------------------------------------------------------------//
    bool isCornerCoordInside(const VecI<n> &corner_coord) const
    {
      for (uint i = 0; i < n; ++i)
        if (corner_coord[i] >= m_corner_resolution[i])
          return false;
      return true;
    }
    //--------------------------------------------------------------------------//
    VecR<n> cellMinPos(const VecI<n> &cell_coord) const
    {
      assert(isCellCoordInside(cell_coord));
      VecR<n> result;
      for (uint i = 0; i < n; ++i)
        result[i] = m_domain[i].min + (m_edge_len[i] * cell_coord[i]);
      return result;
    }
    //--------------------------------------------------------------------------//
    VecR<n> cellMaxPos(const VecI<n> &cell_coord) const
    {
      assert(isCellCoordInside(cell_coord));
      VecR<n> result;
      for (uint i = 0; i < n; ++i)
        result[i] = m_domain[i].min + (m_edge_len[i] * (cell_coord[i] + 1));
      return result;
    }
    //--------------------------------------------------------------------------//
    VecR<n> cellCenterPos(const VecI<n> &cell_coord) const
    {
      return (cellMinPos(cell_coord) + cellMaxPos(cell_coord)) * 0.5;
    }
    //--------------------------------------------------------------------------//
    VecR<n> cornerCoordToPos(const VecI<n> &corner_coord) const
    {
      assert(isCornerCoordInside(corner_coord));
      VecR<n> result;
      for (uint i = 0; i < n; ++i)
        result[i] = m_domain[i].min + (m_edge_len[i] * corner_coord[i]);
      return result;
    }
    //--------------------------------------------------------------------------//
    VecR<n> cornerIndexToPos(uint corner_idx) const
    {
      return cornerCoordToPos(cornerIndexToCoord(corner_idx));
    }
    //--------------------------------------------------------------------------//
    VecI<n> posToCellCoord(const VecR<n> &pos) const
    {
      assert(m_domain.isInside(pos));
      VecI<n> cell_coord;
      for (uint i = 0; i < n; ++i)
        cell_coord[i] = std::min(m_cell_resolution[i] - 1,
                                 (uint)((pos[i] - m_domain[i].min) * m_rec_edge_len[i]));
      return cell_coord;
    }
    //--------------------------------------------------------------------------//
    uint cornerCoordToIndex(const VecI<n> &corner_coord) const
    {
      uint corner_idx = 0;
      for (uint i = 0, mult = 1; i < n; mult *= this->cornerResolution(i++))
        corner_idx += corner_coord[i] * mult;
      return corner_idx;
    }
    //--------------------------------------------------------------------------//
    VecI<n> cornerIndexToCoord(uint corner_idx) const
    {
      VecI<n> coord;
      for (uint i = 0, div = 1; i < n; div *= this->cornerResolution(i++))
        coord[i] = corner_idx / div % this->cornerResolution(i);
      return coord;
    }
    //--------------------------------------------------------------------------//
    uint cellCoordToIndex(const VecI<n> &cell_coord) const
    {
      uint cell_idx = 0;
      for (uint i = 0, mult = 1; i < n; mult *= this->cellResolution(i++))
        cell_idx += cell_coord[i] * mult;
      return cell_idx;
    }
    //--------------------------------------------------------------------------//
    VecI<n> cellIndexToCoord(uint cell_idx) const
    {
      VecI<n> coord;
      for (uint i = 0, div = 1; i < n; div *= this->cellResolution(i++))
        coord[i] = cell_idx / div % this->cellResolution(i);
      return coord;
    }
    //--------------------------------------------------------------------------//
    uint edgeCoordToIndex(const VecI<n + 1> &edge_coord) const
    {
      VecI<n> corner_coord = VC::vecn::first<n>(edge_coord);
      uint edge_dim = VC::vecn::last(edge_coord);
      return cornerCoordToIndex(corner_coord) * 3 + edge_dim;
    }
    //--------------------------------------------------------------------------//
    VecI<n + 1> edgeIndexToCoord(uint edge_idx) const
    {
      VecI<n + 1> coord;
      coord[n] = edge_idx % 3;
      for (uint i = 0, div = 3; i < n; div *= this->cornerResolution(i++))
        coord[i] = edge_idx / div % this->cornerResolution(i);
      return coord;
    }
    //--------------------------------------------------------------------------//
    VecR<n> edgeMinPos(const VecI<n + 1> &edge_coord) const
    {
      VecI<n> corner_coord = VC::vecn::first<n>(edge_coord);
      return cornerCoordToPos(corner_coord);
    }
    //--------------------------------------------------------------------------//
    VecR<n> edgeMaxPos(const VecI<n + 1> &edge_coord) const
    {
      VecI<n> corner_coord = VC::vecn::first<n>(edge_coord);
      corner_coord[edge_coord[n]] += 1;
      return cornerCoordToPos(corner_coord);
    }
    //--------------------------------------------------------------------------//
    int cellCoordDistance(const VecR<n> &pos1, const VecR<n> &pos2) const
    {
      int distance = 0;
      VecI<n> dist_vec = posToCellCoord(pos1) - posToCellCoord(pos2);
      for (uint i = 0; i < n; ++i)
        distance += std::abs((int)dist_vec[i]);
      return distance;
    }
    //--------------------------------------------------------------------------//
    std::vector<VecI<n>> cellNeighborCoords(const VecI<n> &cell_coord, uint neighbor_size = 1) const
    {
      std::vector<VecI<n>> coords;
      uint range = 1 + 2 * neighbor_size;
      uint num_neighbors = pow(range, n);
      coords.reserve(num_neighbors);
      // find all neighbors (also diagonal ones) -> create greater square (2D) / cube (3D) / ...
      for (uint power_id = 0; power_id < num_neighbors; ++power_id)
      {
        // 'unpack' power_id to corresponding cell id
        VecI<n> neighbor_cell_id(cell_coord);
        for (uint i = 0, div = 1; i < n; ++i, div *= range)
          neighbor_cell_id[i] += (power_id / div % range) - neighbor_size; // expression is inside [-neighbor_size, neighbor_size]

        if (isCellCoordInside(neighbor_cell_id))
          coords.push_back(neighbor_cell_id);
      }
      return coords;
    }
    //--------------------------------------------------------------------------//
    std::vector<VecI<n>> cellNeighborCoords(const VecR<n> &pos, real range) const
    {
      VecI<n> own_cell = posToCellCoord(pos);
      std::vector<VecI<n>> coords;
      // search for all roughly nearby cells
      uint max_neighbor_size = 0;
      for (uint i = 0; i < n; ++i)
        max_neighbor_size = std::max(max_neighbor_size, (uint)(range * m_rec_edge_len[i]) + 1);
      // filter by finding nearest distances between pos and cells
      for (const VecI<n> &nbh_cell : cellNeighborCoords(own_cell, max_neighbor_size))
      {
        VecR<n> min = cellMinPos(nbh_cell), max = cellMaxPos(nbh_cell), nearest_pos;
        for (uint i = 0, div = 1; i < n; ++i, div *= range)
          nearest_pos[i] = nbh_cell[i] == own_cell[i] ? pos[i] : nbh_cell[i] < own_cell[i] ? max[i]
                                                                                           : min[i];
        if (VC::vecn::norm(pos - nearest_pos) <= range)
          coords.push_back(nbh_cell);
      }
      return coords;
    }
    //--------------------------------------------------------------------------//
    virtual void redefine(const VecI<n> &cell_resolution, const Domain<n> &domain)
    {
      m_cell_resolution = cell_resolution;
      m_corner_resolution = cell_resolution + VecI<n>(1);
      m_domain = domain;
      m_domain_len = m_domain.difs();
      for (uint i = 0; i < n; ++i)
      {
        m_edge_len[i] = m_domain_len[i] / m_cell_resolution[i];
        m_rec_edge_len[i] = m_cell_resolution[i] / m_domain_len[i];
      }
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    VecI<n> m_cell_resolution;
    VecI<n> m_corner_resolution;
    Domain<n> m_domain;
    VecR<n> m_domain_len;
    VecR<n> m_edge_len;
    VecR<n> m_rec_edge_len;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n, typename T>
  class CornerDataGrid : public Grid<n>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    CornerDataGrid() : Grid<n>(), m_data() {}
    //--------------------------------------------------------------------------//
    CornerDataGrid(const VecI<n> &cell_resolution,
                   const Domain<n> &domain)
        : Grid<n>(cell_resolution, domain),
          m_data(this->numCorners()) {}
    //--------------------------------------------------------------------------//
    CornerDataGrid(const VecI<n> &cell_resolution,
                   const Domain<n> &domain,
                   const T &value)
        : Grid<n>(cell_resolution, domain),
          m_data(this->numCorners(), value)
    {
      assert(this->numCorners() == m_data.size());
    }
    //--------------------------------------------------------------------------//
    CornerDataGrid(const VecI<n> &cell_resolution,
                   const Domain<n> &domain,
                   const std::vector<T> &data)
        : Grid<n>(cell_resolution, domain),
          m_data(data)
    {
      assert(this->numCorners() == data.size());
    }
    //--------------------------------------------------------------------------//
    CornerDataGrid(const VecI<n> &cell_resolution,
                   const Domain<n> &domain,
                   std::vector<T> &&data)
        : Grid<n>(cell_resolution, domain),
          m_data(data)
    {
      assert(this->numCorners() == data.size());
    }
    //--------------------------------------------------------------------------//
    virtual ~CornerDataGrid() {}
    //--------------------------------------------------------------------------//
    const std::vector<T> &cornerData() const { return m_data; }
    std::vector<T> &cornerData() { return m_data; }
    //--------------------------------------------------------------------------//
    const T &cornerData(uint corner_idx) const { return m_data[corner_idx]; }
    T &cornerData(uint corner_idx) { return m_data[corner_idx]; }
    //--------------------------------------------------------------------------//
    const T &cornerData(const VecI<n> &corner_coord) const { return m_data[this->cornerCoordToIndex(corner_coord)]; }
    T &cornerData(const VecI<n> &corner_coord) { return m_data[this->cornerCoordToIndex(corner_coord)]; }
    //--------------------------------------------------------------------------//
    virtual void redefine(const VecI<n> &cell_resolution, const Domain<n> &domain) override
    {
      Grid<n>::redefine(cell_resolution, domain);
      m_data.resize(this->numCorners());
    }
    //--------------------------------------------------------------------------//
    virtual void writeFile(const std::string &file_path) const
    {
      using namespace std;

      filesystem::create_directories(filesystem::absolute(file_path).parent_path());
      ofstream out(file_path, ios::binary);
      if (!out)
        throw runtime_error("CornerDataGrid: could not write to file '" + file_path + "'");

      // write dimension, cell resolution, domain and stepped time
      uint dim = n;
      out.write((char *)&dim, sizeof(dim));
      out.write((char *)&this->cellResolution(), sizeof(this->cellResolution()));
      out.write((char *)&this->domain(), sizeof(this->domain()));
      // write data section
      out.write((char *)&m_data[0], this->numCorners() * sizeof(T));
      out.close();
    }
    //--------------------------------------------------------------------------//
    virtual void readFile(const std::string &file_path)
    {
      using namespace std;

      ifstream in(file_path, ios::binary);
      if (!in)
        throw runtime_error("CornerDataGrid: could not read from file '" + file_path + "'");

      // check dimension
      uint dim;
      in.read((char *)&dim, sizeof(dim));
      if (dim != n)
        throw runtime_error("CornerDataGrid: expected dimension " + to_string(n) + ", but was " + to_string(dim));
      // get cell resolution and domain
      VecI<n> cell_res;
      Domain<n> domain;
      in.read((char *)&cell_res, sizeof(cell_res));
      in.read((char *)&domain, sizeof(domain));
      // apply meta data which also resizes the data vector
      redefine(cell_res, domain);
      // get data section
      in.read((char *)&m_data[0], this->numCorners() * sizeof(T));
      in.close();
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    std::vector<T> m_data;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n, typename T>
  class CellDataGrid : public Grid<n>
  {
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    CellDataGrid() : Grid<n>(), m_data() {}
    //--------------------------------------------------------------------------//
    CellDataGrid(const VecI<n> &cell_resolution,
                 const Domain<n> &domain)
        : Grid<n>(cell_resolution, domain),
          m_data(this->numCells()) {}
    //--------------------------------------------------------------------------//
    CellDataGrid(const VecI<n> &cell_resolution,
                 const Domain<n> &domain,
                 const T &value)
        : Grid<n>(cell_resolution, domain),
          m_data(this->numCells(), value)
    {
      assert(this->numCells() == m_data.size());
    }
    //--------------------------------------------------------------------------//
    CellDataGrid(const VecI<n> &cell_resolution,
                 const Domain<n> &domain,
                 const std::vector<T> &data)
        : Grid<n>(cell_resolution, domain),
          m_data(data)
    {
      assert(this->numCells() == data.size());
    }
    //--------------------------------------------------------------------------//
    CellDataGrid(const VecI<n> &cell_resolution,
                 const Domain<n> &domain,
                 std::vector<T> &&data)
        : Grid<n>(cell_resolution, domain),
          m_data(data)
    {
      assert(this->numCells() == data.size());
    }
    //--------------------------------------------------------------------------//
    virtual ~CellDataGrid() {}
    //--------------------------------------------------------------------------//
    const std::vector<T> &cellData() const { return m_data; }
    std::vector<T> &cellData() { return m_data; }
    //--------------------------------------------------------------------------//
    const T &cellData(uint cell_idx) const { return m_data[cell_idx]; }
    T &cellData(uint cell_idx) { return m_data[cell_idx]; }
    //--------------------------------------------------------------------------//
    const T &cellData(const VecI<n> &cell_coord) const { return m_data[this->cellCoordToIndex(cell_coord)]; }
    T &cellData(const VecI<n> &cell_coord) { return m_data[this->cellCoordToIndex(cell_coord)]; }
    //--------------------------------------------------------------------------//
    virtual void redefine(const VecI<n> &cell_resolution, const Domain<n> &domain) override
    {
      Grid<n>::redefine(cell_resolution, domain);
      m_data.resize(this->numCells());
    }
    //--------------------------------------------------------------------------//
    virtual void writeFile(const std::string &file_path) const
    {
      std::filesystem::create_directories(std::filesystem::absolute(file_path).parent_path());
      std::ofstream out(file_path, std::ios::binary);
      if (!out)
        throw std::runtime_error("CellDataGrid: could not write to file '" + file_path + "'");

      // META DATA (ASCII)
      // write dimension and cell resolution
      out << n << ' ' << this->m_cell_resolution;
      // write domain
      for (uint i = 0; i < n; ++i)
        out << ' ' << this->m_domain[i].min << ' ' << this->m_domain[i].max;
      // delimiter for data section
      out << '\n';

      // VECTOR DATA (BINARY)
      out.write((char *)&m_data[0], this->numCells() * sizeof(T));
      out.close();
    }
    //--------------------------------------------------------------------------//
    virtual void readFile(const std::string &file_path)
    {
      std::ifstream in(file_path, std::ios::binary);
      if (!in)
        throw std::runtime_error("CellDataGrid: could not read from file '" + file_path + "'");

      // META DATA (ASCII)
      // check dimension
      uint dimension;
      in >> dimension;
      if (dimension != n)
        throw std::runtime_error("CornerDataGrid: expected dimension " + std::to_string(n) + ", but was " + std::to_string(dimension));
      // get cell resolution
      VecI<n> cell_res;
      in >> cell_res;
      // get domain
      Domain<n> domain;
      for (uint i = 0; i < n; ++i)
        in >> domain[i].min >> domain[i].max;
      // apply meta data which also resizes the data vector
      redefine(cell_res, domain);
      // skip delimiter '\n'
      in.ignore(1);

      // VECTOR DATA (BINARY)
      in.read((char *)&m_data[0], this->numCells() * sizeof(T));
      in.close();
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    std::vector<T>
        m_data;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//