#pragma once
//--------------------------------------------------------------------------//
#include "scalarfield.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim
{
  //--------------------------------------------------------------------------//
  template <uint n>
  class NumericDiff
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using ScFieldBase_ptr_t = std::shared_ptr<ScalarFieldBase<n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    NumericDiff(ScFieldBase_ptr_t p_scalar, real delta)
        : NumericDiff(p_scalar, VecR<n>(delta)) {}
    //--------------------------------------------------------------------------//
    NumericDiff(ScFieldBase_ptr_t p_scalar, const VecR<n> &deltas)
        : mp_scalar(p_scalar),
          m_deltas(deltas) {}
    //--------------------------------------------------------------------------//
    virtual ~NumericDiff() {}
    //--------------------------------------------------------------------------//
    const Domain<n> &domain() const { return mp_scalar->domain(); }
    const Range &domain(uint i) const { return mp_scalar->domain(i); }
    //--------------------------------------------------------------------------//
    const VecR<n> &deltas() const { return m_deltas; };
    real delta(uint i) const { return m_deltas[i]; };
    //--------------------------------------------------------------------------//
    virtual VecR<n> grad(const VecR<n> &pos) const
    {
      VecR<n> gradient;
      for (uint i = 0; i < n; ++i)
      {
        // depending on being too near to the border, use central, forward or backward differences
        VecR<n> pos_m1(pos), pos_p1(pos);
        if (pos[i] - m_deltas[i] >= mp_scalar->domain(i).min)
          pos_m1[i] -= m_deltas[i];
        if (pos[i] + m_deltas[i] <= mp_scalar->domain(i).max)
          pos_p1[i] += m_deltas[i];
        // in case of central differences: h_i is already already doubled
        real h_i = pos_p1[i] - pos_m1[i];
        // central:  (g f(x)) / (d x_i) = (f(x + h_i e_i) - f(x - h_i e_i)) / 2h_i
        // forward:  (g f(x)) / (d x_i) = (f(x + h_i e_i) - f(x)) / h_i
        // backward: (g f(x)) / (d x_i) = (f(x) - f(x - h_i e_i)) / h_i
        gradient[i] = (mp_scalar->value(pos_p1) - mp_scalar->value(pos_m1)) / h_i;
      }
      return gradient;
    }
    //--------------------------------------------------------------------------//
    virtual EMat<n, n> hessian(const VecR<n> &pos) const
    {
      EMat<n, n> H;
      for (uint i = 0; i < n; ++i)
      {
        real h_i = m_deltas[i];
        // depending on being too near to the border, use central, forward or backward differences
        VecR<n> pos_m1(pos), pos_p1(pos);
        if (pos[i] - h_i >= mp_scalar->domain(i).min)
          pos_m1[i] -= h_i;
        if (pos[i] + h_i <= mp_scalar->domain(i).max)
          pos_p1[i] += h_i;
        // in case of central differences: delta is already already double of h_i - times 2 is correct anyway
        real delta = 2 * (pos_p1[i] - pos_m1[i]);
        H.col(i) = EVec<n>((grad(pos_p1) - grad(pos_m1)).data.data()) / delta;
      }
      return (H + H.transpose()) * 0.5;
    }
    //--------------------------------------------------------------------------//
  protected:
    //--------------------------------------------------------------------------//
    ScFieldBase_ptr_t mp_scalar;
    VecR<n> m_deltas;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//

  //--------------------------------------------------------------------------//
  template <uint n>
  class SmoothNumericDiff : public NumericDiff<n>
  //--------------------------------------------------------------------------//
  {
    //--------------------------------------------------------------------------//
    using ScFieldBase_ptr_t = std::shared_ptr<ScalarFieldBase<n>>;
    //--------------------------------------------------------------------------//
  public:
    //--------------------------------------------------------------------------//
    SmoothNumericDiff(ScFieldBase_ptr_t p_scalar, real delta)
        : NumericDiff<n>(p_scalar, VecR<n>(delta)) {}
    //--------------------------------------------------------------------------//
    SmoothNumericDiff(ScFieldBase_ptr_t p_scalar, const VecR<n> &deltas)
        : NumericDiff<n>(p_scalar, deltas) {}
    //--------------------------------------------------------------------------//
    virtual ~SmoothNumericDiff() {}
    //--------------------------------------------------------------------------//
    virtual VecR<n> grad(const VecR<n> &pos) const override
    {
      Eigen::MatrixXd A;
      Eigen::VectorXd b;
      setupGradSystem(pos, A, b);
      Eigen::VectorXd weights = A.completeOrthogonalDecomposition().solve(b);
      VecR<n> result;
      for (uint i = 0; i < n; ++i)
        result[i] = weights(i + 1);
      return result;
    }
    //--------------------------------------------------------------------------//
    void setupGradSystem(const VecR<n> &pos, Eigen::MatrixXd &A, Eigen::VectorXd &b) const;
    //--------------------------------------------------------------------------//
  };
  //--------------------------------------------------------------------------//
  template <>
  inline void SmoothNumericDiff<2>::setupGradSystem(const VecR<2> &pos,
                                                 Eigen::MatrixXd &A,
                                                 Eigen::VectorXd &b) const
  {
    A.resize(13, 3);
    b.resize(13);
    VecR<2> temp(pos);
    int line_id = 0;
    for (int i = -2; i <= 2; ++i)
    {
      temp[0] = pos[0] + i * this->m_deltas[0];
      if (temp[0] < this->mp_scalar->domain(0).min || temp[0] > this->mp_scalar->domain(0).max)
        continue;
      for (int j = -2 + abs(i); j <= 2 - abs(i); ++j)
      {
        temp[1] = pos[1] + j * this->m_deltas[1];
        if (temp[1] < this->mp_scalar->domain(1).min || temp[1] > this->mp_scalar->domain(1).max)
          continue;
        A.row(line_id) = Eigen::RowVector3d{1, temp[0], temp[1]};
        b(line_id) = this->mp_scalar->value(temp);
        ++line_id;
      }
    }
    A.conservativeResize(line_id, Eigen::NoChange);
    b.conservativeResize(line_id);
  }
  //--------------------------------------------------------------------------//
  template <>
  inline void SmoothNumericDiff<3>::setupGradSystem(const VecR<3> &pos,
                                                 Eigen::MatrixXd &A,
                                                 Eigen::VectorXd &b) const
  {
    A.resize(25, 4);
    b.resize(25);
    VecR<3> temp(pos);
    int line_id = 0;
    for (int i = -2; i <= 2; ++i)
    {
      temp[0] = pos[0] + i * this->m_deltas[0];
      if (temp[0] < this->mp_scalar->domain(0).min || temp[0] > this->mp_scalar->domain(0).max)
        continue;
      for (int j = -2 + abs(i); j <= 2 - abs(i); ++j)
      {
        temp[1] = pos[1] + j * this->m_deltas[1];
        if (temp[1] < this->mp_scalar->domain(1).min || temp[1] > this->mp_scalar->domain(1).max)
          continue;
        for (int k = -2 + abs(i) + abs(j); k <= 2 - abs(i) - abs(j); ++k)
        {
          temp[2] = pos[2] + k * this->m_deltas[2];
          if (temp[2] < this->mp_scalar->domain(2).min || temp[2] > this->mp_scalar->domain(2).max)
            continue;
          A.row(line_id) = Eigen::RowVector4d{1, temp[0], temp[1], temp[2]};
          b(line_id) = this->mp_scalar->value(temp);
          ++line_id;
        }
      }
    }
    A.conservativeResize(line_id, Eigen::NoChange);
    b.conservativeResize(line_id);
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//