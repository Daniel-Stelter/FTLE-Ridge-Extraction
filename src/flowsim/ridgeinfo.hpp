#pragma once
//--------------------------------------------------------------------------//
#include "numericdiff.hpp"
//--------------------------------------------------------------------------//
namespace dst::flowsim::ridgeinfo
{
  //--------------------------------------------------------------------------//
  namespace kindlmann
  {
    //--------------------------------------------------------------------------//
    template <uint n>
    real strength(const NumericDiff<n> &num_diff, const VecR<n> &pos)
    {
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(num_diff.hessian(pos));
      return -solver.eigenvalues()(0);
    }
    //--------------------------------------------------------------------------//
    template <uint n>
    EMat<n, n> tangent(const NumericDiff<n> &num_diff, const VecR<n> &pos)
    {
      // use eigenvalues of Hessian matrix
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(num_diff.hessian(pos));
      EMat<n, n> tangent = Eigen::MatrixXd::Zero(n, n);
      for (uint i = 1; i < n; ++i)
        tangent += solver.eigenvectors().col(i) * solver.eigenvectors().col(i).transpose();
      return tangent;
    }
    //--------------------------------------------------------------------------//
    template <uint n>
    VecR<n> direction(const NumericDiff<n> &num_diff, const VecR<n> &pos)
    {
      static const EMat<n, n> I = Eigen::MatrixXd::Identity(n, n);
      return VecR<n>(EVec<n>((I - tangent<n>(num_diff, pos)) * EVec<n>(num_diff.grad(pos).data.data())).data());
    }
    //--------------------------------------------------------------------------//
    template <uint n>
    void all(const NumericDiff<n> &num_diff,
             const VecR<n> &pos,
             real *p_strength = nullptr,
             EMat<n, n> *p_tangent = nullptr,
             VecR<n> *p_dir = nullptr)
    {
      // use eigenvalues of Hessian matrix
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(num_diff.hessian(pos));

      // RIDGE STRENGTH
      if (p_strength)
        *p_strength = -solver.eigenvalues()(0);

      // TANGENT || RIDGE DIRECTION
      if (p_tangent || p_dir)
      {
        // either directly use p_tangent or use a temporary matrix
        EMat<n, n> tmp_tangent;
        EMat<n, n> *p_used_tangent = p_tangent ? p_tangent : &tmp_tangent;
        *p_used_tangent = Eigen::MatrixXd::Zero(n, n);
        for (uint i = 1; i < n; ++i)
          *p_used_tangent += solver.eigenvectors().col(i) * solver.eigenvectors().col(i).transpose();
        // RIDGE DIRECTION
        if (p_dir)
        {
          static const EMat<n, n> I = Eigen::MatrixXd::Identity(n, n);
          *p_dir = VecR<n>(EVec<n>((I - *p_used_tangent) * EVec<n>(num_diff.grad(pos).data.data())).data());
        }
      }
    }
    //--------------------------------------------------------------------------//
  }
  //--------------------------------------------------------------------------//
  namespace stelter
  {
    //--------------------------------------------------------------------------//
    template <uint n>
    real strength(const NumericDiff<n> &num_diff, const VecR<n> &pos)
    {
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(num_diff.hessian(pos));
      return -solver.eigenvalues()(0);
    }
    //--------------------------------------------------------------------------//
    template <uint n>
    EMat<n, n> tangent(const NumericDiff<n> &num_diff, const VecR<n> &pos)
    {
      // use eigenvalues of Hessian matrix
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(num_diff.hessian(pos));
      // compare smallest and greatest eigenvalues
      real abs_eig_val_dif = abs(solver.eigenvalues()[n - 1]) - abs(solver.eigenvalues()[0]);

      EMat<n, n> tangent = Eigen::MatrixXd::Zero(n, n);
      if (abs(abs_eig_val_dif) > Globals::ZERO)
      {
        // use all but the eigenvector of greatest magnitude (with most extreme eigenvalue) for the tangent
        if (abs_eig_val_dif > 0)
          for (uint i = 0; i < n - 1; ++i)
            tangent += solver.eigenvectors().col(i) * solver.eigenvectors().col(i).transpose();
        else
          for (uint i = 1; i < n; ++i)
            tangent += solver.eigenvectors().col(i) * solver.eigenvectors().col(i).transpose();
      }
      return tangent;
    }
    //--------------------------------------------------------------------------//
    template <uint n>
    VecR<n> direction(const NumericDiff<n> &num_diff, const VecR<n> &pos)
    {
      static const EMat<n, n> I = Eigen::MatrixXd::Identity(n, n);
      return VecR<n>(EVec<n>((I - tangent<n>(num_diff, pos)) * EVec<n>(num_diff.grad(pos).data.data())).data());
    }
    //--------------------------------------------------------------------------//
    template <uint n>
    void all(const NumericDiff<n> &num_diff,
             const VecR<n> &pos,
             real *p_strength = nullptr,
             EMat<n, n> *p_tangent = nullptr,
             VecR<n> *p_dir = nullptr)
    {
      // use eigenvalues of Hessian matrix
      Eigen::SelfAdjointEigenSolver<EMat<n, n>> solver;
      solver.computeDirect(num_diff.hessian(pos));
      // compare smallest and greatest eigenvalues
      real abs_eig_val_dif = abs(solver.eigenvalues()[n - 1]) - abs(solver.eigenvalues()[0]);

      // RIDGE STRENGTH
      if (p_strength)
        *p_strength = -solver.eigenvalues()(0);

      // TANGENT || RIDGE DIRECTION
      if (p_tangent || p_dir)
      {
        // either directly use p_tangent or use a temporary matrix
        EMat<n, n> tmp_tangent;
        EMat<n, n> *p_used_tangent = p_tangent ? p_tangent : &tmp_tangent;
        *p_used_tangent = Eigen::MatrixXd::Zero(n, n);
        if (abs(abs_eig_val_dif) > Globals::ZERO)
        {
          // use all but the eigenvector of greatest magnitude (with most extreme eigenvalue) for the tangent
          if (abs_eig_val_dif > 0)
            for (uint i = 0; i < n - 1; ++i)
              *p_used_tangent += solver.eigenvectors().col(i) * solver.eigenvectors().col(i).transpose();
          else
            for (uint i = 1; i < n; ++i)
              *p_used_tangent += solver.eigenvectors().col(i) * solver.eigenvectors().col(i).transpose();
        }
        // RIDGE DIRECTION
        if (p_dir)
        {
          static const EMat<n, n> I = Eigen::MatrixXd::Identity(n, n);
          *p_dir = VecR<n>(EVec<n>((I - *p_used_tangent) * EVec<n>(num_diff.grad(pos).data.data())).data());
        }
      }
    }
    //--------------------------------------------------------------------------//
  }
  //--------------------------------------------------------------------------//
}
//--------------------------------------------------------------------------//