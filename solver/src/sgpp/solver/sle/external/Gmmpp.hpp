// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/solver/sle/common/MySLESolver.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace solver {

/**
 * Linear system solver using Gmm++ (iterative sparse solver).
 */
class Gmmpp : public MySLESolver {
 public:
  /**
   * Destructor.
   */
  ~Gmmpp() override;

  /**
   * @param       system  system to be solved
   * @param       b       right-hand side
   * @param[out]  x       solution to the system
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  bool solve(base::SLE& system, base::DataVector& b, base::DataVector& x) const override;
};
}  // namespace solver
}  // namespace sgpp
