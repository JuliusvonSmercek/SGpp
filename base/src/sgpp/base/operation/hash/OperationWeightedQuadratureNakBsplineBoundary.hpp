// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationWeightedQuadrature.hpp>
#include <sgpp/globaldef.hpp>
#include "common/basis/NakBsplineBoundaryBasis.hpp"

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, not a knot B-spline grid created by transformation of a not a knot
 * B-spline combigrid
 */
class OperationWeightedQuadratureNakBsplineBoundary : public OperationWeightedQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakBsplineBoundary
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationWeightedQuadratureNakBsplineBoundary(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationWeightedQuadratureNakBsplineBoundary() override {}

  /**
   * Quadrature for not a knot B-spline basis functions w.r.t. a probability density function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdf			probability density function
   * @parm quadOrder	order for the gauss Legendre quadrature
   */
  double doWeightedQuadrature(DataVector& alpha, sgpp::base::DistributionsVector pdfs,
                              size_t quadOrder);

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// NakBsplineBoundaryCombigrid Basis object
  SNakBsplineBoundaryBase base;
};

}  // namespace base
}  // namespace sgpp
