// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/globaldef.hpp>
#include "common/basis/NakPBsplineBasis.hpp"

namespace sgpp {
namespace base {

/**
 * Operation for evaluating nak polynomial B-spline linear combinations on Noboundary grids.
 */
class OperationEvalNakPBsplineNaive : public OperationEval {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationEvalNakPBsplineNaive(GridStorage& storage, size_t degree)
      : storage(storage), base(degree), pointInUnitCube(storage.getDimension()) {}

  /**
   * Destructor.
   */
  ~OperationEvalNakPBsplineNaive() override {}

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @return          value of the linear combination
   */
  double eval(const DataVector& alpha, const DataVector& point) override;

  /**
   * @param      alpha  coefficient matrix (each column is a coefficient vector)
   * @param      point  evaluation point
   * @param[out] value  values of linear combination
   */
  void eval(const DataMatrix& alpha, const DataVector& point, DataVector& value) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SNakPBsplineBase base;
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
};

}  // namespace base
}  // namespace sgpp
