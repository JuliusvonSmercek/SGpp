// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALHESSIANMODFUNDAMENTALSPLINE_HPP
#define OPERATIONEVALHESSIANMODFUNDAMENTALSPLINE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessian.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace sgpp {
namespace base {

/**
 * Operation for evaluating modified B-spline linear combinations on Noboundary grids,
 * their gradients and their Hessians.
 */
class OperationEvalHessianModFundamentalSplineNaive : public
  OperationEvalHessian {
 public:
  /**
   * Constructor.
   *
   * @param storage   storage of the sparse grid
   * @param degree    B-spline degree
   */
  OperationEvalHessianModFundamentalSplineNaive(GridStorage& storage, size_t degree) :
    storage(storage),
    base(degree),
    pointInUnitCube(storage.getDimension()),
    innerDerivative(storage.getDimension()) {
  }

  /**
   * Destructor.
   */
  ~OperationEvalHessianModFundamentalSplineNaive() override {
  }

  /**
   * @param       alpha     coefficient vector
   * @param       point     evaluation point
   * @param[out]  gradient  gradient vector of the linear combination
   * @param[out]  hessian   Hessian matrix of the linear combination
   * @return                value of the linear combination
   */
  double evalHessian(const DataVector& alpha,
                      const DataVector& point,
                      DataVector& gradient,
                      DataMatrix& hessian) override;

  /**
   * @param       alpha     coefficient matrix (each column is a coefficient vector)
   * @param       point     evaluation point
   * @param[out]  value     values of the linear combination
   * @param[out]  gradient  Jacobian of the linear combination (each row is a gradient vector)
   * @param[out]  hessian   vector of Hessians of the linear combination
   */
  void evalHessian(const DataMatrix& alpha,
                   const DataVector& point,
                   DataVector& value,
                   DataMatrix& gradient,
                   std::vector<DataMatrix>& hessian) override;

 protected:
  /// storage of the sparse grid
  GridStorage& storage;
  /// 1D B-spline basis
  SFundamentalSplineModifiedBase base;
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
  /// inner derivative (temporary vector)
  DataVector innerDerivative;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALHESSIANMODFUNDAMENTALSPLINE_HPP */
