// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTestPrewavelet.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

double OperationTestPrewavelet::test(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                                      sgpp::base::DataVector& classes) {
  base::PrewaveletBasis<unsigned int, unsigned int> base;
  return test_dataset(this->storage, base, alpha, data, classes);
}

double OperationTestPrewavelet::testMSE(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataMatrix& data,
                                         sgpp::base::DataVector& refValues) {
  base::PrewaveletBasis<unsigned int, unsigned int> base;
  return test_dataset_mse(this->storage, base, alpha, data, refValues);
}

double OperationTestPrewavelet::testWithCharacteristicNumber(
    sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes,
    sgpp::base::DataVector& charaNumbers) {
  base::PrewaveletBasis<unsigned int, unsigned int> base;
  return test_datasetWithCharacteristicNumber(this->storage, base, alpha, data, classes,
                                              charaNumbers, 0.0);
}

void OperationTestPrewavelet::calculateROCcurve(sgpp::base::DataVector& alpha,
                                                sgpp::base::DataMatrix& data,
                                                sgpp::base::DataVector& classes,
                                                sgpp::base::DataVector& thresholds,
                                                sgpp::base::DataMatrix& ROC_curve) {
  base::PrewaveletBasis<unsigned int, unsigned int> base;
  test_calculateROCcurve(this->storage, base, alpha, data, classes, thresholds, ROC_curve);
}

}  // namespace datadriven
}  // namespace sgpp
