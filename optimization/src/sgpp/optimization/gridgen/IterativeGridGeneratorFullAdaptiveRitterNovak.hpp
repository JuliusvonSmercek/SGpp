// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <cstddef>

// TODO: delete, when compiles
// TODO: doc when compiles -> thesis abschnitt gem
// #include <cassert>
// #include <iomanip>
// #include <iostream>
// #include <ostream>
// #include <sgpp_base.hpp>
// #include <sgpp_optimization.hpp>

namespace sgpp {
namespace optimization {

/**
 * Iterative grid generation based on a modified, full adaptive Ritter/Novak's refinement criterion.
 * A large number of hyperparameters allow various criteria to be weighted differently. Caution:
 * This class uses HashRefinementMultiple, so it generates grids that don't meet the "hierarchical
 * ancestors" requirement!
 *
 * Literature: Julius von Smercek, Bachelor Thesis ``Critical Event Estimation using Gradient-based
 * Refinement Methods for Adaptive Sparse Grids'', 2025
 *
 * @see HashRefinementMultiple
 */
class IterativeGridGeneratorFullAdaptiveRitterNovak : public IterativeGridGenerator {
 public:
  /// default level of initial regular sparse grid
  static const base::level_t DEFAULT_INITIAL_LEVEL = 3;
  /// default maximal level of grid points
  static const base::level_t DEFAULT_MAX_LEVEL = 20;

  class Hyperparameters {
   private:
    double levelDegree, threshold, extendedThreshold, gradient, secondGradient, gaussianProcess;

   public:
    Hyperparameters(const double levelDegree, const double threshold,
                    const double extendedThreshold, const double gradient,
                    const double secondGradient, const double gaussianProcess);

    Hyperparameters();

    double getLevelDegree() const;

    double getThreshold() const;

    double getExtendedThreshold() const;

    double getGradient() const;

    double getSecondGradient() const;

    double getGaussianProcess() const;
  };

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f                     objective function
   * @param grid                  grid (should be empty)
   * @param N                     maximal number of grid points
   * @param hyperparameters       hyperparameters for adaptive weight function
   * @param useStoppingCriterion  whether to use the stopping criterion
   * @param refineLeftOrRight     whether to only refine left or right child in each dimension,
   * activates the full adaptive approach
   * @param initialLevel          level of initial regular sparse grid
   * @param maxLevel              maximal level of grid points
   */
  IterativeGridGeneratorFullAdaptiveRitterNovak(
      base::ScalarFunction& f, base::Grid& grid, const size_t N,
      const Hyperparameters hyperparameters, const bool useStoppingCriterion,
      const bool refineLeftOrRight = true, const base::level_t initialLevel = DEFAULT_INITIAL_LEVEL,
      const base::level_t maxLevel = DEFAULT_MAX_LEVEL);

  /**
   * Set the Gaussian process functions.
   *
   * @param addPattern         function to add patterns to Gaussian process
   * @param evaluate           function to evaluate the Gaussian process at given points
   */
  void setGaussianProcess(
      const std::function<void(const base::DataVector& x, const double y)>& addPattern,
      const std::function<std::vector<std::vector<std::array<std::pair<double, double>, 2>>>(
          const std::vector<std::vector<std::array<base::GridPoint, 2>>>& refineableGridPoints)>&
          evaluate) {}

  /**
   * Destructor.
   */
  ~IterativeGridGeneratorFullAdaptiveRitterNovak() override;

  /**
   * Generate the grid.
   *
   * @return true on success, otherwise false
   */
  bool generate() override;

  /**
   * @return used hyperparameters
   */
  Hyperparameters getHyperparameters() const;

  /**
   * @param hyperparameters   hyperparameters for adaptive weight function
   */
  void setHyperparameters(const Hyperparameters hyperparameters);

  /**
   * @return              level of initial regular sparse grid
   */
  base::level_t getInitialLevel() const;

  /**
   * @param initialLevel  level of initial regular sparse grid
   */
  void setInitialLevel(base::level_t initialLevel);

  /**
   * @return          maximal level of grid points
   */
  base::level_t getMaxLevel() const;

  /**
   * @param maxLevel  maximal level of grid points
   */
  void setMaxLevel(base::level_t maxLevel);

  /**
   * @return whether to use the stopping criterion
   */
  bool getUseStoppingCriterion() const;

  /**
   * @param useStoppingCriterion  whether to use the stopping criterion
   */
  void setUseStoppingCriterion(const bool useStoppingCriterion);

  /**
   * @return coefficients of the sparse grid surrogate
   */
  base::DataVector getAlpha() const;

 protected:
  /// hyperparameters for adaptive weight function based on ritter-novak
  Hyperparameters hyperparameters;
  /// level of initial regular sparse grid
  base::level_t initialLevel;
  /// maximal level of grid points
  base::level_t maxLevel;
  /// whether to use the stopping criterion
  bool useStoppingCriterion;
  /// whether to only refine left or right child in each dimension
  bool refineLeftOrRight;
  /// coefficients of the sparse grid surrogate
  base::DataVector alpha;
  /// function to add patterns to Gaussian process
  std::function<void(const base::DataVector& x, const double y)> addPattern;
  /// function to evaluate the Gaussian process at given points
  std::function<std::vector<std::vector<std::array<std::pair<double, double>, 2>>>(
      const std::vector<std::vector<std::array<base::GridPoint, 2>>>& refineableGridPoints)>
      evaluate;

  // void setDebugFunction(
  //     const std::function<void(base::Grid& grid, const base::DataVector& functionValues)>&
  //         debugFunction) {
  //   this->debugFunction = debugFunction;
  // }

  // void setDebugFunction(
  //     std::function<void(base::Grid& grid, const base::DataVector& functionValues)>&&
  //         debugFunction) {
  //   this->debugFunction = std::move(debugFunction);
  // }

  // void setDebugStoppingCriterionFunction(
  //     const std::function<void(const size_t, const size_t, const UncertaintyBound, const size_t,
  //                              const PointCollection&, std::vector<double>&,
  //                              base::InterpolantScalarFunction&)>&
  //                              debugStoppingCriterionFunction) {
  //   this->debugStoppingCriterionFunction = debugStoppingCriterionFunction;
  // }

  // void setDebugStoppingCriterionFunction(
  //     std::function<void(const size_t, const size_t, const UncertaintyBound, const size_t,
  //                        const PointCollection&, std::vector<double>&,
  //                        base::InterpolantScalarFunction&)>&& debugStoppingCriterionFunction) {
  //   this->debugStoppingCriterionFunction = std::move(debugStoppingCriterionFunction);
  // }

  // std::function<void(base::Grid& grid, const base::DataVector& functionValues)> debugFunction =
  //     [](base::Grid& grid, const base::DataVector& functionValues) {};

  // std::function<void(const size_t, const size_t, const UncertaintyBound, const size_t,
  //                    const PointCollection&, std::vector<double>&,
  //                    base::InterpolantScalarFunction&)>
  //     debugStoppingCriterionFunction =
  //         [](const size_t, const size_t, const UncertaintyBound, const size_t,
  //            const PointCollection&, std::vector<double>&, base::InterpolantScalarFunction&) {};
};

std::ostream& operator<<(
    std::ostream& os,
    const IterativeGridGeneratorFullAdaptiveRitterNovak::Hyperparameters& hyperparameters);
}  // namespace optimization
}  // namespace sgpp
