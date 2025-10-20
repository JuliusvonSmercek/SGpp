// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <boost/math/distributions/normal.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace sgpp {
namespace optimization {
/**
 * @brief Helper classes and functions for the IterativeGridGeneratorFullAdaptiveRitterNovak.
 */
namespace IGGFARNHelper {

/**
 * @brief Stores and normalizes the hyperparameters for the refinement criterion.
 */
class Hyperparameters {
 private:
  double levelDegree;
  double threshold;
  double extendedThreshold;
  double gradient;
  double secondGradient;
  double gaussianProcess;

 public:
  /**
   * @brief Constructor that takes individual weights for each criterion.
   * The weights are normalized to sum to 1.
   *
   * @param levelDegree       Weight for the level-based criterion.
   * @param threshold         Weight for the threshold criterion (value at parent).
   * @param extendedThreshold Weight for the extended threshold criterion (value at children).
   * @param gradient          Weight for the gradient-based criterion.
   * @param secondGradient    Weight for the second-order gradient criterion.
   * @param gaussianProcess   Weight for the Gaussian Process uncertainty criterion.
   */
  Hyperparameters(double levelDegree, double threshold, double extendedThreshold, double gradient,
                  double secondGradient, double gaussianProcess);

  /**
   * @brief Default constructor, initializes all weights to zero.
   */
  Hyperparameters();

  double getLevelDegree() const;
  double getThreshold() const;
  double getExtendedThreshold() const;
  double getGradient() const;
  double getSecondGradient() const;
  double getGaussianProcess() const;
};

/// Stream operator for printing hyperparameters.
std::ostream& operator<<(std::ostream& os, const Hyperparameters& hyperparameters);

/**
 * @brief A class to represent a value with uncertainty bounds, typically from sampling.
 */
class UncertaintyBound {
 private:
  double _lower;
  double _value;
  double _upper;
  double _coV;  // Coefficient of Variation

 public:
  /**
   * @brief Construct from explicit lower, value, upper, and CoV.
   *
   * @param lower Lower bound.
   * @param value The estimated value.
   * @param upper Upper bound.
   * @param coV   Coefficient of variation.
   */
  UncertaintyBound(double lower, double value, double upper, double coV);

  /**
   * @brief Default constructor, initializes to zero.
   */
  UncertaintyBound();

  /**
   * @brief Construct from a number of positive samples out of a total number of samples.
   * Calculates the value and confidence interval.
   *
   * @param samples Number of positive samples.
   * @param N       Total number of samples.
   */
  UncertaintyBound(size_t samples, size_t N);

  /**
   * @brief Divides this uncertainty bound by another. Used for relative error estimation.
   *
   * @param other The divisor.
   * @return A new UncertaintyBound representing the quotient.
   */
  UncertaintyBound div(const UncertaintyBound& other) const;

  double lower() const;
  double value() const;
  double upper() const;
  double coV() const;
};

/// Stream operator for printing uncertainty bounds.
std::ostream& operator<<(std::ostream& os, const UncertaintyBound& bounds);

/**
 * @brief Implements a stopping criterion based on the convergence of the surrogate model.
 *
 * It compares the current grid's surrogate with surrogates from previous iterations. If the
 * relative mismatch probability between them is consistently below a tolerance, it signals to stop.
 */
class StoppingCriterion {
 protected:
  static const size_t stopping_criterion_count_difference = 50;
  static constexpr double stopping_criterion_tolerance = 0.025;

  const std::vector<sgpp::base::DataVector>& monteCarloPoints;
  std::vector<
      std::tuple<size_t, std::vector<double>, std::unique_ptr<base::Grid>, base::DataVector>>
      last_grid_values;

  /**
   * @brief Measures the failure and mismatch probabilities between two functions.
   *
   * @param points             Points to sample from.
   * @param surrogate_values   Cache for values of the first function.
   * @param original_values    Cache for values of the second function.
   * @param surrogate_f        First function (e.g., current surrogate).
   * @param original_f         Second function (e.g., previous surrogate).
   * @param early_stop_criterion A predicate to decide if sampling can stop early.
   * @return A tuple containing (number of samples used, failure probability, mismatch
   * probability).
   */
  std::tuple<size_t, UncertaintyBound, UncertaintyBound> measure_error(
      const std::vector<base::DataVector>& points, std::vector<double>& surrogate_values,
      std::vector<double>& original_values,
      const std::function<double(const sgpp::base::DataVector&)>& surrogate_f,
      const std::function<double(const sgpp::base::DataVector&)>& original_f,
      const std::function<bool(const UncertaintyBound&, const UncertaintyBound&, size_t)>&
          early_stop_criterion);

 public:
  /**
   * @brief Constructor.
   *
   * @param monteCarloPoints Points to use for validation.
   */
  explicit StoppingCriterion(const std::vector<sgpp::base::DataVector>& monteCarloPoints);

  /**
   * @brief Destructor.
   */
  ~StoppingCriterion();

  /**
   * @brief Checks if the grid generation should stop.
   *
   * @param grid  The current grid.
   * @param alpha The current surplus vector.
   * @return True if the stopping criterion is met, false otherwise.
   */
  bool shouldStop(base::Grid& grid, const base::DataVector& alpha);
};

// --- Implementations of helper classes ---

namespace {
// Confidence level for Wilson score interval
constexpr double accuracy = 0.999;

/**
 * @brief Computes the quantile of the standard normal distribution (inverse of the CDF).
 * @param p Probability (must be in (0, 1)).
 * @return The corresponding z-score.
 */
static double inverse_phi(double p) {
  boost::math::normal_distribution<> normal_dist(0.0, 1.0);
  return boost::math::quantile(normal_dist, p);
}

/**
 * @brief Calculates the Wilson score confidence interval for a binomial proportion.
 * It is more reliable than the Wald interval, especially for p near 0 or 1.
 *
 * @param p     Observed proportion.
 * @param N     Total number of samples.
 * @param gamma Confidence level (e.g., 0.999).
 * @return A pair containing the lower and upper bounds of the confidence interval.
 */
static std::pair<double, double> confidence_interval_wilson(double p, size_t N, double gamma) {
  assert(0.0 <= p && p <= 1.0);
  assert(0.0 <= gamma && gamma <= 1.0);
  assert(0 < N);

  const double alpha = 1.0 - gamma;
  const double c = inverse_phi(1.0 - alpha / 2.0);
  const double cSquared = c * c;

  const double denominator = 2.0 * (static_cast<double>(N) + cSquared);
  const double numerator_center = 2.0 * static_cast<double>(N) * p + cSquared;
  const double numerator_margin =
      c * std::sqrt(cSquared + 4.0 * static_cast<double>(N) * p * (1.0 - p));

  const double lower = (numerator_center - numerator_margin) / denominator;
  const double upper = (numerator_center + numerator_margin) / denominator;

  return std::make_pair(std::max(0.0, lower), std::min(1.0, upper));
}
}  // namespace

inline UncertaintyBound::UncertaintyBound(double lower, double value, double upper, double coV)
    : _lower(lower), _value(value), _upper(upper), _coV(coV) {}

inline UncertaintyBound::UncertaintyBound() : UncertaintyBound(0.0, 0.0, 0.0, 0.0) {}

inline UncertaintyBound::UncertaintyBound(const size_t samples, const size_t N)
    : _value(static_cast<double>(samples) / static_cast<double>(N)) {
  assert(samples <= N);
  assert(0 < N);

  // For a Bernoulli distribution, the coefficient of variation is sqrt((1-p)/(p*N)).
  // If p is 0, CoV is infinite. We return 1.0 as a placeholder.
  this->_coV =
      (0.0 == this->_value) ? 1.0 : std::sqrt((1.0 - _value) / (_value * static_cast<double>(N)));
  std::tie(this->_lower, this->_upper) = confidence_interval_wilson(this->_value, N, accuracy);
}

inline UncertaintyBound UncertaintyBound::div(const UncertaintyBound& other) const {
  // Division by zero is handled by returning a safe but potentially wide interval.
  const double lower = (0.0 == other._upper) ? 0.0 : this->_lower / other._upper;
  const double value = (0.0 == other._value) ? 1.0 : this->_value / other._value;
  const double upper = (0.0 == other._lower) ? 1.0 : this->_upper / other._lower;

  return UncertaintyBound(lower, value, upper, -1.0);  // CoV not meaningful after division
}

inline double UncertaintyBound::lower() const { return _lower; }
inline double UncertaintyBound::value() const { return _value; }
inline double UncertaintyBound::upper() const { return _upper; }
inline double UncertaintyBound::coV() const { return _coV; }

inline std::ostream& operator<<(std::ostream& os, const UncertaintyBound& bounds) {
  os << "[" << 100.0 * bounds.lower() << "%, " << 100.0 * bounds.value() << "%, "
     << 100.0 * bounds.upper() << "%]";
  return os;
}

}  // namespace IGGFARNHelper

/**
 * @brief Iterative grid generation based on a modified, fully adaptive Ritter/Novak's refinement
 * criterion.
 *
 * This generator uses a combination of weighted criteria to decide which grid point to refine.
 * The criteria include the grid point's level, function value (threshold), gradient, and more.
 * A large number of hyperparameters allow these criteria to be weighted differently.
 *
 * @warning This class uses HashRefinementMultiple, so it generates grids that may not meet the
 * "hierarchical ancestors" requirement.
 *
 * Literature: Julius von Smercek, Bachelor Thesis ``Critical Event Estimation using Gradient-based
 * Refinement Methods for Adaptive Sparse Grids'', 2025
 *
 * @see HashRefinementMultiple
 */
class IterativeGridGeneratorFullAdaptiveRitterNovak : public IterativeGridGenerator {
 public:
  /// Default level of the initial regular sparse grid.
  static const base::level_t DEFAULT_INITIAL_LEVEL = 3;
  /// Default maximal level of grid points.
  static const base::level_t DEFAULT_MAX_LEVEL = 20;

  /**
   * Constructor.
   * Do not destruct the grid object before this generator object!
   *
   * @param f                 Objective function.
   * @param grid              Grid object (should be empty).
   * @param N                 Maximal number of grid points.
   * @param hyperparameters   Hyperparameters for the adaptive weight function.
   * @param refineLeftOrRight Whether to only refine the left or right child in each dimension,
   * activating the "full adaptive" approach.
   * @param initialLevel      Level of the initial regular sparse grid.
   * @param maxLevel          Maximal level permitted for any grid point.
   */
  IterativeGridGeneratorFullAdaptiveRitterNovak(
      base::ScalarFunction& f, base::Grid& grid, size_t N,
      const IGGFARNHelper::Hyperparameters& hyperparameters, bool refineLeftOrRight = true,
      const base::level_t initialLevel = DEFAULT_INITIAL_LEVEL,
      const base::level_t maxLevel = DEFAULT_MAX_LEVEL);

  /**
   * Destructor.
   */
  ~IterativeGridGeneratorFullAdaptiveRitterNovak() override;

  /**
   * Sets the functions for interacting with a Gaussian Process model.
   * This is optional and only needed if the `gaussianProcess` hyperparameter is non-zero.
   *
   * @param init        Function to initialize the Gaussian process with the problem's dimension.
   * @param addPattern  Function to add a new training point (coordinate and value) to the GP.
   * @param evaluate    Function to evaluate the GP at potential new grid points, returning
   * mean and variance.
   */
  void setGaussianProcess(
      const std::function<void(const size_t dimension)>& init,
      const std::function<void(const base::DataVector& x, double y)>& addPattern,
      const std::function<std::vector<std::vector<std::array<std::pair<double, double>, 2>>>(
          const std::vector<std::vector<std::array<base::GridPoint, 2>>>& refineableGridPoints)>&
          evaluate);

  /**
   * Generates the grid adaptively until the maximum number of points `N` is reached
   * or a stopping criterion is met.
   *
   * @return True on success, false otherwise.
   */
  bool generate() override;

  /**
   * @return The currently used hyperparameters.
   */
  const IGGFARNHelper::Hyperparameters& getHyperparameters() const;

  /**
   * @param hyperparameters New hyperparameters for the adaptive weight function.
   */
  void setHyperparameters(const IGGFARNHelper::Hyperparameters& hyperparameters);

  /**
   * @return Level of the initial regular sparse grid.
   */
  base::level_t getInitialLevel() const;

  /**
   * @param initialLevel New level for the initial regular sparse grid.
   */
  void setInitialLevel(base::level_t initialLevel);

  /**
   * @return Maximal permitted level for grid points.
   */
  base::level_t getMaxLevel() const;

  /**
   * @param maxLevel New maximal level for grid points.
   */
  void setMaxLevel(base::level_t maxLevel);

  /**
   * Activates the stopping criterion.
   * The generation will stop if the surrogate model's predictions stabilize between iterations.
   *
   * @param validationPoints A set of points used to compare surrogate models.
   */
  void activateStoppingCriterion(const std::vector<sgpp::base::DataVector>& validationPoints);

  /**
   * Disables the stopping criterion.
   */
  void disableStoppingCriterion();

  /**
   * @return True if the stopping criterion is enabled, false otherwise.
   */
  bool isStoppingCriterionEnabled() const;

  /**
   * @return The coefficients (surpluses) of the final sparse grid interpolant.
   */
  const base::DataVector& getAlpha() const;

 protected:
  /// Hyperparameters for the adaptive weight function.
  IGGFARNHelper::Hyperparameters hyperparameters;
  /// Level of the initial regular sparse grid.
  base::level_t initialLevel;
  /// Maximal level of grid points.
  base::level_t maxLevel;
  /// Validation points for the stopping criterion. Empty if disabled.
  std::vector<sgpp::base::DataVector> stoppingCriterionValidationPoints;
  /// Whether to refine only the left or right child in each dimension.
  bool refineLeftOrRight;
  /// Coefficients of the sparse grid surrogate.
  base::DataVector alpha;
  /// Function to initialize the Gaussian process.
  std::function<void(const size_t dimension)> gpInit;
  /// Function to add a pattern to the Gaussian process.
  std::function<void(const base::DataVector& x, double y)> gpAddPattern;
  /// Function to evaluate the Gaussian process at given points.
  std::function<std::vector<std::vector<std::array<std::pair<double, double>, 2>>>(
      const std::vector<std::vector<std::array<base::GridPoint, 2>>>& refineableGridPoints)>
      gpEvaluate;
};
}  // namespace optimization
}  // namespace sgpp
