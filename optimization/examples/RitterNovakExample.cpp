// #define GAUSSIAN_PROCESS_AVAILABLE
#ifdef GAUSSIAN_PROCESS_AVAILABLE
#include <gp.h>
#include <rprop.h>
#endif

#include <cassert>
#include <cmath>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorFullAdaptiveRitterNovak.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp>
#include <vector>

void gridGenerationRitterNovak() {
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;
  const double adaptivity = 0.15;
  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(*function, *grid, maxGridPoints, adaptivity);

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  sgpp::base::HierarchisationSLE hierSLE(*grid);
  const size_t N = hierSLE.getDimension();

  sgpp::base::sle_solver::Auto sleSolver;
  sgpp::base::DataVector alpha(N);
  sgpp::base::DataVector b = gridGen.getFunctionValues();

  if (!sleSolver.solve(hierSLE, b, alpha)) {
    throw std::runtime_error("Solving failed");
  }

  // if success print all point on the grid with their coefficients
  for (size_t i = 0; i < gridStorage.getSize(); ++i) {
    sgpp::base::DataVector point(dimension);
    gridStorage[i].getStandardCoordinates(point);
    std::cout << "Point " << i << ": " << point.toString() << ", Coefficient: " << alpha[i] << std::endl;
  }
}

void gridGenerationFullAdaptiveRitterNovakMinimal() {
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;
  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  std::vector<std::pair<double, double>> domain(dimension);
  for (size_t t = 0; t < dimension; ++t) {
    domain[t] = std::make_pair(0.0, 1.0);
  }

  const sgpp::optimization::IGGFARNHelper::Hyperparameters hyperparameters(0.49825808, 0.00042716, 0.00003296,
                                                                           0.14564321, 0.11427497, 0);

  sgpp::optimization::IterativeGridGeneratorFullAdaptiveRitterNovak gridGen(*function, *grid, maxGridPoints,
                                                                            hyperparameters, domain);

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  const sgpp::base::DataVector alpha = gridGen.getAlpha();

  // if success print all point on the grid with their coefficients
  for (size_t i = 0; i < gridStorage.getSize(); ++i) {
    sgpp::base::DataVector point(dimension);
    gridStorage[i].getStandardCoordinates(point);
    std::cout << "Point " << i << ": " << point.toString() << ", Coefficient: " << alpha[i] << std::endl;
  }
}

#ifdef GAUSSIAN_PROCESS_AVAILABLE
void gridGenerationFullAdaptiveRitterNovak() {
  const size_t dimension = 2;
  const size_t bSplineDegree = 3;
  const size_t maxGridPoints = 50;
  const size_t stoppingCriterionMonteCarloPointCount = 50'000;
  const bool useStoppingCriterion = true;
  std::unique_ptr<sgpp::base::ScalarFunction> function =
      std::make_unique<sgpp::optimization::test_problems::RastriginObjective>(dimension);

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createBsplineGrid(dimension, bSplineDegree));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  std::vector<std::pair<double, double>> domain(dimension);
  for (size_t t = 0; t < dimension; ++t) {
    domain[t] = std::make_pair(0.0, 1.0);
  }

  const sgpp::optimization::IGGFARNHelper::Hyperparameters hyperparameters(0.49825808, 0.00042716, 0.00003296,
                                                                           0.14564321, 0.11427497, 0.24136362);

  sgpp::optimization::IterativeGridGeneratorFullAdaptiveRitterNovak gridGen(*function, *grid, maxGridPoints,
                                                                            hyperparameters, domain);

  if (useStoppingCriterion) {
    std::vector<sgpp::base::DataVector> validationPoints;
    for (size_t i = 0; i < stoppingCriterionMonteCarloPointCount; ++i) {
      sgpp::base::DataVector point(dimension);
      for (size_t t = 0; t < dimension; ++t) {
        const double rand01 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        point[t] = domain[t].first + rand01 * (domain[t].second - domain[t].first);
      }
      validationPoints.push_back(point);
    }
    gridGen.activateStoppingCriterion(validationPoints);
  }

  std::unique_ptr<libgp::GaussianProcess> gp;
  gridGen.setGaussianProcess(
      [&](const size_t dim) {
        assert(dim == dimension);
        // initialize Gaussian process for 2-D input using the squared exponential
        // covariance function with additive white noise.
        // Define the covariance function (kernel)
        // libgp uses a string format: CovSum( CovMatern3iso, CovNoise )
        // CovMatern3iso: Matern kernel with smoothness 3/2
        // CovNoise: White noise kernel
        gp = std::make_unique<libgp::GaussianProcess>(dimension, "CovSum( CovMatern3iso, CovNoise)");

        // How quickly the function changes
        const double length_scale = 0.05;
        // Overall variance of the function
        const double signal_variance = 1.0;
        // Assumed noise in the observations
        const double noise_variance = 1e-3;

        Eigen::VectorXd params(gp->covf().get_param_dim());
        params << std::log(length_scale), std::log(signal_variance), std::log(noise_variance);

        // These are initial guesses - ideally, they should be optimized!
        gp->covf().set_loghyper(params);

        std::cout << "Using Kernel: " << gp->covf().to_string() << std::endl;
        std::cout << "Initial hyperparameters: [" << length_scale << ", " << signal_variance << ", " << noise_variance
                  << "]" << std::endl;
      },
      [&](const sgpp::base::DataVector& x, const double y) { gp->add_pattern(x.data(), y); },
      [&](const std::vector<std::vector<std::array<sgpp::base::GridPoint, 2>>>& refineableGridPoints,
          const std::vector<std::vector<std::array<bool, 2>>>& ignore)
          -> std::vector<std::vector<std::array<std::pair<double, double>, 2>>> {
        libgp::RProp rProp;
        rProp.init();
        rProp.maximize(&*gp, 5, false);

        Eigen::VectorXd params{gp->covf().get_loghyper()};
        params[2] = std::log(1e-3);
        gp->covf().set_loghyper(params);

        // // TODO: gp_params[2] should always be 1e-3???
        // const Eigen::VectorXd params2 = gp->covf().get_loghyper();
        // std::cout << gp->get_sampleset_size() << ": " << gp->log_likelihood() << " - [" << std::exp(params2[0])
        // <<
        // ", "
        //           << std::exp(params2[1]) << ", " << std::exp(params2[2]) << "]" << std::endl;

        std::vector<std::vector<std::array<std::pair<double, double>, 2>>> result(refineableGridPoints.size());
        for (size_t t = 0; t < refineableGridPoints.size(); ++t) {
          result[t].resize(refineableGridPoints[t].size());
          for (size_t i = 0; i < refineableGridPoints[t].size(); ++i) {
            result[t][i] = std::array<std::pair<double, double>, 2>();
            for (size_t lr = 0; lr < 2; ++lr) {
              if (ignore[i][t][lr]) {
                continue;
              }

              sgpp::base::DataVector new_point;
              refineableGridPoints[t][i][lr].getStandardCoordinates(new_point);

              const double predicted_mean = gp->f(new_point.data());
              // Warning: somehow predicted_variance can be small negative -> numerical instabilities?
              const double predicted_variance = std::abs(gp->var(new_point.data()));
              const double predicted_stddev = std::sqrt(predicted_variance);

              result[t][i][lr] = std::make_pair(predicted_mean, predicted_stddev);
            }
          }
        }

        return result;
      });

  if (!gridGen.generate()) {
    throw std::runtime_error("Grid generation failed, exiting");
  }

  const sgpp::base::DataVector alpha = gridGen.getAlpha();

  // if success print all point on the grid with their coefficients
  for (size_t i = 0; i < gridStorage.getSize(); ++i) {
    sgpp::base::DataVector point(dimension);
    gridStorage[i].getStandardCoordinates(point);
    std::cout << "Point " << i << ": " << point.toString() << ", Coefficient: " << alpha[i] << std::endl;
  }
}
#endif

int main(int argc, const char* argv[]) {
  gridGenerationRitterNovak();
  gridGenerationFullAdaptiveRitterNovakMinimal();
#ifdef GAUSSIAN_PROCESS_AVAILABLE
  gridGenerationFullAdaptiveRitterNovak();
#endif
  return 0;
}
