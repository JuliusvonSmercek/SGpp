// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/sle/CloneableSLE.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/solver/sle/external/Gmmpp.hpp>

#include <sgpp/globaldef.hpp>

#ifdef USE_GMMPP
#include <gmm/gmm.h>
#endif /* USE_GMMPP */

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace sgpp {
namespace solver {

#ifdef USE_GMMPP
/**
 * Gmm++ status printing callback.
 *
 * @param iter  iteration information
 */
void callback(const gmm::iteration& iter) {
  base::Printer::getInstance().printStatusUpdate(
      "solving with Gmm++ (k = " + std::to_string(iter.get_iteration()) +
      ", residual norm = " + std::to_string(iter.get_res()) + ")");
}

/**
 * @param       A   coefficient matrix
 * @param       b   right-hand side
 * @param[out]  x   solution of the linear system
 * @return          whether all went well (false if errors occurred)
 */
bool solveInternal(const gmm::csr_matrix<double>& A, base::DataVector& b, base::DataVector& x) {
  // allow warnings
  gmm::warning_level::level(1);
  const size_t n = b.getSize();
  std::vector<double> bVec(b.getPointer(), b.getPointer() + n);
  std::vector<double> xVec(n, 0.0);

  // ILU preconditioning
  base::Printer::getInstance().printStatusUpdate("constructing preconditioner");
  gmm::ilu_precond<gmm::csr_matrix<double>> P(A);

  base::Printer::getInstance().printStatusNewLine();
  base::Printer::getInstance().printStatusUpdate("solving with Gmm++");

  gmm::iteration iter(1e-6, 0, 10000);
  iter.set_callback(&callback);

  try {
    // call GMRES
    gmm::gmres(A, xVec, bVec, P, 50, iter);

    double res = iter.get_res();

    if (iter.converged() && (res < 1e3)) {
      // GMRES converged
      x = base::DataVector(xVec);
      base::Printer::getInstance().printStatusUpdate(
          "solving with Gmm++ (k = " + std::to_string(iter.get_iteration()) +
          ", residual norm = " + std::to_string(res) + ")");
      base::Printer::getInstance().printStatusEnd();
      return true;
    } else {
      // GMRES didn't converge ==> try again without preconditioner
      gmm::identity_matrix P;

      base::Printer::getInstance().printStatusNewLine();
      base::Printer::getInstance().printStatusUpdate(
          "solving with preconditioner failed, trying again without one");
      base::Printer::getInstance().printStatusNewLine();

      // call GMRES again
      gmm::gmres(A, xVec, bVec, P, 50, iter);
      res = iter.get_res();

      if (iter.converged() && (res < 1e3)) {
        x = base::DataVector(xVec);
        base::Printer::getInstance().printStatusUpdate(
            "solving with Gmm++ (k = " + std::to_string(iter.get_iteration()) +
            ", residual norm = " + std::to_string(res) + ")");
        base::Printer::getInstance().printStatusEnd();
        return true;
      } else {
        base::Printer::getInstance().printStatusEnd(
            "error: Could not solve linear system, "
            "method didn't converge");
        return false;
      }
    }
  } catch (std::exception& e) {
    base::Printer::getInstance().printStatusEnd("error: Could not solve linear system, what(): " +
                                                std::string(e.what()));
    return false;
  }
}
#endif /* USE_GMMPP */

Gmmpp::~Gmmpp() {}

bool Gmmpp::solve(base::SLE& system, base::DataVector& b, base::DataVector& x) const {
#ifdef USE_GMMPP
  base::Printer::getInstance().printStatusBegin("Solving linear system (Gmm++)...");

  const size_t n = system.getDimension();
  size_t nnz = 0;
  size_t rowsDone = 0;
  gmm::csr_matrix<double> A2;

  {
    gmm::row_matrix<gmm::rsvector<double>> A(n, n);

// parallelize only if the system is cloneable
#pragma omp parallel if (system.isCloneable()) shared(system, A, nnz, rowsDone, n) default(none)
    {
      base::SLE* system2 = &system;
#ifdef _OPENMP
      std::unique_ptr<base::CloneableSLE> clonedSLE;

      if (system.isCloneable() && (omp_get_max_threads() > 1)) {
        dynamic_cast<base::CloneableSLE&>(system).clone(clonedSLE);
        system2 = clonedSLE.get();
      }

#endif /* _OPENMP */

// copy system matrix to Gmm++ matrix object
#pragma omp for ordered schedule(static)

      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
          double entry = system2->getMatrixEntry(i, j);

          if (entry != 0) {
#pragma omp critical
            {
              A(i, j) = entry;
              nnz++;
            }
          }
        }

#pragma omp atomic
        rowsDone++;

        // status message
        if (rowsDone % 100 == 0) {
          char str[10];
          snprintf(str, sizeof(str), "%.1f%%",
                   static_cast<double>(rowsDone) / static_cast<double>(n) * 100.0);
          base::Printer::getInstance().printStatusUpdate("constructing sparse matrix (" +
                                                         std::string(str) + ")");
        }
      }
    }

    // the Gmm++ manual said to do so... no idea if that's necessary
    gmm::clean(A, 1e-12);
    gmm::copy(A, A2);
  }

  base::Printer::getInstance().printStatusUpdate("constructing sparse matrix (100.0%)");
  base::Printer::getInstance().printStatusNewLine();

  // print ratio of nonzero entries
  {
    char str[10];
    double nnz_ratio = static_cast<double>(nnz) / (static_cast<double>(n) * static_cast<double>(n));
    snprintf(str, sizeof(str), "%.1f%%", nnz_ratio * 100.0);
    base::Printer::getInstance().printStatusUpdate("nnz ratio: " + std::string(str));
    base::Printer::getInstance().printStatusNewLine();
  }

  x.resize(n);
  bool result = solveInternal(A2, b, x);
  return result;
#else
  std::cerr << "Error in solver::Gmmpp::solve: "
            << "SG++ was compiled without Gmm++ support!\n";
  return false;
#endif /* USE_GMMPP */
}
}  // namespace solver
}  // namespace sgpp
