// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/sle/CloneableSLE.hpp>
#include <sgpp/base/tools/Printer.hpp>

#include <sgpp/solver/sle/external/Armadillo.hpp>

#include <sgpp/globaldef.hpp>

#ifdef USE_ARMADILLO
#include <armadillo>

typedef arma::vec ArmadilloVector;
typedef arma::mat ArmadilloMatrix;
#endif /* USE_ARMADILLO */

#include <cstddef>
#include <iostream>
#include <string>

namespace sgpp {
namespace solver {

Armadillo::~Armadillo() {}

bool Armadillo::solve(base::SLE& system, base::DataVector& b, base::DataVector& x) const {
  base::DataMatrix B(b.getPointer(), b.getSize(), 1);
  base::DataMatrix X(B.getNrows(), B.getNcols());

  // call version for multiple RHSs
  if (solve(system, B, X)) {
    x.resize(X.getNrows());
    X.getColumn(0, x);
    return true;
  } else {
    return false;
  }
}

bool Armadillo::solve(base::SLE& system, base::DataMatrix& B, base::DataMatrix& X) const {
#ifdef USE_ARMADILLO
  base::Printer::getInstance().printStatusBegin("Solving linear system (Armadillo)...");

  const arma::uword n = static_cast<arma::uword>(system.getDimension());
  ArmadilloMatrix A(n, n);
  size_t nnz = 0;
  size_t rowsDone = 0;

  A.zeros();

// parallelize only if the system is cloneable
#pragma omp parallel if (system.isCloneable()) shared(system, A, nnz, rowsDone, n)  // default(none)

  {
    base::SLE* system2 = &system;
#ifdef _OPENMP
    std::unique_ptr<base::CloneableSLE> clonedSLE;

    if (system.isCloneable() && (omp_get_max_threads() > 1)) {
      dynamic_cast<base::CloneableSLE&>(system).clone(clonedSLE);
      system2 = clonedSLE.get();
    }

#endif /* _OPENMP */

// copy system matrix to Armadillo matrix object
#pragma omp for ordered schedule(static)

    for (arma::uword i = 0; i < n; i++) {
      for (arma::uword j = 0; j < n; j++) {
        A(i, j) = system2->getMatrixEntry(i, j);

        // count nonzero entries
        // (not necessary, you can also remove that if you like)
        if (A(i, j) != 0) {
#pragma omp atomic
          nnz++;
        }
      }

#pragma omp atomic
      rowsDone++;

      // status message
      if (rowsDone % 100 == 0) {
        char str[10];
        snprintf(str, sizeof(str), "%.1f%%",
                 static_cast<double>(rowsDone) / static_cast<double>(n) * 100.0);
        base::Printer::getInstance().printStatusUpdate("constructing matrix (" + std::string(str) +
                                                       ")");
      }
    }
  }

  base::Printer::getInstance().printStatusUpdate("constructing matrix (100.0%)");
  base::Printer::getInstance().printStatusNewLine();

  // print ratio of nonzero entries
  {
    char str[10];
    double nnzRatio = static_cast<double>(nnz) / (static_cast<double>(n) * static_cast<double>(n));
    snprintf(str, sizeof(str), "%.1f%%", nnzRatio * 100.0);
    base::Printer::getInstance().printStatusUpdate("nnz ratio: " + std::string(str));
    base::Printer::getInstance().printStatusNewLine();
  }

  if (B.getNcols() == 1) {
    // only one RHS ==> use vector version of arma::solve
    ArmadilloVector bArmadillo(B.getPointer(), n);
    ArmadilloVector xArmadillo(n);

    base::Printer::getInstance().printStatusUpdate("solving with Armadillo");

    if (arma::solve(xArmadillo, A, bArmadillo)) {
      base::DataVector x(xArmadillo.memptr(), n);
      X.resize(n, 1);
      X.setColumn(0, x);
      base::Printer::getInstance().printStatusEnd();
      return true;
    } else {
      base::Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
      return false;
    }
  } else {
    // multiple RHSs ==> use matrix version of arma::solve
    const arma::uword B_count = static_cast<arma::uword>(B.getNcols());
    ArmadilloMatrix BArmadillo(n, B_count);
    ArmadilloMatrix XArmadillo(n, B_count);
    base::DataVector b(n);

    // copy RHSs to Armadillo matrix
    for (arma::uword i = 0; i < B_count; i++) {
      B.getColumn(i, b);
      BArmadillo.col(i) = ArmadilloVector(b.getPointer(), n);
    }

    base::Printer::getInstance().printStatusUpdate("solving with Armadillo");

    if (arma::solve(XArmadillo, A, BArmadillo)) {
      X.resize(n, B_count);

      // convert solutions to base::DataVector
      for (arma::uword i = 0; i < B_count; i++) {
        base::DataVector x(XArmadillo.colptr(i), n);
        X.setColumn(i, x);
      }

      base::Printer::getInstance().printStatusEnd();
      return true;
    } else {
      base::Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
      return false;
    }
  }

#else
  std::cerr << "Error in solver::Armadillo::solve: "
            << "SG++ was compiled without Armadillo support!\n";
  return false;
#endif /* USE_ARMADILLO */
}
}  // namespace solver
}  // namespace sgpp
