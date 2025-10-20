// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/IterativeGaussianElimination.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace sgpp {
namespace base {
namespace sle_solver {

namespace {
// Use full pivoting for maximum numerical stability.
const bool FULL_PIVOTING = true;

/**
 * @brief Performs an iterative in-place LU decomposition with full pivoting.
 *
 * Updates an existing LU decomposition of size `oldn x oldn` to a new decomposition
 * of size `N x N`.
 *
 * @param A         Matrix to be decomposed. Initially contains the old LU factors and new raw data.
 * @param pivotRow  Row pivot vector to be updated.
 * @param pivotCol  Column pivot vector to be updated.
 * @param oldn      The size of the existing decomposition.
 * @param tolerance A small tolerance for detecting singularity.
 * @return          True on success, false if the matrix is singular.
 */
bool iterative_lu_decomposition(DataMatrix& A, std::vector<size_t>& pivotRow,
                                std::vector<size_t>& pivotCol, const size_t oldn,
                                const double tolerance) {
  const size_t N = A.getNcols();

  for (size_t i = 0; i < oldn; ++i) {
    for (size_t k = oldn; k < N; ++k) {
      A(k, i) /= A(i, i);

      for (size_t j = i + 1; j < oldn; ++j) {
        A(k, j) -= A(k, i) * A(i, j);
        A(j, k) -= A(j, i) * A(i, k);
      }

      for (size_t j = oldn; j < N; ++j) {
        A(k, j) -= A(k, i) * A(i, j);
      }
    }
  }

  for (size_t k = oldn; k < N - 1; ++k) {
    double maxA = 0.0;
    size_t imax = k;
    size_t jmax = k;

    for (size_t i = k; i < N; ++i) {
      if (FULL_PIVOTING) {
        for (size_t j = k; j < N; ++j) {
          const double absA = std::abs(A(i, j));
          if (maxA < absA) {
            maxA = absA;
            imax = i;
            jmax = j;
          }
        }
      } else {  // Partial pivoting
        const double absA = std::abs(A(i, k));
        if (maxA < absA) {
          maxA = absA;
          imax = i;
        }
      }
    }

    if (maxA < tolerance) {
      return false;  // Matrix is singular.
    }

    if (imax != k) {
      std::swap(pivotRow[k], pivotRow[imax]);
      for (size_t i = 0; i < N; ++i) {
        std::swap(A(k, i), A(imax, i));
      }
    }
    if (jmax != k) {
      std::swap(pivotCol[k], pivotCol[jmax]);
      for (size_t j = 0; j < N; ++j) {
        std::swap(A(j, k), A(j, jmax));
      }
    }

    if (std::abs(A(k, k)) < tolerance) {
      return false;  // Pivoting led to a zero pivot, matrix is singular.
    }

    for (size_t i = k + 1; i < N; ++i) {
      A(i, k) /= A(k, k);
      for (size_t j = k + 1; j < N; ++j) {
        A(i, j) -= A(i, k) * A(k, j);
      }
    }
  }

  return true;
}

/**
 * @brief Prepares the LU matrix and pivot vectors for an iterative decomposition.
 *
 * @param A         The full, updated source matrix.
 * @param LU        The matrix to prepare (will contain old LU + new raw data).
 * @param oldn      The size of the existing decomposition.
 * @param pivotRow  Row pivot vector.
 * @param pivotCol  Column pivot vector.
 * @param tolerance A small tolerance for detecting singularity.
 * @return          True on success, false if the matrix is singular.
 */
bool lu_decomposition_update(const DataMatrix& A, DataMatrix& LU, const size_t oldn,
                             std::vector<size_t>& pivotRow, std::vector<size_t>& pivotCol,
                             const double tolerance) {
  const size_t N = A.getNcols();

  size_t oldPivotSize = pivotRow.size();
  pivotRow.resize(N);
  if (oldPivotSize < N) {
    std::iota(pivotRow.begin() + oldPivotSize, pivotRow.end(), oldPivotSize);
  }

  oldPivotSize = pivotCol.size();
  pivotCol.resize(N);
  if (oldPivotSize < N) {
    std::iota(pivotCol.begin() + oldPivotSize, pivotCol.end(), oldPivotSize);
  }

  LU.resizeQuadratic(N);
  // Copy new data from A into LU, applying old pivot permutations.
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = oldn; j < N; ++j) {
      LU(i, j) = A(pivotRow[i], pivotCol[j]);
      LU(j, i) = A(pivotRow[j], pivotCol[i]);
    }
  }

  return iterative_lu_decomposition(LU, pivotRow, pivotCol, oldn, tolerance);
}

/**
 * @brief Solves Ax=b using a pre-computed LU decomposition (PAQ = LU).
 *
 * @param LU        The decomposed matrix.
 * @param pivotRow  Row pivot vector P.
 * @param pivotCol  Column pivot vector Q.
 * @param b         Right-hand side vector.
 * @return          The solution vector x.
 */
DataVector lu_decomposition_solve(const DataMatrix& LU, const std::vector<size_t>& pivotRow,
                                  const std::vector<size_t>& pivotCol, const DataVector& b) {
  const size_t N = LU.getNcols();
  DataVector y(N);

  // Forward substitution: Ly = Pb
  for (size_t i = 0; i < N; ++i) {
    y[i] = b[pivotRow[i]];
    for (size_t k = 0; k < i; ++k) {
      y[i] -= LU(i, k) * y[k];
    }
  }

  // Backward substitution: Ux = y
  for (int i = N - 1; i >= 0; --i) {
    for (size_t k = i + 1; k < N; ++k) {
      y[i] -= LU(i, k) * y[k];
    }
    y[i] /= LU(i, i);
  }

  // Undo column permutations: x = Q*x_perm
  DataVector x(N);
  for (size_t i = 0; i < N; ++i) {
    x[pivotCol[i]] = y[i];
  }

  return x;
}
}  // namespace

IterativeGaussianElimination::IterativeGaussianElimination() : A(), LU(), pivotRow(), pivotCol() {}

IterativeGaussianElimination::~IterativeGaussianElimination() {}

bool IterativeGaussianElimination::iterativeSolve(SLE& system, DataVector& b, DataVector& x) {
  Printer::getInstance().printStatusBegin(
      "Solving linear system (Iterative Gaussian elimination)...");

  // old size of the system
  const size_t oldn = A.getNcols();
  // size of the system
  const size_t n = system.getDimension();

  if (b.getSize() != n) {
    Printer::getInstance().printStatusEnd("Error: right-hand side vector has incorrect dimension.");
    return false;
  } else if (0 == n) {
    x.resize(0);
    Printer::getInstance().printStatusEnd();
    return true;
  }

  // Update internal copy of the system matrix A.
  // The old part of the matrix is assumed to be unchanged for efficiency.
  A.resizeQuadratic(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < oldn; j++) {
      if (system.getMatrixEntry(i, j) != A(i, j) || system.getMatrixEntry(j, i) != A(j, i)) {
        Printer::getInstance().printStatusEnd("Error: Matrix is different than last time.");
      }
    }
    for (size_t j = oldn; j < n; ++j) {
      A(i, j) = system.getMatrixEntry(i, j);
      if (i != j) {
        A(j, i) = system.getMatrixEntry(j, i);
      }
    }
  }

  bool success = lu_decomposition_update(A, LU, oldn, pivotRow, pivotCol, DEGENERATION_TOLERANCE);

  if (!success) {
    Printer::getInstance().printStatusUpdate("Iterative update failed. Re-solving from scratch...");
    // Resetting LU and pivot information by passing oldn=0, which re-initializes everything.
    success = lu_decomposition_update(A, LU, 0, pivotRow, pivotCol, DEGENERATION_TOLERANCE);
  }

  if (success) {
    x = lu_decomposition_solve(LU, pivotRow, pivotCol, b);
  } else {
    Printer::getInstance().printStatusUpdate("Error: Matrix is singular. Could not solve the system.");
  }

  Printer::getInstance().printStatusEnd();
  return success;
}
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
