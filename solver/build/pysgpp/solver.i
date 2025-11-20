// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%rename(SolverModuleSLESolver)          sgpp::solver::IterativeSLESolver;
%rename(SolverModuleBiCGStab)           sgpp::solver::BiCGStab;

// The Good, i.e. without any modifications
%include "solver/src/sgpp/solver/common/IterativeSGSolver.hpp"
%include "solver/src/sgpp/solver/sle/common/IterativeSLESolver.hpp"
%include "solver/src/sgpp/solver/ode/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "solver/src/sgpp/solver/sle/native/iterative/ConjugateGradients.hpp"
%include "solver/src/sgpp/solver/sle/native/iterative/BiCGStab.hpp"
%include "solver/src/sgpp/solver/ode/integrators/Euler.hpp"
%include "solver/src/sgpp/solver/ode/integrators/CrankNicolson.hpp"
%include "solver/src/sgpp/solver/sle/common/TypesSolver.hpp"
%include "solver/src/sgpp/solver/sle/common/SLESolverTypeParser.hpp"

%include "solver/src/sgpp/solver/pde/OperationParabolicPDESolverSystem.hpp"

%include "solver/src/sgpp/solver/sle/common/SLESolver.hpp"
%include "solver/src/sgpp/solver/sle/external/Armadillo.hpp"
%rename(AutoSLESolver) sgpp::solver::Auto;
%include "solver/src/sgpp/solver/sle/external/Auto.hpp"
%include "solver/src/sgpp/solver/sle/native/iterative/MyBiCGStab.hpp"
%include "solver/src/sgpp/solver/sle/external/Eigen.hpp"
%include "solver/src/sgpp/solver/sle/native/direct/GaussianElimination.hpp"
%include "solver/src/sgpp/solver/sle/native/direct/IterativeGaussianElimination.hpp"
%include "solver/src/sgpp/solver/sle/external/Gmmpp.hpp"
%include "solver/src/sgpp/solver/sle/external/UMFPACK.hpp"

// global variables for the support of SLE solver libaries (set at compile-time)
const bool ARMADILLO_ENABLED;
const bool EIGEN_ENABLED;
const bool GMMPP_ENABLED;
const bool UMFPACK_ENABLED;

%{
#ifdef USE_ARMADILLO
    const bool ARMADILLO_ENABLED = true;
#else
    const bool ARMADILLO_ENABLED = false;
#endif

#ifdef USE_EIGEN
    const bool EIGEN_ENABLED = true;
#else
    const bool EIGEN_ENABLED = false;
#endif

#ifdef USE_GMMPP
    const bool GMMPP_ENABLED = true;
#else
    const bool GMMPP_ENABLED = false;
#endif

#ifdef USE_UMFPACK
    const bool UMFPACK_ENABLED = true;
#else
    const bool UMFPACK_ENABLED = false;
#endif
%}

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
