// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%rename(SolverModuleSLESolver)          sgpp::solver::SLESolver;
%rename(SolverModuleBiCGStab)           sgpp::solver::BiCGStab;

// The Good, i.e. without any modifications
%include "solver/src/sgpp/solver/common/SGSolver.hpp"
%include "solver/src/sgpp/solver/sle/common/SLESolver.hpp"
%include "solver/src/sgpp/solver/ode/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "solver/src/sgpp/solver/sle/native/iterative/ConjugateGradients.hpp"
%include "solver/src/sgpp/solver/sle/native/iterative/BiCGStab.hpp"
%include "solver/src/sgpp/solver/ode/integrators/Euler.hpp"
%include "solver/src/sgpp/solver/ode/integrators/CrankNicolson.hpp"
%include "solver/src/sgpp/solver/sle/common/TypesSolver.hpp"
%include "solver/src/sgpp/solver/sle/common/SLESolverTypeParser.hpp"

%include "solver/src/sgpp/solver/pde/OperationParabolicPDESolverSystem.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
