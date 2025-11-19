// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <sgpp/solver/sle/native/iterative/ConjugateGradients.hpp>
#include <sgpp/solver/sle/native/iterative/BiCGStab.hpp>
#include <sgpp/solver/ode/integrators/Euler.hpp>
#include <sgpp/solver/ode/integrators/CrankNicolson.hpp>
#include <sgpp/solver/ode/integrators/AdamsBashforth.hpp>
#include <sgpp/solver/ode/step_control/VarTimestep.hpp>
#include <sgpp/solver/ode/step_control/StepsizeControl.hpp>
#include <sgpp/solver/ode/step_control/StepsizeControlEJ.hpp>
#include <sgpp/solver/ode/step_control/StepsizeControlH.hpp>
#include <sgpp/solver/ode/step_control/StepsizeControlMC.hpp>
#include <sgpp/solver/ode/step_control/StepsizeControlBDF.hpp>
#include <sgpp/solver/sle/common/TypesSolver.hpp>
#include <sgpp/solver/sle/common/SLESolverTypeParser.hpp>
#include <sgpp/solver/pde/OperationParabolicPDESolverSystem.hpp>

#endif /* SOLVER_HPP */
