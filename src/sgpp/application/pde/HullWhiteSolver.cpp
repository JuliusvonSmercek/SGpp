/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

//#include "algorithm/pde/HullWhiteODESolverSystem.hpp"
#include "application/pde/HullWhiteSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace sg
{

HullWhiteSolver::HullWhiteSolver() : ParabolicPDESolver()
{
	//this->bStochasticDataAlloc = false;
	this->bGridConstructed = false;
	this->myScreen = NULL;
	this->useCoarsen = false;
	this->coarsenThreshold = 0.0;
	this->coarsenPercent = 0.0;
	this->numExecCoarsen = 0;
}


HullWhiteSolver::~HullWhiteSolver()
{
	/*if (this->bStochasticDataAlloc)
	{
		delete this->mus;
		delete this->sigmas;
		delete this->rhos;
	}*/
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void HullWhiteSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
{
	this->dim = BoundingBox.getDimensions();
	this->levels = level;

	this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

	GridGenerator* myGenerator = this->myGrid->createGridGenerator();
	myGenerator->regular(this->levels);
	delete myGenerator;

	this->myBoundingBox = this->myGrid->getBoundingBox();
	this->myGridStorage = this->myGrid->getStorage();

	//std::string serGrid;
	//myGrid->serialize(serGrid);
	//std::cout << serGrid << std::endl;

	this->bGridConstructed = true;
}

/*void BlackScholesSolver::setStochasticData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r)
{
	this->mus = new DataVector(mus);
	this->sigmas = new DataVector(sigmas);
	this->rhos = new DataVector(rhos);
	this->r = r;

	bStochasticDataAlloc = true;
}*/

void HullWhiteSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	/*if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
		BlackScholesODESolverSystem* myBSSystem = new BlackScholesODESolverSystem(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ExEul", false, this->useCoarsen, this->coarsenThreshold, this->coarsenPercent, this->numExecCoarsen);
		SGppStopwatch* myStopwatch = new SGppStopwatch();
		double execTime;

		std::cout << "Using Explicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
		myStopwatch->start();
		myEuler->solve(*myCG, *myBSSystem, true, verbose);
		execTime = myStopwatch->stop();

		std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
		std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myBSSystem;
		delete myCG;
		delete myEuler;
		delete myStopwatch;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveExplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	}*/
}


void HullWhiteSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	/*if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
		BlackScholesODESolverSystem* myBSSystem = new BlackScholesODESolverSystem(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", false, this->useCoarsen, this->coarsenThreshold, this->coarsenPercent, this->numExecCoarsen);
		SGppStopwatch* myStopwatch = new SGppStopwatch();
		double execTime;

		std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
		myStopwatch->start();
		myEuler->solve(*myCG, *myBSSystem, true, verbose);
		execTime = myStopwatch->stop();

		std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
		std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myBSSystem;
		delete myCG;
		delete myEuler;
		delete myStopwatch;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	}*/
}

void HullWhiteSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul)
{
	/*if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
#ifdef USEOMPTHREE
		BlackScholesODESolverSystemParallelOMP* myBSSystem = new BlackScholesODESolverSystemParallelOMP(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "CrNic", this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->coarsenPercent, this->numExecCoarsen);
#else
		BlackScholesODESolverSystem* myBSSystem = new BlackScholesODESolverSystem(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "CrNic", this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->coarsenPercent, this->numExecCoarsen);
#endif
		SGppStopwatch* myStopwatch = new SGppStopwatch();
		double execTime;

		size_t numCNSteps;
		size_t numIESteps;

		numCNSteps = numTimesteps;
		if (numTimesteps > NumImEul)
		{
			numCNSteps = numTimesteps - NumImEul;
		}
		numIESteps = NumImEul;

		Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
		CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize, this->myScreen);

		myStopwatch->start();
		if (numIESteps > 0)
		{
			std::cout << "Using Implicit Euler to solve " << numIESteps << " timesteps:" << std::endl;
			myBSSystem->setODESolver("ImEul");
			myEuler->solve(*myCG, *myBSSystem, false, false);
		}
		myBSSystem->setODESolver("CrNic");
		std::cout << "Using Crank Nicolson to solve " << numCNSteps << " timesteps:" << std::endl << std::endl << std::endl << std::endl;
		myCN->solve(*myCG, *myBSSystem, true, false);
		execTime = myStopwatch->stop();

		std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
		std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		std::cout << "Average Grid size: " << static_cast<double>(myBSSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
		std::cout << "Average Grid size (Inner): " << static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myBSSystem;
		delete myCG;
		delete myCN;
		delete myEuler;
		delete myStopwatch;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
	}*/
}


void HullWhiteSolver::initGridWithPayoff(DataVector& alpha, double strike, std::string payoffType, double sigma, double a, double t, double T)
{
	double tmp;

	if (this->bGridConstructed)
	{
		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);
	        coordsStream >> tmp;
           /* double	dblFuncValues= exp((0.04*(t-T)) + 0.04*(1-exp(-a*(T-t)))/a - 1/(4*pow(a,3)) * pow(sigma,2)*pow((exp(-a*T)-exp(-a*t)),2)*(exp(2*a*t)-1)- tmp/a * (1-exp(-a*(T-t))));


			if (payoffType == "std_euro_call")
			{
				alpha[i] = std::max((dblFuncValues-strike), 0.0);
			}
			else if (payoffType == "std_euro_put")
			{
				alpha[i] = std::max(strike-dblFuncValues, 0.0);
			}*/
	        double* dblFuncValues = new double[1];

	        			for (size_t j = 0; j < 1; j++)
	        			{
	        				coordsStream >> tmp;

	        				dblFuncValues[j] = exp((0.04*(t-T)) + 0.04*(1-exp(-a*(T-t)))/a - 1/(4*pow(a,3)) * pow(sigma,2)*pow((exp(-a*T)-exp(-a*t)),2)*(exp(2*a*t)-1)- tmp/a * (1-exp(-a*(T-t))));
	        			}

	        			if (payoffType == "std_euro_call")
	        			{
	        				tmp = 0.0;
	        				for (size_t j = 0; j < 1; j++)
	        				{
	        					tmp += dblFuncValues[j];
	        				}
	        				alpha[i] = std::max(((tmp)-strike), 0.0);
	        			}
	        			else if (payoffType == "std_euro_put")
	        			{
	        				tmp = 0.0;
	        				for (size_t j = 0; j < dim; j++)
	        				{
	        					tmp += dblFuncValues[j];
	        				}
	        				alpha[i] = std::max(strike-((tmp)), 0.0);
	        			}
			else
			{
				throw new application_exception("HullWhiteSolver::initGridWithPayoff : An unknown payoff-type was specified!");
			}

	        delete[] dblFuncValues;

			//delete dblFuncValues;
		}

		OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HullWhiteSolver::initGridWithPayoff : A grid wasn't constructed before!");
	}
}
/*
double BlackScholesSolver::get1DEuroCallPayoffValue(double assetValue, double strike)
{
	if (assetValue <= strike)
	{
		return 0.0;
	}
	else
	{
		return assetValue - strike;
	}
}*/
/*
void BlackScholesSolver::solve1DAnalytic(std::vector< std::pair<double, double> >& premiums, double maxStock, double StockInc, double strike, double t)
{
	if (bStochasticDataAlloc)
	{
		double stock = 0.0;
		double vola = this->sigmas->get(0);
		StdNormalDistribution* myStdNDis = new StdNormalDistribution();

		for (stock = 0.0; stock <= maxStock; stock += StockInc)
		{
			double dOne = (log((stock/strike)) + ((this->r + (vola*vola*0.5))*(t)))/(vola*sqrt(t));
			double dTwo = dOne - (vola*sqrt(t));
			double prem = (stock*myStdNDis->getCumulativeDensity(dOne)) - (strike*myStdNDis->getCumulativeDensity(dTwo)*(exp((-1.0)*this->r*t)));

			premiums.push_back(std::make_pair(stock, prem));
		}

		delete myStdNDis;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solve1DAnalytic : Stochastic parameters weren't set!");
	}
} */
/*
void BlackScholesSolver::print1DAnalytic(std::vector< std::pair<double, double> >& premiums, std::string tfilename)
{
	typedef std::vector< std::pair<double, double> > printVector;
	std::ofstream fileout;

	fileout.open(tfilename.c_str());
	for(printVector::iterator iter = premiums.begin(); iter != premiums.end(); iter++)
	{
		fileout << iter->first << " " << iter->second << " " << std::endl;
	}
	fileout.close();
}
*/
/*
std::vector<size_t> BlackScholesSolver::getAlgorithmicDimensions()
{
	return this->myGrid->getAlgorithmicDimensions();
}*/

/*
void BlackScholesSolver::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims)
{
	this->myGrid->setAlgorithmicDimensions(newAlgoDims);
}
*/
void HullWhiteSolver::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Hull White Solver, 1.3.0", "TUM (C) 2009-2010, by Chao qi");
	this->myScreen->writeStartSolve("One dimensional Hull White Solver");
}
}