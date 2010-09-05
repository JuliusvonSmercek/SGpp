/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BLACKSCHOLESHULLWHITESOLVER_HPP
#define BLACKSCHOLESHULLWHITESOLVER_HPP

#include "sgpp.hpp"

#include "application/pde/ParabolicPDESolver.hpp"

#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/common/BoundingBox.hpp"

#include "tools/common/StdNormalDistribution.hpp"

#include "application/common/ScreenOutput.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace sg
{

/**
 * This class provides a simple-to-use solver of the multi dimensional Black
 * Scholes Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the Black Scholes
 * Equation with Sparse Grids!
 *
 * @version $HEAD$
 */
class BlackScholesHullWhiteSolver : public ParabolicPDESolver
{
private:
	/// vector that contains the assets' weight
	DataVector* mus;
	/// vector that contains the standard deviations
	DataVector* sigmas;
	/// Matrix that contains the correlations
	DataMatrix* rhos;
	/// the riskfree rate
	double r;
	/// stores if the stochastic asset data was passed to the solver
	bool bStochasticDataAlloc;
	/// screen object used in this solver
	ScreenOutput* myScreen;
	/// use coarsening between timesteps in order to reduce gridsize
	bool useCoarsen;
	/// Threshold used to decide if a grid point should be deleted
	double coarsenThreshold;
	/// Threshold used to decide if a grid point should be refined
	double refineThreshold;
	/// adaptive mode during solving Black Scholes Equation: none, coarsen, refine, coarsenNrefine
	std::string adaptSolveMode;
	/// refine mode during solving Black Scholes Equation: classic or maxLevel
	std::string refineMode;
	/// number of points the are coarsened in each coarsening-step
	int numCoarsenPoints;
	/// identifies if the Black Scholes Equation should be solved on a log-transformed grid
	bool useLogTransform;
	/// max. level for refinement during solving
	size_t refineMaxLevel;
	/// variable to store needed solving iterations
	size_t nNeededIterations;
	/// variable to store the solving time
	double dNeededTime;
	/// variable to store start grid size (Inner Grid)
	size_t staInnerGridSize;
	/// variable to store final grid size (Inner Grid)
	size_t finInnerGridSize;
	/// variable to store average grid size (Inner Grid)
	size_t avgInnerGridSize;
	/// Percent how many of the removable points should be tested for deletion
	double coarsenPercent;
	/// denotes the number of coarsening procedures within one timestep
	size_t numExecCoarsen;
	double theta;
	/// the standard deviations
	double sigma;
	///
	double a;
	/**
	 * returns the option value (payoff value) for an European call option
	 *
	 * @param assetValue the current asset's value
	 * @param strike the strike price of the option
	 *
	 * @return the call premium
	 */
	//double get1DEuroCallPayoffValue(double assetValue, double strike);

	/**
	 * Inits the alpha vector with a payoff function of an European call option or put option.
	 * The grid is initialized based on Cartesian coordinates!
	 *
	 * @param alpha the coefficient vector of the grid's ansatzfunctions
	 * @param strik the option's strike
	 * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
	 */
	//void initCartesianGridWithPayoff(DataVector& alpha, double strike, std::string payoffType);

	/**
	 * Inits the alpha vector with a payoff function of an European call option or put option
	 * The grid is initialized based on log-transformed coordinates!
	 *
	 * @param alpha the coefficient vector of the grid's ansatzfunctions
	 * @param strik the option's strike
	 * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
	 */
	//void initLogTransformedGridWithPayoff(DataVector& alpha, double strike, std::string payoffType);

public:
	/**
	 * Std-Constructor of the solver
	 */
	BlackScholesHullWhiteSolver(bool useLogTransform = false);

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~BlackScholesHullWhiteSolver();

	void constructGrid(BoundingBox& myBoundingBox, size_t level);

	/**
	 * This function tries to refine the grid such that
	 * most of the grid points are used for interpolation of the singularity. So this grid
	 * is able to approximate the start solution better.
	 *
	 * After refining the grid the payoff function is applied to the grid.
	 *
	 * Only on Cartesian grids!
	 *
	 * @param alpha reference to a DataVector object that contains the gird ansatzfunction's coefficients
	 * @param strike containing the option's strike
	 * @param payoffType the type of payoff Function used ONLY supported: avgM
	 * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
	 */
	//void refineInitialGridWithPayoff(DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance);

	/**
	 * This function tries to refine the grid such that
	 * most of the grid points are used for interpolation of the singularity. So this grid
	 * is able to approximate the start solution better. Refining is done only if the max
	 * refinement level hasn't be reached.
	 *
	 * After refining the grid the payoff function is applied to the grid.
	 *
	 * Only on Cartesian grids!
	 *
	 * @param alpha reference to a DataVector object that contains the gird ansatzfunction's coefficients
	 * @param strike containing the option's strike
	 * @param payoffType the type of payoff Function used ONLY supported: avgM
	 * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
	 * @param maxLevel maximum level of refinement
	 */
	//void refineInitialGridWithPayoffToMaxLevel(DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance, size_t maxLevel);

	/**
	 * In order to solve the multi dimensional Black Scholes Equation you have to provided
	 * some statistical data about the underlying (assets' weight, standard deviation
	 * and the correlation between them). This function allows you to set this data.
	 *
	 * @param mus a DataVector that contains the underlyings' weight
	 * @param sigmas a DataVector that contains the underlyings' standard deviations
	 * @param rhos a DataMatrix that contains the correlations between the underlyings
	 * @param r the riskfree rate used in the market model
	 */
	void setStochasticData(DataVector& mus, DataVector& sigmas, DataMatrix& rhos, double r);

	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul = 0);

	/**
	 * Solves the closed form of the Black Scholes equation, the Black Scholes
	 * formular. It evaluates the Black Scholes formular in a Stock Price Range
	 * from 0.0 to maxStock, by increasing the stock price in every step by
	 * a given (small) values, so the analytical solution of the PDE can
	 * be determined and compared.
	 *
	 * @param premiums the result vector, here the combinations of stock price and premium are stored
	 * @param maxStock the maximum stock regarded in these calculations
	 * @param StockInc the increase of the stockprice in one step
	 * @param strike the strike price of the Option
	 * @param t time to maturity
	 */
	//void solve1DAnalytic(std::vector< std::pair<double, double> >& premiums, double maxStock, double StockInc, double strike, double t);

	/**
	 * Writes the premiums into a file that can be easily plot with gnuplot
	 *
	 * @param premiums the result vector, here the combinations of stock price and premium are stored
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
//	void print1DAnalytic(std::vector< std::pair<double, double> >& premiums, std::string tfilename);

	/**
	 * Inits the alpha vector with a payoff function of an European call option or put option
	 *
	 * @param alpha the coefficient vector of the grid's ansatzfunctions
	 * @param strik the option's strike
	 * @param payoffType specifies the type of the combined payoff function; std_euro_call or std_euro_put are available
	 */
	//void initGridWithPayoff(DataVector& alpha, double strike, std::string payoffType);
	void initGridWithPayoffBSHW(DataVector& alpha, double strike, std::string payoffType);
	/**
	 * Inits the screen object
	 */
	void initScreen();

	/**
	 * returns the algorithmic dimensions (the dimensions in which the Up Down
	 * operations (need for space discretization) should be applied)
	 *
	 * @return the algorithmic dimensions
	 */
	std::vector<size_t> getAlgorithmicDimensions();

	/**
	 * sets the algorithmic dimensions (the dimensions in which the Up Down
	 * operations (need for space discretization) should be applied)
	 *
	 * @param algoDims std::vector containing the algorithmic dimensions
	 */
	void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

	/**
	 *	enables coarsening of grid during solving the Black Scholes
	 *	Equation. The coarsening settings have to be specified in order to
	 *	enable coarsening.
	 *
	 *	@param coarsenThreshold Threshold needed to determine if a grid point should be removed
	 *	@param refineMode the Mode used for refining the grid: classic or maxLevel
	 *	@param refineMaxLevel max. level for refinement during solving
	 *	@param adaptSolveMode adaptive mode during solving equation: coarsen, refine, coarsenNrefine
	 *	@param numCoarsenPoints number of points coarsened, -1 all coarsenable points are coarsened
	 *	@param refineThreshold Threshold needed to determine if a grid point should be refined
	 */
	void setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, size_t refineMaxLevel, int numCoarsenPoints, double coarsenThreshold, double refineThreshold);

	/**
	 * prints the 2D interpolation error @money into a file. This file is plotable via gnuplot. A bounding
	 * box [0,x] X [0,y] is assumed.
	 *
	 * Only on Cartesian grids!
	 *
	 * @param alpha the sparse grid's coefficients
	 * @param tFilename the name of file contain the interpolation error
	 * @param numTestpoints Number of equal distribute testpoints @money
	 * @param strike the option's strike
	 */
	//void printPayoffInterpolationError2D(DataVector& alpha, std::string tFilename, size_t numTestpoints, double strike);

	/**
	 * gets the number of gridpoints @money
	 *
	 * Only on Cartesian grids!
	 *
	 * @param payoffType the payoff type
	 * @param strike the option's strike
	 * @param eps epsilon to determine the gridpoints, use if @money is not exactly on grid
	 *
	 * @param number of gridpoints @money
	 */
	size_t getGridPointsAtMoney(std::string payoffType, double strike, double eps = 0.0);

	/**
	 * gets the number needed iterations to solve Black Scholes Equation
	 *
	 * @return number of iterations needed to solve Black Scholes Equation, if called before solving 0 is returned
	 */
	size_t getNeededIterationsToSolve();

	/**
	 * gets needed time in seconds to solve Black Scholes Equation
	 *
	 * @return needed time in seconds to solve Black Scholes Equation, if called before solving 0 is returned
	 */
	double getNeededTimeToSolve();

	/**
	 * gets the number of points in start grid
	 *
	 * @returns the number of points in start grid, if called before constructing grid, 0 is returned
	 */
	size_t getStartInnerGridSize();

	/**
	 * gets the number of points in final grid
	 *
	 * @returns the number of points in final grid, if called before solving, 0 is returned
	 */
	size_t getFinalInnerGridSize();

	/**
	 * gets the number of average gridpoints
	 *
	 * @returns the number of average gridpoints, if called before solving, 0 is returned
	 */
	size_t getAverageInnerGridSize();
};

}

#endif /* BLACKSCHOLESSOLVER_HPP */