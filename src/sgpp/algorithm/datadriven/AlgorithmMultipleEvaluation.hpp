/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ALGORITHMMULTIPLEEVALUATION_HPP
#define ALGORTIHMMULTIPLEEVALUATION_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "algorithm/common/AlgorithmEvaluation.hpp"
#include "algorithm/common/AlgorithmEvaluationTransposed.hpp"
#include "algorithm/common/AlgorithmEvaluationIterative.hpp"

#include <vector>
#include <utility>
#include <iostream>

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg {

/**
 * Basic multiplaction with B and B^T on grids with no boundaries.
 * If there are @f$N@f$ basis functions @f$\varphi(\vec{x})@f$ and @f$m@f$ data points, then B is a (Nxm) matrix, with
 * @f[ (B)_{i,j} = \varphi_i(x_j). @f]
 *
 * @todo (blank) check if it is possible to have some functor for the BASIS type
 */
template<class BASIS>
class AlgorithmMultipleEvaluation
{
public:

	/**
	 * Performs the DGEMV Operation on the grid
	 *
	 * This operation can be executed in parallel by setting the USEOMP define
	 *
	 * @todo (heinecke, nice) add mathematical description
	 *
	 * @param storage GridStorage object that contains the grid's points information
	 * @param basis a reference to a class that implements a specific basis
	 * @param source the coefficients of the grid points
	 * @param x the d-dimensional vector with data points (row-wise)
	 * @param result the result vector of the matrix vector multiplication
	 */
	void mult(GridStorage* storage, BASIS& basis, DataVector& source, DataVector& x, DataVector& result)
	{
		result.setAll(0.0);
		size_t source_size = source.getSize();

#ifdef USEOMP
		#pragma omp parallel
		{
			DataVector privateResult(result);
			std::vector<double> line;
			AlgorithmEvaluationTransposed<BASIS> AlgoEvalTrans(storage);

			privateResult.setAll(0.0);

			#pragma omp for schedule(static)
			for(size_t i = 0; i < source_size; i++)
			{
				x.getLine(i, line);

				AlgoEvalTrans(basis, line, source[i], privateResult);
			}

			#pragma omp critical
			{
				result.add(privateResult);
			}
		}
#else
		std::vector<double> line;
		AlgorithmEvaluationTransposed<BASIS> AlgoEvalTrans(storage);

		for(size_t i = 0; i < source_size; i++)
		{
			x.getLine(i, line);

			AlgoEvalTrans(basis, line, source[i], result);
		}
#endif /* USEOMP */
	}

	/**
	 * Performs the DGEMV Operation on the grid having a transposed matrix
	 *
	 * This operation can be executed in parallel by setting the USEOMP define
	 *
	 * @todo (heinecke, nice) add mathematical description
	 *
	 * @param storage GridStorage object that contains the grid's points information
	 * @param basis a reference to a class that implements a specific basis
	 * @param source the coefficients of the grid points
	 * @param x the d-dimensional vector with data points (row-wise)
	 * @param result the result vector of the matrix vector multiplication
	 */
	void mult_transpose(GridStorage* storage, BASIS& basis, DataVector& source, DataVector& x, DataVector& result)
	{
		result.setAll(0.0);
		size_t result_size = result.getSize();

//#ifdef USEOMP
//		#pragma omp parallel
//		{
//			std::vector<double> line;
//			AlgorithmEvaluationIterative<BASIS> AlgoEval(storage);
//
//			#pragma omp for schedule (static)
//			for(size_t i = 0; i < result_size; i++)
//			{
//				x.getLine(i, line);
//
//				result[i] = AlgoEval(basis, line, source);
//			}
//		}
//#else
//		std::vector<double> line;
//		AlgorithmEvaluationIterative<BASIS> AlgoEval(storage);
//
//		for(size_t i = 0; i < result_size; i++)
//		{
//			x.getLine(i, line);
//
//			result[i] = AlgoEval(basis, line, source);
//		}
//#endif /* USEOMP */

#ifdef USEOMP
		#pragma omp parallel
		{
			std::vector<double> line;
			AlgorithmEvaluation<BASIS> AlgoEval(storage);

			#pragma omp for schedule (static)
			for(size_t i = 0; i < result_size; i++)
			{
				x.getLine(i, line);

				result[i] = AlgoEval(basis, line, source);
			}
		}
#else
		std::vector<double> line;
		AlgorithmEvaluation<BASIS> AlgoEval(storage);

		for(size_t i = 0; i < result_size; i++)
		{
			x.getLine(i, line);

			result[i] = AlgoEval(basis, line, source);
		}
#endif /* USEOMP */
	}
};

}

#endif /* ALGORTIHMMULTIPLEEVALUATION_HPP */