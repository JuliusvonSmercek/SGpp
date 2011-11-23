/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef OPERATIONMULTIPLEEVALPREWAVELET_HPP
#define OPERATIONMULTIPLEEVALPREWAVELET_HPP

#include "datadriven/operation/OperationMultipleEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg
{
namespace base
{

/**
 * This class implements OperationMultipleEval for a grids with prewavelet ansatzfunctions without boundaries
 */
class OperationMultipleEvalPrewavelet : public OperationMultipleEval
{
public:
	/**
	 * Constructor of OperationMultipleEvalPrewavelet
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationMultipleEvalPrewavelet(GridStorage* storage, DataMatrix* dataset) : OperationMultipleEval(dataset) {
		this->storage = storage;
	}

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalPrewavelet() {}

	virtual void mult(DataVector& alpha, DataVector& result);
	virtual void multTranspose(DataVector& source, DataVector& result);

protected:
	/// Pointer to the grid's GridStorage object
	GridStorage* storage;
};

}
}

#endif /* OPERATIONMULTIPLEEVALPREWAVELET_HPP */