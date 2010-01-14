/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef OPERATIONTESTMODPOLY_HPP
#define OPERATIONTESTMODPOLY_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * This class implements OperationTest for a grids with mod poly basis ansatzfunctions with
 *
 * @version $HEAD$
 */
class OperationTestModPoly : public OperationTest
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationTestModPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestModPoly() {}

	virtual double test(DataVector& alpha, DataVector& data, DataVector& classes);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataVector& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Poly Basis object
	SModPolyBase base;
};

}

#endif /* OPERATIONTESTMODPOLY_HPP */