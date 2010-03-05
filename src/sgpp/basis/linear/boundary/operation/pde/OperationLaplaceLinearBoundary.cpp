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

#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include "grid/common/BoundingBox.hpp"

namespace sg
{

OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(GridStorage* storage) : UnidirGradient(storage)
{
}

OperationLaplaceLinearBoundary::~OperationLaplaceLinearBoundary()
{
}

void OperationLaplaceLinearBoundary::mult(DataVector& alpha, DataVector& result)
{
	this->updown(alpha, result);
}

void OperationLaplaceLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// Bounding Box handling
	BoundingBox* boundingBox = this->storage->getBoundingBox();
	double q = boundingBox->getIntervalWidth(dim);
	double Qqout = 1.0/(q*q);

	// init the coefficients of the ansatz functions with boundary
	result.setAll(0.0);

	if (q != 1.0)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				if (index == 0)
				{
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += Qqout*alpha[i];

						// down
						if (index == 0)
						{
							GridIndex index_one = (*storage)[i];
							index_one.set(dim, 0, 1);
							if (!boundingBox->hasDirichletBoundaryRight(dim))
							{
								result[(*storage)[&index_one]] += ((-1.0*Qqout) * alpha[i]);
							}
						}
					}
				}
				if (index == 1)
				{
					if (!boundingBox->hasDirichletBoundaryRight(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += Qqout*alpha[i];
					}
				}
			}
			//only affects the diagonal of the stiffness matrix
			else
			{
				result[i] = alpha[i]*(Qqout*pow(2.0, static_cast<int>(level+1)));
			}
		}
	}
	else
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				if (index == 0)
				{
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += alpha[i];

						// down
						if (index == 0)
						{
							GridIndex index_one = (*storage)[i];
							index_one.set(dim, 0, 1);
							if (!boundingBox->hasDirichletBoundaryRight(dim))
							{
								result[(*storage)[&index_one]] += ((-1.0) * alpha[i]);
							}
						}
					}
				}
				if (index == 1)
				{
					if (!boundingBox->hasDirichletBoundaryRight(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += alpha[i];
					}
				}
			}
			//only affects the diagonal of the stiffness matrix
			else
			{
				result[i] = alpha[i]*pow(2.0, static_cast<int>(level+1));
			}
		}
	}
}

void OperationLaplaceLinearBoundary::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// Bounding Box handling
	BoundingBox* boundingBox = this->storage->getBoundingBox();
	double q = boundingBox->getIntervalWidth(dim);
	double Qqout = 1.0/(q*q);

	// init the coefficients of the ansatz functions with boundary
	result.setAll(0.0);

	if (q != 1.0)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				// up
				if (index == 1)
				{
					GridIndex index_zero = (*storage)[i];
					index_zero.set(dim, 0, 0);
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						result[(*storage)[&index_zero]] += ((-1.0*Qqout) * alpha[i]);
					}
				}
			}
		}
	}
	else
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				// up
				if (index == 1)
				{
					GridIndex index_zero = (*storage)[i];
					index_zero.set(dim, 0, 0);
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						result[(*storage)[&index_zero]] += ((-1.0) * alpha[i]);
					}
				}
			}
		}
	}
}

}