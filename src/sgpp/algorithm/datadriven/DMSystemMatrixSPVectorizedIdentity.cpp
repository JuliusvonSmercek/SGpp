/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "exception/operation_exception.hpp"

namespace sg
{

DMSystemMatrixSPVectorizedIdentity::DMSystemMatrixSPVectorizedIdentity(Grid& SparseGrid, DataMatrixSP& trainData, float lambda, std::string vecMode)
{
	// handle unsupported vector extensions
	if (vecMode != "SSE" && vecMode != "AVX")
	{
		throw new operation_exception("DMSystemMatrixVectorizedIdentity : Only SSE or AVX are supported vector extensions!");
	}

	// create the operations needed in ApplyMatrix
	this->vecMode = vecMode;
	this->lamb = lambda;
	this->B = SparseGrid.createOperationBVectorizedSP(this->vecMode);
	this->data = new DataMatrixSP(trainData);

	if (this->vecMode == "SSE")
	{
		this->vecWidth = 4;
	}
	else if (this->vecMode == "AVX")
	{
		this->vecWidth = 8;
	}
	// should not happen because this exception should have been thrown some lines upwards!
	else
	{
		throw new operation_exception("DMSystemMatrixVectorizedIdentity : Only SSE or AVX are supported vector extensions!");
	}

	numTrainingInstances = data->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = data->getNrows() % this->vecWidth;
	size_t loopCount = this->vecWidth - remainder;

	if (loopCount != this->vecWidth)
	{
		DataVectorSP lastRow(data->getNcols());
		for (size_t i = 0; i < loopCount; i++)
		{
			data->getRow(data->getNrows()-1, lastRow);
			data->resize(data->getNrows()+1);
			data->setRow(data->getNrows()-1, lastRow);
		}
	}
	data->transpose();
}

DMSystemMatrixSPVectorizedIdentity::~DMSystemMatrixSPVectorizedIdentity()
{
	delete this->B;
	delete this->data;
}

void DMSystemMatrixSPVectorizedIdentity::mult(DataVectorSP& alpha, DataVectorSP& result)
{
	DataVectorSP temp((*data).getNcols());

    // Operation B
    this->B->multTransposeVectorized(alpha, (*data), temp);
    // patch result -> set additional entries zero
    if (numTrainingInstances != temp.getSize())
    {
    	for (size_t i = 0; i < (temp.getSize()-numTrainingInstances); i++)
    	{
    		temp.set(temp.getSize()-(i+1), 0.0);
    	}
    }
    this->B->multVectorized(temp, (*data), result);

    result.axpy(numTrainingInstances*this->lamb, alpha);
}

void DMSystemMatrixSPVectorizedIdentity::generateb(DataVectorSP& classes, DataVectorSP& b)
{
	DataVectorSP myClasses(classes);

	// Apply padding
	size_t numCols = (*data).getNcols();
	if (numCols != myClasses.getSize())
	{
		myClasses.resizeZero(numCols);
	}
	this->B->multVectorized(myClasses, (*data), b);
}

void DMSystemMatrixSPVectorizedIdentity::rebuildLevelAndIndex()
{
	this->B->rebuildLevelAndIndex();
}

}