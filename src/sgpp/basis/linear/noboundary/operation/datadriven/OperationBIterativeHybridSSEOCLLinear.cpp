/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeHybridSSEOCLLinear.hpp"
#include "exception/operation_exception.hpp"
#include "common/AlignedMemory.hpp"

#ifdef USEOMP
#include "omp.h"
#endif

// This value is adjusted for a 2 socket Intel Westmere System (X5650) (SMT on) with 2 NVidia Fermis (GTX470)
#define PERCENT_CPUS 32
#define PERCENT_CPUS_MULT 20

#ifdef USEICCINTRINSICS
// include SSE3 intrinsics
#include <pmmintrin.h>

union doubleAbsMaskHybrid
{
   const double d;
   const __int64 i;

   doubleAbsMaskHybrid() : i(0x7FFFFFFFFFFFFFFF) {}
};

_MM_ALIGN16 const doubleAbsMaskHybrid absMaskHybridDouble;
static const __m128d abs2MaskHybridDouble = _mm_load1_pd( &absMaskHybridDouble.d );
#endif

namespace sg
{

OperationBIterativeHybridSSEOCLLinear::OperationBIterativeHybridSSEOCLLinear(GridStorage* storage) : storage(storage)
{
	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);

	myTimer = new SGppStopwatch();
	myOCLKernels = new OCLKernels();
}

OperationBIterativeHybridSSEOCLLinear::~OperationBIterativeHybridSSEOCLLinear()
{
	delete Level;
	delete Index;
	delete myTimer;
	delete myOCLKernels;
}

void OperationBIterativeHybridSSEOCLLinear::rebuildLevelAndIndex()
{
	delete Level;
	delete Index;

	Level = new DataMatrix(storage->size(), storage->dim());
	Index = new DataMatrix(storage->size(), storage->dim());

	storage->getLevelIndexArraysForEval(*Level, *Index);

	myOCLKernels->resetKernels();
}

double OperationBIterativeHybridSSEOCLLinear::multVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	size_t source_size = alpha.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrSource = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrLevel = this->Level->getPointer();
    double* ptrIndex = this->Index->getPointer();
    double* ptrGlobalResult = result.getPointer();

    if (data.getNrows() % 128 != 0 || source_size != data.getNrows())
    {
    	throw operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
    }

#ifdef USEOMPTHREE
    // split result into GPU and CPU partition
    size_t gpu_partition;
    if (storageSize < 40000)
    {
    	gpu_partition = (storageSize * (100-(PERCENT_CPUS_MULT)))/100;
    }
    else
    {
    	gpu_partition = (storageSize * (100-(PERCENT_CPUS)))/100;
    }
    size_t gpu_pad = gpu_partition % 128;
    gpu_partition -= gpu_pad;

    // Do on-demand transpose
	double* ptrTransData = new double[dims*source_size];

	#pragma omp parallel
    {
		#pragma omp single nowait
    	{
			#pragma omp task
    		{
    			myOCLKernels->multOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, gpu_partition);
    		}

			#pragma omp task
    		{
#ifdef USEICCINTRINSICS
    			for (size_t n = 0; n < source_size; n++)
    			{
    				for(size_t d = 0; d < dims; d++)
    				{
    					ptrTransData[(d*source_size)+n] = ptrData[(n*dims)+d];
    				}
    			}

    			for (size_t j = gpu_partition; j < storageSize; j++)
				{
					#pragma omp task firstprivate(j)
					{
						__m128d res = _mm_set1_pd(0.0f);

						for (size_t i = 0; i < source_size; i+=8)
						{
							__m128d support_0 = _mm_load_pd(&(ptrSource[i+0]));
							__m128d support_1 = _mm_load_pd(&(ptrSource[i+2]));
							__m128d support_2 = _mm_load_pd(&(ptrSource[i+4]));
							__m128d support_3 = _mm_load_pd(&(ptrSource[i+6]));

							__m128d one = _mm_set1_pd(1.0);
							__m128d zero = _mm_set1_pd(0.0);

							for (size_t d = 0; d < dims; d++)
							{
								__m128d eval_0 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+0]));
								__m128d eval_1 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrTransData[(d*source_size)+i+6]));;

								__m128d level = _mm_load1_pd(&(ptrLevel[(j*dims)+d]));
								__m128d index = _mm_load1_pd(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm_mul_pd(eval_0, level);
								eval_1 = _mm_mul_pd(eval_1, level);
								eval_2 = _mm_mul_pd(eval_2, level);
								eval_3 = _mm_mul_pd(eval_3, level);

								eval_0 = _mm_sub_pd(eval_0, index);
								eval_1 = _mm_sub_pd(eval_1, index);
								eval_2 = _mm_sub_pd(eval_2, index);
								eval_3 = _mm_sub_pd(eval_3, index);

								eval_0 = _mm_and_pd(abs2MaskHybridDouble, eval_0);
								eval_1 = _mm_and_pd(abs2MaskHybridDouble, eval_1);
								eval_2 = _mm_and_pd(abs2MaskHybridDouble, eval_2);
								eval_3 = _mm_and_pd(abs2MaskHybridDouble, eval_3);

								eval_0 = _mm_sub_pd(one, eval_0);
								eval_1 = _mm_sub_pd(one, eval_1);
								eval_2 = _mm_sub_pd(one, eval_2);
								eval_3 = _mm_sub_pd(one, eval_3);

								eval_0 = _mm_max_pd(zero, eval_0);
								eval_1 = _mm_max_pd(zero, eval_1);
								eval_2 = _mm_max_pd(zero, eval_2);
								eval_3 = _mm_max_pd(zero, eval_3);

								support_0 = _mm_mul_pd(support_0, eval_0);
								support_1 = _mm_mul_pd(support_1, eval_1);
								support_2 = _mm_mul_pd(support_2, eval_2);
								support_3 = _mm_mul_pd(support_3, eval_3);
							}

							support_0 = _mm_add_pd(support_0, support_1);
							support_2 = _mm_add_pd(support_2, support_3);
							support_0 = _mm_add_pd(support_0, support_2);

							res = _mm_add_pd(res, support_0);
						}

						res = _mm_hadd_pd(res, res);

						_mm_store_sd(&(ptrGlobalResult[j]), res);
					}
				}
#else
				for (size_t j = gpu_partition; j < storageSize; j++)
				{
					#pragma omp task firstprivate(j)
					{
						ptrGlobalResult[j] = 0.0;

						for (size_t i = 0; i < source_size; i++)
						{
							double curSupport = ptrSource[i];

							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max<double>(last, 0.0);
								curSupport *= localSupport;
							}

							ptrGlobalResult[j] += curSupport;
						}
					}
				}
#endif
    		}
    	}
    }

    double time = 0.0;
    //cleanup
    delete[] ptrTransData;
#else
    double time = myOCLKernels->multOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims, storageSize);

    // do the rest...
	size_t numWGs = storageSize/OCL_MULT_N_DATAPREFETCH_BLOCKSIZE_SP;
    size_t global = numWGs*OCL_MULT_N_DATAPREFETCH_BLOCKSIZE_SP;

#ifdef USEOMP
	#pragma omp parallel for
#endif
	for (size_t j = global; j < storageSize; j++)
	{
		ptrGlobalResult[j] = 0.0;

		for (size_t i = 0; i < source_size; i++)
		{
			double curSupport = ptrSource[i];

			for (size_t d = 0; d < dims; d++)
			{
				double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
				double index_calc = eval - (ptrIndex[(j*dims)+d]);
				double abs = fabs(index_calc);
				double last = 1.0 - abs;
				double localSupport = std::max<double>(last, 0.0);
				curSupport *= localSupport;
			}

			ptrGlobalResult[j] += curSupport;
		}
	}
#endif

	return time;
}

double OperationBIterativeHybridSSEOCLLinear::multTransposeVectorized(DataVector& alpha, DataMatrix& data, DataVector& result)
{
	size_t result_size = result.getSize();
    size_t dims = storage->dim();
    size_t storageSize = storage->size();

    result.setAll(0.0);

    double* ptrAlpha = alpha.getPointer();
    double* ptrData = data.getPointer();
    double* ptrResult = result.getPointer();
    double* ptrLevel = this->Level->getPointer();
    double* ptrIndex = this->Index->getPointer();

    if (data.getNrows() % 128 != 0 || result_size != data.getNrows())
    {
    	throw operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
    }

#ifdef USEOMPTHREE
    // split result into GPU and CPU partition
    size_t cpu_partition = (result_size * PERCENT_CPUS)/100;
    size_t cpu_pad = cpu_partition % 128;
    cpu_partition -= cpu_pad;
    size_t gpu_partition = result_size - cpu_partition;

	#pragma omp parallel
    {
		#pragma omp single nowait
    	{
			#pragma omp task
    		{
    			myOCLKernels->multTransOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, gpu_partition);
    		}

			#pragma omp task
    		{
#ifdef USEICCINTRINSICS
				for (size_t i = gpu_partition; i < result_size; i+=8)
				{
					#pragma omp task firstprivate(i)
					{
						__m128d res_0 = _mm_load_pd(&(ptrResult[i+0]));
						__m128d res_1 = _mm_load_pd(&(ptrResult[i+2]));
						__m128d res_2 = _mm_load_pd(&(ptrResult[i+4]));
						__m128d res_3 = _mm_load_pd(&(ptrResult[i+6]));

						// Do on-demand transpose
						double* ptrTransData = new double[dims*8];
						for (size_t n = 0; n < 8; n++)
						{
							for(size_t d = 0; d < dims; d++)
							{
								ptrTransData[(d*8)+n] = ptrData[((i+n)*dims)+d];
							}
						}

						for (size_t j = 0; j < storageSize; j++)
						{
							__m128d support_0 = _mm_load1_pd(&(ptrAlpha[j]));
							__m128d support_1 = _mm_load1_pd(&(ptrAlpha[j]));
							__m128d support_2 = _mm_load1_pd(&(ptrAlpha[j]));
							__m128d support_3 = _mm_load1_pd(&(ptrAlpha[j]));

							__m128d one = _mm_set1_pd(1.0);
							__m128d zero = _mm_set1_pd(0.0);

							for (size_t d = 0; d < dims; d++)
							{
								__m128d eval_0 = _mm_load_pd(&(ptrTransData[(d*8)+0]));
								__m128d eval_1 = _mm_load_pd(&(ptrTransData[(d*8)+2]));
								__m128d eval_2 = _mm_load_pd(&(ptrTransData[(d*8)+4]));
								__m128d eval_3 = _mm_load_pd(&(ptrTransData[(d*8)+6]));;

								__m128d level = _mm_load1_pd(&(ptrLevel[(j*dims)+d]));
								__m128d index = _mm_load1_pd(&(ptrIndex[(j*dims)+d]));

								eval_0 = _mm_mul_pd(eval_0, level);
								eval_1 = _mm_mul_pd(eval_1, level);
								eval_2 = _mm_mul_pd(eval_2, level);
								eval_3 = _mm_mul_pd(eval_3, level);

								eval_0 = _mm_sub_pd(eval_0, index);
								eval_1 = _mm_sub_pd(eval_1, index);
								eval_2 = _mm_sub_pd(eval_2, index);
								eval_3 = _mm_sub_pd(eval_3, index);

								eval_0 = _mm_and_pd(abs2MaskHybridDouble, eval_0);
								eval_1 = _mm_and_pd(abs2MaskHybridDouble, eval_1);
								eval_2 = _mm_and_pd(abs2MaskHybridDouble, eval_2);
								eval_3 = _mm_and_pd(abs2MaskHybridDouble, eval_3);

								eval_0 = _mm_sub_pd(one, eval_0);
								eval_1 = _mm_sub_pd(one, eval_1);
								eval_2 = _mm_sub_pd(one, eval_2);
								eval_3 = _mm_sub_pd(one, eval_3);

								eval_0 = _mm_max_pd(zero, eval_0);
								eval_1 = _mm_max_pd(zero, eval_1);
								eval_2 = _mm_max_pd(zero, eval_2);
								eval_3 = _mm_max_pd(zero, eval_3);

								support_0 = _mm_mul_pd(support_0, eval_0);
								support_1 = _mm_mul_pd(support_1, eval_1);
								support_2 = _mm_mul_pd(support_2, eval_2);
								support_3 = _mm_mul_pd(support_3, eval_3);
							}

							res_0 = _mm_add_pd(res_0, support_0);
							res_1 = _mm_add_pd(res_1, support_1);
							res_2 = _mm_add_pd(res_2, support_2);
							res_3 = _mm_add_pd(res_3, support_3);
						}

						delete[] ptrTransData;

						_mm_store_pd(&(ptrResult[i+0]), res_0);
						_mm_store_pd(&(ptrResult[i+2]), res_1);
						_mm_store_pd(&(ptrResult[i+4]), res_2);
						_mm_store_pd(&(ptrResult[i+6]), res_3);
					}
				}
#else
				for (size_t i = gpu_partition; i < result_size; i++)
				{
					#pragma omp task firstprivate(i)
					{
						for (size_t j = 0; j < storageSize; j++)
						{
							double curSupport = ptrAlpha[j];

							for (size_t d = 0; d < dims; d++)
							{
								double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(i*dims)+d]));
								double index_calc = eval - (ptrIndex[(j*dims)+d]);
								double abs = fabs(index_calc);
								double last = 1.0 - abs;
								double localSupport = std::max<double>(last, 0.0);
								curSupport *= localSupport;
							}

							ptrResult[i] += curSupport;
						}
					}
				}
#endif
    		}
    	}
    }
    double time = 0.0;
#else
    double time = myOCLKernels->multTransOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims, result_size);
#endif
   	return time;
}

}