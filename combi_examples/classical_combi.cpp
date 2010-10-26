/*
 * clasical_combi.cpp
 *
 *  Created on: Aug 20, 2010
 *      Author: Aliz Nagy
 *
 *      A program presenting some examples with the basic functions introduced with the combination technique and their usage
 *      The concepts exemplified in this piece of code include the construction of a sparse grid(with the create***Grid method of the grid class),
 *      generation of a regular sparse grid, assignment of function values to the gridpoints, the construction of a set of fullgrids coreesponding to
 *      the sparse grid, the decomposition of the sparse grid into the fullgrids(assignment of the right function values to the fullgrid points),
 *      the evaluation of the function on the fullgrids, the result interpolated using the combined values on the grids, the recomposition of the fullgrids
 *      into a sparse grid(also with the combinational formula), and comparing the result obtained on full grids with the result obtained on a sparse grid
 *      and the actual function value
 */

#include <iostream>
#include <math.h>
#include "sgpp.hpp"
using namespace sg;
const double PI = 3.14159265;
/**
 * The function we will try to interpolate
 * */
double f(double x0, double x1, double x2)
{
	return 1.0+(0.25*(x0-0.7)*(x0-0.7)+2.0)+(0.25*(x1-0.7)*(x1-0.7)+2.0)+(0.25*(x2-0.7)*(x2-0.7)+2.0);
}
int main()
{
	int dim=3;
	int level=4;
	/***Create a new Grid(LinearGrid, LinearBoundary,or LinearTrapezoidBoundary)*/
	Grid* grid = Grid::createLinearTrapezoidBoundaryGrid(dim);
    GridStorage* gridStorage = grid->getStorage();
    GridGenerator* gridGen = grid->createGridGenerator();
    /***Generate a new regular grid structure corresponding to the type of the chosen grid*/
    gridGen->regular(level);
    cout<<"Dimensionality: "<<gridStorage->dim()<<endl;
    cout<<"Number of grid points: "<<gridStorage->size()<<endl;
    /***Create coefficient vector of type DataVector*/
    DataVector alpha(gridStorage->size());
    /***Set the function values for every member of the DataVector*/
    for (int i=0; i < gridStorage->size(); i++) {
    		GridIndex* gp = gridStorage->get(i);
    		alpha[i] = f(gp->abs(0),gp->abs(1),gp->abs(2));
    	 }
    /**
     * Create the set of fullgrids into which the grid decomposes, the gridType can be "linear", "linearBoundary" or "linearTrapezoidBoundary
     * */
    FullGridSet fgs(dim,level,grid->getType());
    cout<<"Number of full grids:"<<fgs.getSize()<<endl;
    cout<<"The fullgrids:"<<endl;
    /**
     * This prints the levels of all fullgrids
     */
    fgs.printGridSet();
    /**
     * Assign the function values to every gridpoint of the fullgrids, using the values from the sparse grid
     * */
    /**
     * This step can be replaced by the manual assignment of function values to every gridpoint in the fullgrids
     */
    fgs.deCompose(gridStorage,alpha);
    /**
     * Create a new datavector which contains the coordinates of the point we want to interpolate
     */
    DataVector p(dim);
    p[0] = 0.91;
    p[1] = 0.23;
    p[2] = 0.76;
    /**
     *We now verify if the decomposition of the sparse grid was correct,i.e. see if every gridpoint of the fullgrids contains the correct function value,
     *and evaluate the function
     */
    for (int i=0;i<fgs.getSize();i++)
           {
			size_t m=fgs[i].getSize();
           	FullGrid *fg=fgs.at(i);//works also with &(fgs[i])
           	for (int j=0;j<m;j++)
           	       	{
					 double d=f(fg->getCoord(0,j),fg->getCoord(1,j),fg->getCoord(2,j));
					 /**
					  * If the function value in the fullgrid point is different from the value already assigned we will signal an error(if this happens we did something wrong)
					  * */
					 if (d!=(*fg)[j]) cout<<"error: "<<i<<","<<j<<": "<<d<<"!="<<(*fg)[j]<<"\n";
           	       	}
           	/**
           	 * *Evaluates the fullgrid in an arbitrary point, and assigns the resulting value to the field variable value of the fullgrid,
           	 *  The same value is returned by the function and can be accesed later through the function call fullgrid.val()
             */
			fg->eval(p);
           }
    /**
     * Combines the interpolated results on the fullgrids into one value, which in case of function interpolation equals the value on the sparse grid itself
     * We print the interpolation value as Uc(p)
     */
    cout<<"Uc("<<p[0]<<","<<p[1]<<","<<p[2]<<")="<<fgs.combinedResult()<<endl;
    /**
     * We will now verify if the same value is obtained on the sparse grid
     * */
    DataVector beta(gridStorage->size());
    /**
     * We will recompose the values of the gridpoints into the sparse grid structure using the same coefficients as with the combinedResult function
     * This step can also be replaced by the direct assignment of function values to the sparse grid points as we showed it before
     * */
    fgs.reCompose(gridStorage,&beta);
    OperationHierarchisation* oh=grid->createOperationHierarchisation();
    oh->doHierarchisation(beta);
    OperationEval* opEval = grid->createOperationEval();
    /**
     * We print the value interpolated on the sparse grid and also the real value of the function in the given point
     * */
    cout<<"Usg("<<p[0]<<","<<p[1]<<","<<p[2]<<")="<< opEval->eval(beta, p) << endl;
    cout<<"f("<<p[0]<<","<<p[1]<<","<<p[2]<<")="<<f(p[0],p[1],p[2])<<endl;
    delete grid;
    delete gridGen;
    delete opEval;
    delete oh;
    return 0;
}