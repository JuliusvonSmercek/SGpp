/*****************************************************************************/
/* This file is part of jsgpp, a program package making use of spatially     */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)                         */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* jsgpp is free software; you can redistribute it and/or modify             */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* jsgpp is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with jsgpp; if not, write to the Free Software                      */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

%newobject sg::Grid::createLinearGrid(size_t dim);
%newobject sg::Grid::createLinearBoundaryGrid(size_t dim);
%newobject sg::Grid::createLinearTrapezoidBoundaryGrid(size_t dim);
%newobject sg::Grid::createLinearTrapezoidBoundaryGrid(BoudingBox& BB);
%newobject sg::Grid::createModLinearGrid(size_t dim);
%newobject sg::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject sg::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject sg::Grid::createModWaveletGrid(size_t dim);
%newobject sg::Grid::createModBsplineGrid(size_t dim, size_t degree);

%newobject sg::Grid::unserialize(std::string& istr);

%newobject sg::Grid::createOperationB();
%newobject sg::Grid::createGridGenerator();
%newobject sg::Grid::createOperationLaplace();
%newobject sg::Grid::createOperationEval();
%newobject sg::Grid::createOperationTest();
%newobject sg::Grid::createOperationHierarchisation();

%include "stl.i"
%include "typemaps.i"

%apply std::string *OUTPUT { std::string& ostr };
%apply std::string *INPUT { std::string& istr };



//void getMemento();


namespace sg
{

class Grid
{
public:
	static Grid* createLinearGrid(size_t dim);
	static Grid* createLinearBoundaryGrid(size_t dim);
	static Grid* createLinearTrapezoidBoundaryGrid(size_t dim);
	static Grid* createModLinearGrid(size_t dim);
	static Grid* createPolyGrid(size_t dim, size_t degree);
	static Grid* createModPolyGrid(size_t dim, size_t degree);
	static Grid* createModWaveletGrid(size_t dim);
	static Grid* createModBsplineGrid(size_t dim, size_t degree);
	
	static Grid* unserialize(std::string& istr);
	
protected:
	Grid();
	Grid(Grid& o);

public:
	virtual ~Grid();

public:	
	virtual GridGenerator* createGridGenerator() = 0;
	virtual OperationB* createOperationB() = 0;
	virtual OperationEval* createOperationEval() = 0;
	virtual OperationTest* createOperationTest() = 0;
	virtual OperationMatrix* createOperationLaplace() = 0;
	virtual OperationMatrix* createOperationIdentity() = 0;
	virtual OperationHierarchisation* createOperationHierarchisation() = 0;
	
	// @todo remove this when done
	virtual OperationMatrix* createOperationUpDownTest() = 0;
	
	virtual OperationMatrix* createOperationDelta(DataVector& coef) = 0;
	virtual OperationMatrix* createOperationGamma(DataVector& coef) = 0;
	
	virtual GridStorage* getStorage();
	virtual BoundingBox* getBoundingBox();

	virtual const char* getType() = 0;	
	virtual void serialize(std::string& ostr);
	void refine(DataVector* vector, int num);
	double eval(DataVector& alpha, DataVector& point);
	void insertPoint(size_t dim, unsigned int levels[], unsigned int indeces[], bool isLeaf);
	int getSize();
	
};



}

//these are just two new interfaces for consistency with Memento design pattern
%extend sg::Grid{
	Grid* createMemento(){
		return $self;
	}
	
	static sg::Grid* setMemento(std::string& istr){
		return sg::Grid::unserialize(istr);
	}
};

	