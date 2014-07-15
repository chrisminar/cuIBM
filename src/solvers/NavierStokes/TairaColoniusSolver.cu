/***************************************************************************//**
* \file TairaColoniusSolver.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the methods of the class TairaColoniusSolver
*/

#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <sys/stat.h>

/********************//**
* \brief Initialize the solver and assemble the matrices
*/
template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialise()
{	
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	int totalPoints  = NSWithBody<memoryType>::B.totalPoints; 
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP+2*totalPoints);
	if(totalPoints > 0)
	{
		// for each Lagragian point, 24 velocity values could be used
		// when caluclating the discrete delta function
		// 12*2 in 2D
		E.resize(2*totalPoints, numUV, 24*totalPoints);
	}
	NavierStokesSolver<memoryType>::assembleMatrices();
}

/********************//**
* \brief Write the time and forces for each body in the flow fieldi
*
* The forces are calculated using the method of Taira and Colonius (2007).
*/
template <typename memoryType>
void TairaColoniusSolver<memoryType>::writeData()
{	
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	NSWithBody<memoryType>::writeCommon();

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	real         dt  = db["simulation"]["dt"].get<real>();
	int          numBodies  = NSWithBody<memoryType>::B.numBodies;

	// calculate forces using the T&C method
	calculateForce();
	
	// print to file
	NSWithBody<memoryType>::forceFile << NavierStokesSolver<memoryType>::timeStep*dt << '\t';
	for(int l=0; l<numBodies; l++)
	{
		NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::B.forceX[l] << '\t' << NSWithBody<memoryType>::B.forceY[l] << '\t';
	}
	NSWithBody<memoryType>::forceFile << std::endl;
	
	// print forces calculated using the CV approach
	//NSWithBody<memoryType>::calculateForce();
	//NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;
	
	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

/********************//**
* \brief Update the location of the bodies and re-compute matrices
*/
template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateSolverState()
{
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		
		NSWithBody<memoryType>::updateBodies();
		updateQT();
		NavierStokesSolver<memoryType>::generateC();
		
		NavierStokesSolver<memoryType>::logger.startTimer("preconditioner2");
		NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C);
		NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner2");
	}
}

/********************//**
* \brief Constructor of the class TairaColoniusSolver.
*
* Get the simulation parameters of the database.
* Get the computational grid.
*/
template <typename memoryType>
TairaColoniusSolver<memoryType>::TairaColoniusSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

// include inline files located in "./TairaColonius/"
#include "TairaColonius/generateQT.inl"
#include "TairaColonius/generateBC2.inl"
#include "TairaColonius/calculateForce.inl"

// specialization of the class TairaColoniusSolver
template class TairaColoniusSolver<host_memory>;
template class TairaColoniusSolver<device_memory>;
