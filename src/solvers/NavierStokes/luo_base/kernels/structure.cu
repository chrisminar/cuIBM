/***************************************************************************//**
 * \file structure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 */


#include "structure.h"

namespace kernels
{

/*
 * Updates all the velocities and positions of the body nodes
 * param double y y positions of the nodes
 * param double vB v velocities of body nodes
 * param double dy change in y location
 * param vnew new v velocity
 * param totalPoints total number of body points
 */
__global__
void update_body_viv(double *By, double *vB, double *Bdy, double vnew, double midY, int totalPoints)
{
	int i	= threadIdx.x + (blockDim.x * blockIdx.x);
	if (i > totalPoints)
		return;
	vB[i] = vnew;
	By[i] = midY + Bdy[i];
}

__global__
void initialise_old(double *uB0, double unew, int totalPoints)
{
	int i	= threadIdx.x + (blockDim.x * blockIdx.x);
	if (i > totalPoints)
		return;
	uB0[i] = unew;
}

}
