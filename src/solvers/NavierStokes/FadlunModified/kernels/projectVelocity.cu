/***************************************************************************//**
 * \file projectVelocity.cu
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to update the velocity field
 */


#include "projectVelocity.h"

namespace kernels
{
__global__
void project_velocity_X(double *u, double *uhat, double *uold, double *pressure, int *ghostTagsP, double *dx, double dt, int nx, int ny)
{
	int i	= threadIdx.x + (blockDim.x * blockIdx.x),
		I	= i % (nx-1),
		J 	= i / (nx-1),
		ip  = J*nx + I,
		numU = (nx-1)*ny;

	if (i >= numU)
		return;

	u[i] = uhat[i] - (ghostTagsP[ip+1] ==-1 && ghostTagsP[ip] == -1) * dt*(pressure[ip+1]-pressure[ip]) / (0.5*dx[I+1]+0.5*dx[I]);
}

__global__
void project_velocity_Y(double *u, double *uhat, double *uold, double *pressure, int *ghostTagsP, double *dy, double dt, int nx, int ny)
{
	int numU= (nx-1)*ny,
		i	= threadIdx.x + (blockDim.x * blockIdx.x),
		I	= i % nx,
		J	= i / nx,
		ip	= J*nx + I,
		numUV	= (ny-1)*nx + numU;

	i += numU;

	if (i >= numUV)
		return;

	u[i] = uhat[i] - (ghostTagsP[ip+nx] == -1 && ghostTagsP[ip] == -1) * dt*(pressure[ip+nx]-pressure[ip]) / (0.5*dy[J+1]+0.5*dy[J]);
}

}//end namespace kernels
