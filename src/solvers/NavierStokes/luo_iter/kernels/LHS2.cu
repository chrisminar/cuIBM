/***************************************************************************//**
 * \file LHS2.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the poission solve
 */

#include "LHS2.h"

namespace kernels
{
__global__
void LHS2_mid_iter(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt,
					int *count, double *ns_rhs, double *interp_rhs, int *hybridTagsP, int *ghostTagsP,
					double *alpha, double *dpdn,
					int *index1, int *index2, int *index3, int *index4,
					double *q1coef, double *q2coef, double *q3coef, double *q4coef,
					double *q1, double *q2, double *q3, double *q4, int timeStep)
{
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x;
	if (ip >= nx*ny)
		return;
	int	I	= ip % nx,
		J	= ip / nx;

	if (I == 0 || I == nx-1 || J == 0 || J == ny-1)
		return;

	int numE = nx*4-2 + (J-1)*(nx*5-2) + I*5-1;
	double temp = 0;


	if (hybridTagsP[ip]>0)//if were at hybrid node
	{
		int interp_index[4] = {index1[ip], index2[ip], index3[ip], index4[ip]};
		//int nx_index[5] = {ip+nx, ip+1, ip-nx, ip-1, ip};//n e s w p
		double CInterp[4];
		double Cns[5];
		double q[4] = {q1[ip], q2[ip], q3[ip], q4[ip]};

		//calculate the pressure coefficients for the stencil pressure calculation
		Cns[0] = -dt/(dy[J]*(dy[J]+dy[J+1])*0.5); //N
		Cns[1] = -dt/(dx[I]*(dx[I]+dx[I+1])*0.5); //E
		Cns[2] = -dt/(dy[J]*(dy[J]+dy[J-1])*0.5); //s
		Cns[3] = -dt/(dx[I]*(dx[I]+dx[I-1])*0.5); //w
		Cns[4] = -Cns[0] - Cns[1] - Cns[2] - Cns[3]; //P

		//set pressure coefficients for the interpolation pressure calculation
		CInterp[0] = q1coef[ip];
		CInterp[1] = q2coef[ip];
		CInterp[2] = q3coef[ip];
		CInterp[3] = q4coef[ip];

		//multiply by alpha
		for (int i=0;i<4;i++)
		{
			Cns[i] = (1-alpha[ip])*Cns[i]/Cns[4];
			CInterp[i] = alpha[ip]*CInterp[i];
		}

		/*   0  1  2		NW  N   NE
		 *   3  4  5		W   P   E
		 *   6  7  8		SW  S   SE
		 */
		int stencil_index[9]    = {ip + nx - 1, ip + nx, ip + nx + 1,
								   ip - 1     , ip     , ip + 1,
								   ip - nx - 1, ip - nx, ip - nx + 1};
		double stencil[9] = {0, Cns[0], 0, Cns[3], 1, Cns[1], 0, Cns[2], 0};
		//combine ns and interp stencils
		bool stencil_used[9] = {false, true, false, true, true, true, false, true, false};
		for (int n=0;n<4;n++)
		{
			for (int m=0;m<9;m++)
			{
				if (stencil_index[m] == interp_index[n] && m != 4)
				{
					stencil[m] -= CInterp[n]; //this should be minus
				}
			}
		}
		//add ns to sparse matrix
		for (int m = 0; m<9; m++)
		{
			if (stencil_used[m])
			{
				row[numE] = ip;
				col[numE] = stencil_index[m];
				val[numE] = stencil[m];
				numE++;
			}
		}
		ns_rhs[ip] = (1-alpha[ip])/Cns[4];
		interp_rhs[ip] = 0;
		//calc new numE
		numE = ny*nx*5 - ny*2 - nx*2 + count[ip]-1;
		//add interp corner to sparse matrix
		for (int n=0;n<4;n++)
		{
			for (int m=0;m<9;m++)
			{
				if (stencil_index[m] == interp_index[n] && !stencil_used[m])
				{
					row[numE] = ip;
					col[numE] = interp_index[n];
					val[numE] = -CInterp[n];//this should also be minus
				}
				else if(stencil_index[m] == ip && interp_index[n] == ip)
					interp_rhs[ip] += CInterp[n]*q[n]; //this should be addition

			}
		}
	}
	else if (ghostTagsP[ip]>0)
	{
		int interp_index[4] = {index1[ip], index2[ip], index3[ip], index4[ip]};
		bool interp_in[4] = {false, false, false, false};
		int ns_index[5] = {ip + nx, ip + 1, ip - nx, ip -1, ip}; //n e s w p
		bool ns_overlap[5] = {false, false, false, false, true};
		double q[4] = {q1[ip], q2[ip], q3[ip], q4[ip]};
		double CInterp[4];
		CInterp[0] = q1coef[ip];
		CInterp[1] = q2coef[ip];
		CInterp[2] = q3coef[ip];
		CInterp[3] = q4coef[ip];
		//count the number of nodes the interp is using
		//find how which ns nodes are occupied
		int counter = 0;
		temp = 0;
		for (int l=0; l<4; l++)
		{
			if (ghostTagsP[interp_index[l]]>0)
			{
				counter +=1;
				interp_in[l] = true;
			}
			for (int n=0; n<5; n++)
			{
				if (interp_index[l] == ns_index[n])
					ns_overlap[n] = true;
			}
		}
		//add center to matrix
		row[numE] = ip;
		col[numE] = ip;
		val[numE] = 1;
		numE++;
		//add real interp values to matrix
		for (int i=0; i<4; i++)
		{
			if (!interp_in[i] && interp_index[i] != ip)
			{
				row[numE] = ip;
				col[numE] = interp_index[i];
				val[numE] = -CInterp[i];
				numE++;
			}
			else
			{
				temp += CInterp[i] * q[i];
			}
		}
		//fill remainder of values
		int counter2 = 0;
		for (int i=0; i<5; i++)
		{
			if (counter2>=counter)
				break;
			if (ns_overlap[i]==false)
			{
				row[numE] = ip;
				col[numE] = ns_index[i];
				val[numE] = 0;
				numE++;
				counter2++;
			}
		}
		ns_rhs[ip] = 0;
		interp_rhs[ip] = dpdn[ip] + temp;
	}
	else //if were not at a hybrid node
	{
		//temp = dt/(dx[I]*(dx[I]+dx[I+1])*0.5) + dt/(dx[I]*(dx[I]+dx[I-1])*0.5) + dt/(dy[J]*(dy[J]+dy[J+1])*0.5) + dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
		temp = 0;
		//EAST
		row[numE] = ip;
		col[numE] = ip + 1;
		//val[numE] = -dt/(dx[I]*(dx[I]+dx[I+1])*0.5)/temp;
		val[numE] = -dt/(dx[I]*(dx[I]+dx[I+1])*0.5);
		temp += dt/(dx[I]*(dx[I]+dx[I+1])*0.5);
		numE++;

		//WEST
		row[numE] = ip;
		col[numE] = ip - 1;
		//val[numE] = -dt/(dx[I]*(dx[I]+dx[I-1])*0.5)/temp;
		val[numE] = -dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
		temp +=dt/(dx[I]*(dx[I]+dx[I-1])*0.5);
		numE++;

		//NORTH
		row[numE] = ip;
		col[numE] = ip + nx;
		//val[numE] = -dt/(dy[J]*(dy[J]+dy[J+1])*0.5)/temp;
		val[numE] = -dt/(dy[J]*(dy[J]+dy[J+1])*0.5);
		temp += dt/(dy[J]*(dy[J]+dy[J+1])*0.5);
		numE++;

		//SOUTH
		row[numE] = ip;
		col[numE] = ip - nx;
		//val[numE] = -dt/(dy[J]*(dy[J]+dy[J-1])*0.5)/temp;
		val[numE] = -dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
		temp += dt/(dy[J]*(dy[J]+dy[J-1])*0.5);
		numE++;
		//MID
		row[numE] = ip;
		col[numE] = ip;
		//val[numE] = 1;
		val[numE] = temp;
		ns_rhs[ip] = 1;
		interp_rhs[ip] = 0;
		//ns_rhs[ip] = 1/temp;
		//interp_rhs[ip] = 0;
	}
}
}
