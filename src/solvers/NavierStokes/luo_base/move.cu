/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/luo_base.h>

#include <solvers/NavierStokes/luo_base/kernels/structure.h> //update_body_viv

void luo_base::set_movement()
{
	double	t = dt*timeStep,
			f = B.xfrequency,
			xCoeff = B.xCoeff,
			uCoeff = B.uCoeff,
			xPhase = B.xPhase,
			uPhase = B.uPhase,
			totalPoints=B.totalPoints,
			unew,
			xnew;

	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);

	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	//B.uBk = B.uB;
	kernels::update_body_viv<<<grid,block>>>(B.x_r, B.uB_r, B.dx_r, unew, B.midX, totalPoints);
}

void luo_base::viv_movement_LC()
{
	double	Cy	= B.forceY*2.0,
			Mred= 2.0,
			Ured= (*paramDB)["simulation"]["Ured"].get<int>(),
			totalPoints=B.totalPoints,
			vold = B.centerVelocityV0,
			yold = B.midY0,
			vnew,
			ynew;

	double a = dt*M_PI*M_PI*4/(Ured*Ured),
		   b = dt*dt*2*M_PI*M_PI/(Ured*Ured);

	//calc updated velocity
	vnew = (vold - a*(yold+ dt/2*vold) + dt*Cy/2/Mred)/(1+b);
	ynew = yold + dt/2*(vnew + vold);
	B.centerVelocityV = vnew;
	B.midY = ynew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(B.y_r, B.vB_r, B.dy_r, vnew, B.midY, totalPoints);
}

void luo_base::viv_movement_SC()
{
	double	Cy	= B.forceY*2.0,
			Mred= 2.0,
			Ured= (*paramDB)["simulation"]["Ured"].get<int>(),
			totalPoints=B.totalPoints,
			vold= B.centerVelocityV0,
			yold= B.midY0,
			vnew,
			ynew,
			relax = 0.9; //this should be set somewhere else
	SCtol = B.midY;

	double a = dt*M_PI*M_PI*4/(Ured*Ured),
		   b = dt*dt*2*M_PI*M_PI/(Ured*Ured);

	//calc updated velocity
	vnew = (vold - a*(yold+ dt/2*vold) + dt*Cy/2/Mred)/(1+b);
	ynew = relax*(yold + dt/2*(vnew + vold)) + (1-relax)*B.midY;
	B.centerVelocityV = vnew;
	B.midY = ynew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(B.y_r, B.vB_r, B.dy_r, vnew, B.midY, totalPoints);
	//std::cout<<timeStep<<"\t"<<B.midY<<"\t"<<cfl_max<<"\n";

	SCtol = (SCtol-ynew)/SCtol;
}
