#include "SemiNBody.h"
#define RT2 1.414213562373095

void orbel_to_rtheta(const double a,const double lmbda, const double e, const double pomega, double * r, double * theta){
	const double M = lmbda - pomega;
	*(r) = a * (1 + e * cos(M));
	*(theta) = lmbda + 2 * e * sin(M);
}

void interaction_derivs(PhaseState * particle_state, CartesianPhaseStateSimple * planet_state, double * derivs, double * jac){
	const PhaseState p1 = *(particle_state);
	
	// ex = e * cos(-pomega)
	// ey = e * sin(-pomega)
	
	// test-particle orbit
	const double ex = p1.X ;
	const double ey = p1.Y ;
	const double lmbda = p1.lmbda ;
	const double cosl = cos(lmbda) ;
	const double sinl = sin(lmbda) ;
	const double a = p1.Lambda * p1.Lambda;
	const double r = a * (1 - ex*cosl + ey * sinl) ;
	const double theta = lmbda + 2 * sinl * ex + 2 * cosl * ey ;

	// derivatives
	const double dr_dlmbda = a * ( ex*sinl + ey * cosl );
	const double dr_dLambda = 2 * r / p1.Lambda;
	const double dr_dX = -1 * a * cosl;
	const double dr_dY =      a * sinl;
	const double dtheta_dlmbda = 1 + 2 * ex * cosl - 2 * ey * sinl;
	const double dtheta_dLambda = 0;
	const double dtheta_dX = 2 * sinl;
	const double dtheta_dY = 2 * cosl;

	const double Dr[4] = {dr_dlmbda,dr_dY,dr_dLambda,dr_dX};
	const double Dtheta[4] = {dtheta_dlmbda,dtheta_dY,dtheta_dLambda,dtheta_dX};
	// 2nd derivatives
	// r
	const double d2r_dldl = a - r;
	const double d2r_dldY = a * cosl;
	const double d2r_dldL = 2 * dr_dlmbda / p1.Lambda;
	const double d2r_dldX = a * sinl;

	const double d2r_dYdY = 0.;
	const double d2r_dYdL = 2 * dr_dY / p1.Lambda;
	const double d2r_dYdX = 0.;

	const double d2r_dLdL = dr_dLambda / p1.Lambda;
	const double d2r_dLdX = 2 * dr_dX / p1.Lambda;

	const double d2r_dXdX = 0;

	// theta
	const double d2theta_dldl = lmbda - theta;
	const double d2theta_dldY = -2 * sinl;
	const double d2theta_dldL =  0.;
	const double d2theta_dldX =  2 * cosl;

	const double d2theta_dYdY = 0.;
	const double d2theta_dYdL = 0.;
	const double d2theta_dYdX = 0.;

	const double d2theta_dLdL = 0.;
	const double d2theta_dLdX = 0.;

	const double d2theta_dXdX = 0;
	const double D2r[4][4] = {
				{d2r_dldl,d2r_dldY,d2r_dldL,d2r_dldX},
				{d2r_dldY,d2r_dYdY,d2r_dYdL,d2r_dYdX},
				{d2r_dldL,d2r_dYdL,d2r_dLdL,d2r_dLdX},
				{d2r_dldX,d2r_dYdX,d2r_dLdX,d2r_dXdX}
				
	};
	const double D2theta[4][4] = {
				{d2theta_dldl,d2theta_dldY,d2theta_dldL,d2theta_dldX},
				{d2theta_dldY,d2theta_dYdY,d2theta_dYdL,d2theta_dYdX},
				{d2theta_dldL,d2theta_dYdL,d2theta_dLdL,d2theta_dLdX},
				{d2theta_dldX,d2theta_dYdX,d2theta_dLdX,d2theta_dXdX}
				
	};
	// planet-particle separation
	const double x1 = planet_state->x ;
	const double y1 = planet_state->y ;
	const double sin_theta = sin(theta) ;
	const double cos_theta = cos(theta) ;
	const double rho1sq = x1*x1 + y1*y1 + r * r - 2 * x1 * r * cos_theta - 2 * y1 * r * sin_theta;

	const double rho1 = sqrt( rho1sq );
	const double rho1_inv3 = 1./(rho1sq * rho1);
	const double rho1_inv5 = rho1_inv3/rho1sq;
	
	const double drho1sq_dr = 2 * r - 2 * x1 * cos_theta - 2 * y1 * sin_theta;
	const double drho1sq_dtheta = 2 * r * x1 * sin_theta - 2 * r * y1 * cos_theta;
	const double drho1sq_drdr = 2;
	const double drho1sq_drdtheta = 2 * x1 * sin_theta - 2 * y1 * cos_theta ;
	const double drho1sq_dthetadtheta = 2 * r * x1 * cos_theta + 2 * r * y1 * sin_theta;

	// time derivatives
	const double J[NDIM][NDIM] = {{0,0,1,0},{0,0,0,1},{-1,0,0,0},{0,-1,0,0}};
	const double DH_rtheta[2] ={0.5 * drho1sq_dr * rho1_inv3 , 0.5 * drho1sq_dtheta * rho1_inv3};
	double D2H_rtheta[2][2];
	D2H_rtheta[0][0] = -0.75 * drho1sq_dr * drho1sq_dr * rho1_inv5 + 0.5 * drho1sq_drdr * rho1_inv3 ;
	D2H_rtheta[1][1] = -0.75 * drho1sq_dtheta * drho1sq_dtheta * rho1_inv5 + 0.5 * drho1sq_dthetadtheta * rho1_inv3 ;
	D2H_rtheta[0][1] = -0.75 * drho1sq_dtheta * drho1sq_dr * rho1_inv5 + 0.5 * drho1sq_drdtheta * rho1_inv3 ;
	D2H_rtheta[1][0] = D2H_rtheta[0][1];
	
	double d2H[NDIM][NDIM];

	for (int i=0; i<NDIM; i++){
		*(derivs+i)=0;
		for (int j=0; j<NDIM; j++){
			*(derivs+i) += J[i][j] * ( DH_rtheta[0] * Dr[j] + DH_rtheta[1] * Dtheta[j]);
			d2H[i][j]  = D2H_rtheta[0][0] * Dr[i] * Dr[j] + D2H_rtheta[1][1] * Dtheta[i] * Dtheta[j];
			d2H[i][j] += D2H_rtheta[0][1] * Dr[i] * Dtheta[j] + D2H_rtheta[1][0] * Dr[j] * Dtheta[i];
			d2H[i][j] += DH_rtheta[0] * D2r[i][j] + DH_rtheta[1] * D2theta[i][j];
		}
	}
	for(int i=0;i<NDIM;i++){
		for(int j=0;j<NDIM;j++){
			*(jac+INDX(i,j)) = 0;
			for(int k=0;k<4;k++){
				*(jac+INDX(i,j)) += *(*(J+i)+k) * *(*(d2H+k)+j);
			}
		}
	}	
}
