#include "SemiNBody.h"	

static inline double random_real(){
        return -1+2.*((float)rand())/RAND_MAX;
}

void initialize_megno_vars(MEGNO_Auxilary_Variables * megno){
	megno->Y=0;
	megno->W=0;
	megno->megno=2;
}

void initialize_phase_state(PhaseState * Z,double a, double l, double e, double pomega){
	Z->Lambda=sqrt(a);
	Z->lambda=l;
	Z->X = e * cos(pomega);
	Z->Y = e * sin(-1 * pomega);
	// generate random reals	
	double dvec[4];
	double normsq=0;
	for(int i=0; i < 4; i++){
		dvec[i] = random_real();
		normsq += dvec[i]*dvec[i];
	}
	const double EPS = 1.e-14;
	Z->dLambda = EPS*dvec[0]/sqrt(normsq);
	Z->dlambda = EPS*dvec[1]/sqrt(normsq);
	Z->dX = EPS*dvec[2]/sqrt(normsq);
	Z->dY = EPS*dvec[3]/sqrt(normsq);
}
