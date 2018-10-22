#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define NDIM 4 
#define INDX(ROW,COL) NDIM * ROW + COL



typedef struct PhaseState {
  double lmbda,Y,Lambda,X,dlmbda,dY,dLambda,dX;
  double dlmbda_dot,dY_dot,dLambda_dot,dX_dot;
} PhaseState;

typedef struct CartesianPhaseStateSimple {
  double x,y,vx,vy;
} CartesianPhaseStateSimple;

typedef struct CartesianPhaseState {
  double x,y,vx,vy,dx,dy,dvx,dvy;
} CartesianPhaseState;

typedef struct MEGNO_Auxilary_Variables {
  double W, Y, megno;
} MEGNO_Auxilary_Variables;

typedef struct Planet {
	double mu;
	CartesianPhaseStateSimple state;
} Planet;
typedef struct Particle {
	PhaseState state;
	MEGNO_Auxilary_Variables megno;
} Particle;

typedef struct Resonance {
	bool inner_perturber_Q;
	int j,k,l;
	double fCoeff;
} Resonance;

typedef struct ResonancePerturbation {
	double e,lmbda,pomega;
	double mean_motion;
	double mu;
	int N_resonances;
	Resonance * resonances;
	
} ResonancePerturbation;

typedef struct Simulation {
	int N_planets;
	int N_resonance;
	int N_particles;
	Planet * planets;
	ResonancePerturbation * resonances;
	Particle * particles;
	double t,dt;
} Simulation;
double GeneralOrderCoefficient(int , int , int ,double);
void fill_zeroes(double* ,const int );
void initialize_megno_vars(MEGNO_Auxilary_Variables *);
void initialize_phase_state(PhaseState * ,double , double , double , double );
void initialize_resonance_multiplet(Resonance *, bool, int,int, double);
void interaction_step(Simulation *, const double);
void simulation_step(Simulation *);
void free_simulation( Simulation *);
void kepler_2D_advance(  CartesianPhaseState*  ,double);
void kepler_2D_advance_simple(  CartesianPhaseStateSimple*  ,double);
void interaction_derivs(PhaseState *, CartesianPhaseStateSimple * ,double * ,double *);
void add_resonance_derivs(PhaseState* , ResonancePerturbation* ,double* ,double* );
void add_derivs_and_jacobians(double * , double *, double *, double *, double);
int binomialCoeff(int n, int k);
//
//

