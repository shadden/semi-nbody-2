#include "SemiNBody.h"
#define NDIM 4
#define INDX(ROW,COL) NDIM * ROW + COL
void update_particle_megno(Particle * particle,double * delta_dot , const double t, const double dt){
	PhaseState * state = &(particle->state);
	MEGNO_Auxilary_Variables m1 = (particle->megno);
	double delta[NDIM]  = {state->dlambda,state->dY,state->dLambda,state->dX};
	double deltaSq = 0;
	double deltad_delta=0;
	for(int i=0;i<NDIM;i++){
		deltaSq += delta[i] * delta[i];
		deltad_delta += delta[i] * delta_dot[i];
	};
	(particle->megno).Y = m1.Y + deltad_delta / deltaSq * t * dt;
	(particle->megno).W = m1.W + dt * 2 * m1.Y / t;
	(particle->megno).megno = m1.W / t;
}
void kepler_step_particle(Particle * particle,const double t, const double dt){
	PhaseState * state = &(particle->state);
	const double L = state->Lambda;
	const double dL = state->dLambda;
	const double l_dot = 1. /L/L/L;
	const double dl_dot = -3 * l_dot * dL / L;
	(state->lambda) += l_dot * dt;
	(state->dlambda) += dl_dot * dt;
	double delta_dot[4] = {dl_dot,0,0,0};
	update_particle_megno(particle,delta_dot,t,dt);
	
}
void particles_kepler_step(Simulation * sim,const double dt){
	Particle * particle;
	const int Nparticle = sim->N_particles;
	for(int i=0;i<Nparticle;i++){
		particle = ((sim->particles) + i);
		kepler_step_particle(particle,sim->t + dt, dt);
	}
}

void kepler_step(Simulation * sim, const double _dt){
	const int N = sim->N_planets;
	for(int i=0; i<N;i++){
		kepler_2D_advance_simple(&((sim->planets + i)->state),_dt);
	}
	particles_kepler_step(sim , _dt);
}



void add_derivs_and_jacobians(double * tot_deriv, double * tot_jac, double * deriv, double * jac, double mu){
	for(int i=0;i<NDIM;i++){
		*(tot_deriv+i) += mu * *(deriv+i);
		for(int j=0; j<NDIM;j++){
			*(tot_jac +INDX(i,j))+= mu * *(jac +INDX(i,j));
		}
	}
}
void update_particle_state(PhaseState * state,double* vars,double* delta_vars, double* derivs, double* jac, const double dt){
	state->lambda = vars[0] + derivs[0] * dt;
	state->Y      = vars[1] + derivs[1] * dt;
	state->Lambda = vars[2] + derivs[2] * dt;
	state->X      = vars[3] + derivs[3] * dt;
	double delta_dvars[NDIM];
	for(int i=0; i<NDIM;i++){
		delta_dvars[i]=0;
		for(int j=0; j<NDIM;j++){
			delta_dvars[i] += *(jac + INDX(i,j)) * delta_vars[j] * dt;
		}
	}
	state->dlambda = delta_vars[0] + delta_dvars[0];
	state->dY      = delta_vars[1] + delta_dvars[1];
	state->dLambda = delta_vars[2] + delta_dvars[2];
	state->dX      = delta_vars[3] + delta_dvars[3];
}
void fill_zeroes(double* arr,const int N){
	for(int i=0;i<N;i++){
		*(arr+i)=0;
	}
}
void fill_delta_dot_from_jacobian(double * delta_dot, PhaseState * state, double * jac){
	double delta_vars[NDIM] ={state->dlambda,state->dY,state->dLambda,state->dX};
	for(int i=0; i<NDIM;i++){
		delta_dot[i]=0;
		for(int j=0; j<NDIM;j++){
			// Factor of 1/2 included because jac is really 2 x jac!!
			delta_dot[i] += 0.5 * *(jac + INDX(i,j)) * delta_vars[j] ;
		}
	}
}
void interaction_step(Simulation * sim, const double dt){
	// NOTE 'calloc' ensures intial values are all zero!!
	double * derivs = (double *)malloc(NDIM * sizeof(double));
	double * jac = (double *)malloc(NDIM * NDIM * sizeof(double));
	fill_zeroes(derivs,NDIM);
	fill_zeroes(jac,NDIM*NDIM);
	double * planet_derivs = (double *)malloc(NDIM * sizeof(double));
	double * planet_jac = (double *)malloc(NDIM * NDIM * sizeof(double));
	
	Particle * particle;
	PhaseState * particle_state;
	CartesianPhaseStateSimple * planet_state;
	double mu;
	const int Nparticle = sim->N_particles;
	const int Nplanet = sim->N_planets;
	double vars[NDIM],delta_vars[NDIM];
	// loop over particles
	
	for(int i=0; i<Nparticle;i++){
		// sum over perturbing planets
		particle = sim->particles + i;
		particle_state = &(particle->state);
		
		vars[0]=particle_state->lambda;
		vars[1]=particle_state->Y;
		vars[2]=particle_state->Lambda;
		vars[3]=particle_state->X;
		delta_vars[0]=particle_state->dlambda;
		delta_vars[1]=particle_state->dY;
		delta_vars[2]=particle_state->dLambda;
		delta_vars[3]=particle_state->dX;

		for(int j=0;j<Nplanet;j++ ){
			planet_state = &( ((sim->planets)+j)->state);
			mu = ((sim->planets)+j)->mu;
			interaction_derivs(particle_state,planet_state,planet_derivs,planet_jac);
			add_derivs_and_jacobians(derivs,jac,planet_derivs,planet_jac,mu);
		}	

		update_particle_state(particle_state,vars,delta_vars,derivs,jac,dt);

		for(int j=0;j<Nplanet;j++ ){
			planet_state = &( ((sim->planets)+j)->state);
			mu = ((sim->planets)+j)->mu;
			interaction_derivs(particle_state,planet_state,planet_derivs,planet_jac);
			add_derivs_and_jacobians(derivs,jac,planet_derivs,planet_jac,mu);
		}	
		// Reset particle and then update with average of two derivs added together
		particle_state->lambda = vars[0];   particle_state->dlambda = delta_vars[0];
		particle_state->Y      = vars[1];   particle_state->dY      = delta_vars[1];
		particle_state->Lambda = vars[2];   particle_state->dLambda = delta_vars[2];
		particle_state->X      = vars[3];   particle_state->dX      = delta_vars[3];
		update_particle_state(particle_state,vars,delta_vars,derivs,jac, 0.5 * dt);

		double delta_dot[NDIM];
		fill_delta_dot_from_jacobian(delta_dot,particle_state,jac);
		update_particle_megno(particle, delta_dot, sim->t + dt, dt);
	}	

	free(derivs); 
	free(jac);
	free(planet_derivs);
	free(planet_jac);
};

void simulation_step(Simulation * sim, const double dt){
	kepler_step( sim, 0.5 * dt);
	interaction_step( sim, dt);
	kepler_step( sim, 0.5 * dt);
	sim->t += dt;
}
