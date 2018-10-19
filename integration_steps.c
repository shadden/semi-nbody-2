#include "SemiNBody.h"
#define NDIM 4
#define INDX(ROW,COL) NDIM * ROW + COL
void update_particle_megno(Particle * particle , const double t, const double dt){
	PhaseState * state = &(particle->state);
	MEGNO_Auxilary_Variables m1 = (particle->megno);
	double delta[NDIM]  = {state->dlmbda,state->dY,state->dLambda,state->dX};
	double delta_dot[NDIM]  = {state->dlmbda_dot,state->dY_dot,state->dLambda_dot,state->dX_dot};
	double deltaSq = 0;
	double deltad_delta_dt=0 ;
	const double t1 = t + dt;
	for(int i=0;i<NDIM;i++){
		deltaSq += delta[i] * delta[i];
		deltad_delta_dt += delta[i] * delta_dot[i] * dt;
	};

	(particle->megno).Y = m1.Y + deltad_delta_dt / deltaSq * t1 ;
	(particle->megno).W = m1.W + dt * 2 * (particle->megno).Y  / t1 ;
	(particle->megno).megno = (particle->megno).W / t1;
}
void kepler_step_particle(Particle * particle,const double t, const double dt){
	PhaseState * state = &(particle->state);
	const double L = state->Lambda;
	const double dL = state->dLambda;
	const double l_dot = 1. /L/L/L;
	const double dl_dot = -3 * l_dot * dL / L;
	(state->lmbda) += l_dot * dt;
	(state->dlmbda) += dl_dot * dt;
	
}
void particles_kepler_step(Simulation * sim,const double dt){
	Particle * particle;
	const int Nparticle = sim->N_particles;
	for(int i=0;i<Nparticle;i++){
		particle = ((sim->particles) + i);
		kepler_step_particle(particle,sim->t + dt, dt);
	}
}

void kepler_step(Simulation * sim, const double dt){
	
	const int N = sim->N_planets;
	for(int i=0; i<N;i++){
		kepler_2D_advance_simple(&((sim->planets + i)->state),dt);
	}
	particles_kepler_step(sim,dt);
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
	state->lmbda = vars[0] + derivs[0] * dt;
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
	state->dlmbda = delta_vars[0] + delta_dvars[0];
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
	double delta_vars[NDIM] ={state->dlmbda,state->dY,state->dLambda,state->dX};
	for(int i=0; i<NDIM;i++){
		delta_dot[i]=0;
		for(int j=0; j<NDIM;j++){
			// Factor of 1/2 included because jac is really 2 x jac!!
			delta_dot[i] += 0.5 * *(jac + INDX(i,j)) * delta_vars[j] ;
		}
	}
}
void interaction_step(Simulation * sim,const double dt){
	// NOTE 'calloc' ensures intial values are all zero!!
	double * derivs = (double *)malloc(NDIM * sizeof(double));
	double * jac = (double *)malloc(NDIM * NDIM * sizeof(double));
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
		
		vars[0]=particle_state->lmbda;
		vars[1]=particle_state->Y;
		vars[2]=particle_state->Lambda;
		vars[3]=particle_state->X;
		delta_vars[0]=particle_state->dlmbda;
		delta_vars[1]=particle_state->dY;
		delta_vars[2]=particle_state->dLambda;
		delta_vars[3]=particle_state->dX;

		fill_zeroes(derivs,NDIM);
		fill_zeroes(jac,NDIM*NDIM);
		
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
		particle_state->lmbda = vars[0];   particle_state->dlmbda = delta_vars[0];
		particle_state->Y      = vars[1];   particle_state->dY      = delta_vars[1];
		particle_state->Lambda = vars[2];   particle_state->dLambda = delta_vars[2];
		particle_state->X      = vars[3];   particle_state->dX      = delta_vars[3];
		update_particle_state(particle_state,vars,delta_vars,derivs,jac, 0.5 * dt);

	}	

	free(derivs); 
	free(jac);
	free(planet_derivs);
	free(planet_jac);
};
void simple_pendulum_step(Simulation * sim, const double dt){
	const int Nparticle = sim->N_particles;
	// a simple additional perturbation of the form:
	//  cos(lmbda - t)
	double l,L,Delta_L,dl;
	const double t = sim->t;
	const double t1 = t + dt;
	Particle * particle;
	for(int i=0; i<Nparticle; i++){
		particle = (sim->particles + i);
		l = (particle->state).lmbda;
		L = (particle->state).Lambda;
		dl = (particle->state).dlmbda;
		Delta_L =  cos(l-t1) - cos(l-t);
		(particle->state).Lambda = L + Delta_L;
		(particle->state).dLambda += dl * (sin(l-t)-sin(l-t1));
	}
}
void get_var_dot(Simulation * sim){
	
	double L ;
	double dL;
	double L4;
	
	double * derivs = (double *)malloc(NDIM * sizeof(double));
	double * jac = (double *)malloc(NDIM * NDIM * sizeof(double));
	double * planet_derivs = (double *)malloc(NDIM * sizeof(double));
	double * planet_jac = (double *)malloc(NDIM * NDIM * sizeof(double));
	
	Particle * particle;
	PhaseState * particle_state;
	CartesianPhaseStateSimple * planet_state;
	double mu;
	const int Nparticle = sim->N_particles;
	const int Nplanet = sim->N_planets;
	double delta_dot[NDIM];
	// loop over particles
	for(int i=0; i<Nparticle;i++){

		particle = sim->particles + i;
		particle_state = &(particle->state);

		fill_zeroes(derivs,NDIM);
		fill_zeroes(jac,NDIM*NDIM);
		 
		// sum over perturbing planets	
		for(int j=0;j<Nplanet;j++ ){
			planet_state = &( ((sim->planets)+j)->state);
			mu = ((sim->planets)+j)->mu;
			interaction_derivs(particle_state,planet_state,planet_derivs,planet_jac);
			add_derivs_and_jacobians(derivs,jac,planet_derivs,planet_jac,mu);
		}	
		fill_delta_dot_from_jacobian( delta_dot, particle_state, jac) ;
		
		L = particle_state->Lambda;
		L4 = L*L*L*L;
		dL = particle_state->dLambda;

		particle_state->dlmbda_dot = delta_dot[0] - 3 * dL / L4;
		particle_state->dY_dot = delta_dot[1];
		particle_state->dLambda_dot = delta_dot[2];
		particle_state->dX_dot = delta_dot[3];
	}
	free(derivs); 
	free(jac);
	free(planet_derivs);
	free(planet_jac);
}
void simulation_step(Simulation * sim){
	const double dt = sim->dt;
	kepler_step( sim, 0.5 * dt);
	interaction_step(sim,dt);
	kepler_step( sim, 0.5 * dt);

	get_var_dot(sim);
	for(int i=0; i <sim->N_particles;i++){
		update_particle_megno((sim->particles+i), sim->t ,dt);
	}
	sim->t += dt;
}
void integrate_simulation(Simulation * sim, const double t_stop){
	const double t_now = sim->t;
	const double delta_t = t_stop - t_now;
	double dt = sim->dt;
	if(delta_t * dt < 0){
		sim->dt = -1 * dt;
		dt = -1 * dt;
	}
	const int Nstep = (int) delta_t / dt;
	for(int i=0; i<Nstep;i++){
		simulation_step(sim);
	}
}
