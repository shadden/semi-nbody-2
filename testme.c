#include "SemiNBody.h"
void print_particle(Particle * particle){
	double vals[4]= {particle->state.lambda,particle->state.Y,particle->state.Lambda,particle->state.X};
	for(int i =0; i<4;i++){	
		printf("%.8f\t",vals[i]);
	}
}
double get_delta_sq(Particle * particle){
	double delta[4]= {particle->state.dlambda,particle->state.dY,particle->state.dLambda,particle->state.dX};
	double delta_sq = 0;
	for(int i =0; i<4;i++){	
		delta_sq += delta[i] * delta[i];
	}
	return delta_sq;
}
double get_LCE(Particle * particle,const double log_d0, const double time){
	const double dsq = get_delta_sq(particle);
	const double log_d = 0.5 * log(dsq);
	const double LCE = (log_d - log_d0) / time;
	return LCE;
}

double get_MEGNO(Particle * particle){
	return (particle->megno).megno;
}
int main(void){

	const double e0 = 0.1;
	const double a0 = 1;
	const double l0 = M_PI / 3;
	const double pmg0 = 2;
	
	
	Planet planet;
	const double a1 = pow(1.5 , 2./3.);
	planet.mu = 1e-5;
	planet.state.x  = a1;
	planet.state.y  = 0 ;
	planet.state.vx = 0 ;
	planet.state.vy = 1 / sqrt(a1);

	
	Particle particle;
	initialize_phase_state(&(particle.state),a0,l0,e0,pmg0);
	initialize_megno_vars(&(particle.megno));
	const double log_d0 = 0.5 * log( get_delta_sq(&particle) );
	double LCE,megno;
	Simulation sim;
	sim.N_planets = 1;
	sim.N_particles = 1;
	sim.planets = &planet;
	sim.particles= &particle;
	sim.t = 0.;
	const int Nout = 500;
	const double tfin = 2 * M_PI * 1e3;
	const double dt = 2 * M_PI / 30.0;
	int Nsteps = (int) tfin / dt;
	int Nsubstep = Nsteps / Nout;
	for(int i=0; i<Nout ;i++){
		printf("%.8f\t",sim.t);
		print_particle(sim.particles);
		LCE = get_LCE(sim.particles,log_d0,sim.t + dt);
		megno = get_MEGNO(sim.particles);
		printf("%.8g\t%.8g\n",LCE,megno);
		for(int j=0; j<Nsubstep;j++){simulation_step(&sim,dt);}
	}
	return 0;
}
