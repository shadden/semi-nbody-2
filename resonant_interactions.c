#include "SemiNBody.h"
#define NDIM 4
#define INDX(ROW,COL) NDIM * ROW + COL
double mpow(double,int);
double mpow2(double b, int exp){
   if (exp > 2)
   {
      if (exp%2)
         return b*mpow(b*b,(int)exp/2);
      else
         return mpow(b*b,exp/2);
   }
   else if (2 == exp)
      return b*b;
   else if (1 == exp)
      return b;

   return 1.0; // exp == 0
} 
double mpow(double b, int exp){
   if (exp < 0)
   {
      b = 1/b;
      exp *= -1;
   }
   
   return mpow2(b,exp);
}

void get_re_im_power(const double X, const double Y,const int p, double * re_p,double * im_p){
	*re_p=0;
	*im_p=0;
	int sgn = 1;
	int p_choose_k;
	for(int k=0;k<=p;k++){
		p_choose_k = binomialCoeff(p,k);
		if(k%2){
			*(im_p)+= p_choose_k * mpow(X,p-k) * mpow(Y,k) * sgn;
			sgn*= -1;
		}
		else{
			*(re_p)+= p_choose_k * mpow(X,p-k) * mpow(Y,k) * sgn;	
		}
	}
}
void derivs_from_single_resonance(PhaseState* particle,Resonance * res,double e_pert, double lmbda_pert,double pmg_pert, double * derivs,double * jac){
	const int j = res->j;
	const int k = res->k;
	const int l = res->l;
	double lmbda_in,lmbda_out;
	int p,p_pert,dphi_dlmbda;

	if(res->inner_perturber_Q){
		lmbda_in = lmbda_pert;
		lmbda_out = particle->lmbda;
		p = k-l;
		p_pert = l;
		dphi_dlmbda = j;

	}else{
		lmbda_in = particle->lmbda;
		lmbda_out = lmbda_pert;
		p = l;
		p_pert = k-l;
		dphi_dlmbda = k-j;
	}
	const double phi = j * lmbda_out -(j-k) * lmbda_in - p_pert * pmg_pert;
	const double cos_phi = cos(phi);
	const double sin_phi = sin(phi);
	const double e_pert_p_pert = mpow(e_pert,p_pert);
	const double eps = (-1) * (res->fCoeff) * e_pert_p_pert;
	double Xp,Yp;

	const double X = particle->X;
	const double Y = particle->Y;

	get_re_im_power(X,Y,p,&Xp,&Yp);
	
	const double denom = X*X + Y*Y;
	const double Xp_minus_1 = (X * Xp + Y * Yp) / denom;
	const double Yp_minus_1 = (X * Yp - Y * Xp) / denom;
	const double Xp_minus_2 = (X * Xp_minus_1 + Y * Yp_minus_1) / denom;
	const double Yp_minus_2 = (X * Yp_minus_1 - Y * Xp_minus_1) / denom;

	// H ~ eps * ( Xp * cos_phi - Yp * sin_phi)
	
	*(derivs + 0 ) = 0;
	*(derivs + 1 ) = eps * p *  (Xp_minus_1 * cos_phi - Yp_minus_1 * sin_phi);
	*(derivs + 2 ) = eps * dphi_dlmbda * (Xp * sin_phi  + Yp * cos_phi);
	*(derivs + 3 ) = eps * p * (Xp_minus_1 * sin_phi + Yp_minus_1 * cos_phi);

	fill_zeroes(jac,NDIM*NDIM); 

	*(jac+INDX(1,0)) = eps * p * dphi_dlmbda * (-1) * (sin_phi * Xp_minus_1 + cos_phi * Yp_minus_1);
	*(jac+INDX(1,1)) = eps * p * (p-1) * (-1) * (sin_phi * Xp_minus_2 + cos_phi * Yp_minus_2);
	*(jac+INDX(1,3)) = eps * p * (p-1) * (cos_phi * Xp_minus_2 - sin_phi * Yp_minus_2);

	*(jac+INDX(2,0)) = eps * dphi_dlmbda * dphi_dlmbda * (Xp * cos_phi - Yp * sin_phi)  ; 
	*(jac+INDX(2,1)) = eps * p * dphi_dlmbda * (Xp_minus_1 * cos_phi - Yp_minus_1 * sin_phi)  ; 
	*(jac+INDX(2,3)) = eps * p * dphi_dlmbda * (Xp_minus_1 * sin_phi + Yp_minus_1 * cos_phi)  ; 

	*(jac+INDX(3,0)) = eps * p * dphi_dlmbda * (Xp_minus_1 * cos_phi -  Yp_minus_1 * sin_phi);
	*(jac+INDX(3,1)) = eps * p * (p-1) * (Xp_minus_2 * cos_phi - Yp_minus_2 * sin_phi);
	*(jac+INDX(3,3)) = eps * p * (p-1) * (Xp_minus_2 * sin_phi + Yp_minus_2 * cos_phi);

}
void add_resonance_derivs(PhaseState* particle_state, ResonancePerturbation* restrict res_pert ,double* derivs,double* jacobian){
	const int Nres = res_pert->N_resonances;
	Resonance* resonance;
	const double e_pert = res_pert->e;
	const double lmbda_pert = res_pert->lmbda;
	const double pmg_pert = res_pert->pomega;
	const double mu_pert = res_pert->mu;
	double single_res_derivs[NDIM],single_res_jac[NDIM*NDIM];
	for (int i=0;i<Nres;i++){
		resonance = ((res_pert->resonances)+i);
		derivs_from_single_resonance(particle_state, resonance, e_pert, lmbda_pert, pmg_pert,single_res_derivs , single_res_jac);
		add_derivs_and_jacobians(derivs,jacobian,single_res_derivs,single_res_jac,mu_pert);

	}
}


void TestResDerivs(PhaseState* particle_state, const double eps, const int j, const double phi_pert, double * derivs, double* jac){
	const double X0 = particle_state->X;
	const double Y0 = particle_state->Y;
	const double l0 = particle_state->lmbda;

	const double j_psi = j * (l0 - phi_pert); 
	const double eps_cos_j_psi = eps * cos(j_psi); 
	const double eps_sin_j_psi = eps * sin(j_psi); 
	
	*(derivs + 0) = 0;
	*(derivs + 1) = eps_cos_j_psi;
	*(derivs + 2) = j * X0 * eps_sin_j_psi + j * Y0 * eps_cos_j_psi;
	*(derivs + 3) = eps_sin_j_psi;
	
	fill_zeroes(jac,NDIM*NDIM);
	*(jac+INDX(1,0)) = -1 * j * eps_sin_j_psi;
	*(jac+INDX(2,0)) =  j * j * eps_cos_j_psi * X0 - j * j * eps_sin_j_psi * Y0;
	*(jac+INDX(3,0)) =  j * eps_cos_j_psi;
	*(jac+INDX(2,1)) =  j * eps_cos_j_psi;
	*(jac+INDX(2,3)) =  j * eps_sin_j_psi;
}

#if 0
void cosine_sine_array(const double cosx,const double sinx, const int Nmax, double* cosarray, double* sinarray){
	assert(Nmax>=0);
	*(cosarray) = 1;
	*(sinarray) = 0;
	for(int k=1; k<Nmax; k++){
		*(cosarray+k)  = cosx * (*(cosarray+k-1)) - sinx * (*(sinarray+k-1));
		*(sinarray+k)  = sinx * (*(cosarray+k-1)) + cosx * (*(sinarray+k-1));
	}
}
void ResonanceDerivs(PhaseState* particle_state, ResonanceData* restrict resonance_state,double* derivs,double* jacobian){
	const bool inner_Q = rIn->inner_Q;
	int j,k,l;
	double coeff;

	const double X0 = particle_state->X;
	const double Y0 = particle_state->Y;
	const double l0 = particle_state->l;

	const double Esq = 0.5*(X0*X0 + Y0*Y0);
	const double E = sqrt(Esq);
	const double g0 = atan2(Y0,X0);

	double Xdot=0;
	double Ydot=0;
	double Ldot=0;
	double theta,factor;
	double costheta,sintheta;
	double DLdotDl=0,DYdotDl=0,DXdotDl=0,DXdotDY=0,DYdotDX=0,DXdotDX=0;

	// Add up inner planet effects 	
	// get cosine and sine data for

	double cos_l_array[MAX_J]; 
	double sin_l_array[MAX_J]; 
	cosine_sine_array(cos(l0),sin(l0),(rIn->MaxJ)+1,cos_l_array,sin_l_array);

	double cos_g_array[MAX_ORDER+1]; 
	double sin_g_array[MAX_ORDER+1]; 
	cosine_sine_array(cos(g0),sin(g0),(rIn->MaxOrder)+1, cos_g_array,sin_g_array);
	
	double cos_n1t_array[MAX_ORDER+1]; 
	double sin_n1t_array[MAX_ORDER+1]; 
	cosine_sine_array(cos(n1*t),sin(n1*t),(rIn->MaxOrder)+1, cos_n1t_array,sin_n1t_array);

	double c_j_l;
	double s_j_l;
	double c_p_n1t;
	double s_p_n1t;
	double c_op1_g;
	double s_op1_g;
	double c_g = cos_g_array[1];
	double s_g = sin_g_array[1];
	
	for(int i=0; i<NresIn; i++){

		j = *(rIn->ResonanceIndices + 3*i );
		o  = *(rIn->ResonanceIndices + 3*i + 1 );
		p = *(rIn->ResonanceIndices + 3*i + 2 );

		coeff = *( rIn->ResonanceCoefficients + ( MAX_ORDER + 1 )*i + p );

		
		c_j_l = cos_l_array[j];
		s_j_l = sin_l_array[j];
		c_p_n1t = cos_n1t_array[p];
		s_p_n1t = sin_n1t_array[p];
		c_op1_g = o-p-1 >= 0 ? cos_g_array[o-p-1] : c_g;
		s_op1_g = o-p-1 >= 0 ? sin_g_array[o-p-1] : -1*s_g;
		theta = j * l0 + p * n1 * t + (o - p - 1) * g0;			
		costheta = c_j_l * ( c_p_n1t * c_op1_g - s_p_n1t * s_op1_g ) - s_j_l * ( c_p_n1t * s_op1_g + s_p_n1t * c_op1_g);		
		sintheta = c_j_l * ( c_p_n1t * s_op1_g + s_p_n1t * c_op1_g ) + s_j_l * ( c_p_n1t * c_op1_g - s_p_n1t * s_op1_g);
	
#if PRINT
		printf("inner %d %d %d: %g \t %g \n",j,o,p,costheta-cos(theta),sintheta-sin(theta) );
#endif
	
		// Derivatives
		factor  = o >= p+1 ? -RT2 * mu1 * coeff * (o-p) * mpow(e1,p) * mpow(E,o-p-1) : 0;		
		Ydot += factor * costheta;
		Xdot += factor * sintheta;
		Ldot +=  -2 * mu1  * coeff * j * mpow(e1,p) * mpow(E,o-p) * (sintheta * c_g + costheta * s_g );
		
		// Variationals
		DYdotDl += -factor * j * sintheta;
		DXdotDl +=  factor * j * costheta;
		DLdotDl += -2 * mu1  * coeff * j*j * mpow(e1,p) * mpow(E,o-p) * ( costheta*c_g-sintheta*s_g );		

		factor  = o >= p+2 ? -2 * mu1 * coeff * (o-p) * (o-p-1) * mpow(e1,p) * mpow(E,o-p-2) : 0;
		
		DXdotDX +=  0.5 * factor * (sintheta*c_g-costheta*s_g) ;	
 		DYdotDX +=  0.5 * factor * ((costheta * c_g) + (sintheta * s_g)) ;
 		DXdotDY +=  0.5 * factor * ((costheta * c_g) + (sintheta * s_g)) ;
		

	}
	// Coordinates z_i are:
	//	i	coord
	//	-	-----
	//	0	l
	//	1	Y
	//	2	L
	//	3	X
	
	*(derivs) = 0.;
	*(derivs+1) = Ydot;
	*(derivs+2) = Ldot;
	*(derivs+3) = Xdot;
	

	double Hij[4][4];
	Hij[0][0] =  -DLdotDl ; // H_l,l
	Hij[0][1] = -DXdotDl ; // H_l,Y
	Hij[0][2] =  0.	   ; // H_l,L		
	Hij[0][3] =  DYdotDl ; // H_l,X	

	Hij[1][1] = -DXdotDY ; // H_Y,Y
	Hij[1][2] = 0.       ; // H_Y,L
	Hij[1][3] = -DXdotDX  ; // H_Y,X

	Hij[2][2] = 0.       ; // H_L,L
	Hij[2][3] = 0.       ; // H_L,X

	Hij[3][3] = DYdotDX  ; // H_X,X
	
	int row,col;
	
	double jacobian_ij;
	for(row=0;row<4;row++){
		for(col=0;col<4;col++){
		jacobian_ij = 0;
		for(int l=0; l<4; l++){
			if(l>col){
				Hij[l][col] = Hij[col][l];
			}
			jacobian_ij+= symplecticJ[row][l] * Hij[l][col];
		}
				*( jacobian + INDX(row,col) ) = jacobian_ij;
		}
	}

}
#endif
