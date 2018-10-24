from ctypes import CDLL,sizeof,addressof,c_bool,c_double,c_int,Structure,POINTER,pointer,byref,resize,sizeof,cast
import numpy as np
import os
import matplotlib.pyplot as plt

def ctypes_array_resize(array, new_size):
    resize(array, sizeof(array._type_)*new_size)
    return (array._type_*new_size).from_address(addressof(array))

clib = CDLL("/Users/shadden/Projects/SemiNbody2/libSemiNbody.so")

def get_ctype_ptr(dtype,dim,**kwargs):
    return np.ctypeslib.ndpointer(dtype=dtype,ndim=dim,flags='CONTIGUOUS',**kwargs)
p1d=get_ctype_ptr(np.float,1)
p1dInt = get_ctype_ptr(c_int,1)


PhaseState_field_names="lmbda,Y,Lambda,X,dlmbda,dY,dLambda,dX,dlambda_dot,dY_dot,dLambda_dot,dX_dot".split(",")

class PhaseState(Structure):
    _fields_=[(x,c_double) for x in PhaseState_field_names]

CartesianPhaseStateSimple_field_names="x,y,vx,vy".split(",")      

class CartesianPhaseStateSimple(Structure):
    _fields_=[(x,c_double) for x in CartesianPhaseStateSimple_field_names]

class MEGNO_Auxilary_Variables(Structure):
    _fields_ = [("W",c_double),("Y",c_double),("megno",c_double)]

class Planet(Structure):
    _fields_ = [("mu",c_double),("state",CartesianPhaseStateSimple)]

class Particle(Structure):
    _fields_ = [("state",PhaseState),("megno",MEGNO_Auxilary_Variables)]
class Resonance(Structure):
    _fields_ = [("inner_perturber_Q",c_bool),("j",c_int),("k",c_int),("l",c_int),("fCoeff",c_double)]
class ResonancePerturbation(Structure):
    _fields_ = [("e",c_double),("lmbda",c_double),("pomega",c_double),("mean_motion",c_double),("mu",c_double),\
               ("N_resonances",c_int),("resonances",POINTER(Resonance))]

initialize_particle=clib.initialize_particle
initialize_particle.argtypes = [POINTER(Particle),c_double,c_double,c_double,c_double]
initialize_particle.restype = None
from scipy.optimize import brentq
def alew2xv(a,l,e,w):
    E = 0.5 / a
    hsq = a * (1 - e*e)
    M = np.mod(l - w,2 * np.pi)
    f = lambda u: M - u + e * np.sin(u)
    u = brentq(f,0,2*np.pi)
    _x = a * (np.cos(u) - e)
    _y = a * np.sqrt(1 - e*e) * np.sin(u)
    du_dt  = n / (1 - e * np.cos(u))
    _vx = 
    _vy = 

class Simulation(Structure):
    def __init__(self):
        clib.initialize_simulation(byref(self))
        
    def add_particle(self,a=1.,l=0.,e=0.,w=0.):
        self.N_particles+=1
        N = self.N_particles
        # Better way to do this??
        new_array = (N * Particle)()
        for i in range(N-1):
            new_array[i] = self.particles[i]
        initialize_particle(pointer(new_array[N-1]),a,l,e,w)
        self.particles = pointer(new_array[0])

    def add_planet(self,mu=0,a=1.,l=0.,e=0.,w=0.):
        self.N_planets+=1
        N = self.N_planets

        # Better way to do this??
        new_array = (N * Particle)()
        for i in range(N-1):
            new_array[i] = self.planets[i]
        xv = alew2xv(a,l,e,w)
        self.planets = pointer(new_array[0])

#
Simulation._fields_ = [
        ("N_planets",c_int),
        ("N_resonance",c_int),
        ("N_particles",c_int),
        ("planets",POINTER(Planet)),
        ("resonances",POINTER(ResonancePerturbation)),
        ("particles",POINTER(Particle)),
        ("t",c_double),
        ("dt",c_double)
        ]


# initialize_phase_state(PhaseState * Z,double a, double l, double e, double pomega)
initialize_phase_state = clib.initialize_phase_state
initialize_phase_state.argtypes=[POINTER(PhaseState),c_double,c_double,c_double,c_double]
initialize_phase_state.restype=None

initialize_megno_vars = clib.initialize_megno_vars
initialize_megno_vars.argtypes= [POINTER(MEGNO_Auxilary_Variables)]
initialize_megno_vars.restype=None

integrate_simulation = clib.integrate_simulation
integrate_simulation.argtypes=[POINTER(Simulation),c_double]
integrate_simulation.restype=None

initialize_resonance_multiplet=clib.initialize_resonance_multiplet
initialize_resonance_multiplet.argtypes = [POINTER(Resonance),c_bool,c_int,c_int,c_double]
initialize_resonance_multiplet.restype=None

def init_particle(particle,a,l,e,pmg):
    initialize_phase_state(pointer(particle.state),a,l,e,pmg)
    initialize_megno_vars(pointer(particle.megno))
def init_planet(planet,mu,x,y,vx,vy):
    planet.mu = mu
    planet.state.x = x
    planet.state.y = y    
    planet.state.vx = vx
    planet.state.vy = vy
def init_resonance(res,jlist,klist,innerQ,mu,mean_motion,ecc,lmbda,pomega):
    res.mu = mu
    res.mean_motion = mean_motion
    res.e = ecc
    res.lmbda = lmbda
    res.pomega = pomega


    n_subresonance = np.sum(np.array(klist)+1)
    res.N_resonances = n_subresonance

    res_array = (n_subresonance * Resonance)()
    indx=0
    alpha = (mean_motion)**(-2/3)
    if alpha>1:
        alpha=1/alpha
    for j,k in zip(jlist,klist):
        initialize_resonance_multiplet(pointer(res_array[indx]),innerQ,j,k,alpha)
        indx += 1 + k
    res.resonances = pointer(res_array[0])
