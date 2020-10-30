#MATEUS TEIXEIRA MAGALHÃES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import matplotlib.pyplot as plt
import math
import os
import numpy as np
import enum

#class of equations of 'advection_difusion_equation_solver'
class Equation_types(enum.Enum):
    
    Burges = 'burges'
    Boundary_layer = 'boundary_layer'
   
#-----------------------------------ANALITIC SOLUTIONS------------------------

def analitic_boundary_layer_equation(x, v):
    
    return (1 - np.exp(x / v)) / (1 - np.exp(1 / v))

#-----------------------------------------------------------------------------

#normalized variable
def norm_var(PHI, PHI_D, PHI_R):
    
    return (PHI - PHI_R) / (PHI_D - PHI_R)

#TOPUS upwind scheme
def TOPUS(PHI_U, PHI_D, PHI_R, alpha):
    
    PHI_U_norm = norm_var(PHI_U, PHI_D, PHI_R)
    
    if(PHI_U_norm >=0 and PHI_U_norm<=1):
    
        PHI_F_norm = ( alpha * (PHI_U_norm ** 4) )\
        + ( ((-2 * alpha) + 1) * (PHI_U_norm ** 3) )\
        + ( (( (5 * alpha) - 10) / 4 ) * (PHI_U_norm ** 2) )\
        + ( ( (-alpha + 10) / 4) * PHI_U_norm)
        
        PHI_F = PHI_R + ((PHI_D - PHI_R) * PHI_F_norm)
        
        return PHI_F
    
    else:
        
        return PHI_U

#FSFL upwind scheme
def FSFL(PHI_U, PHI_D, PHI_R, beta):
    
    PHI_U_norm = norm_var(PHI_U, PHI_D, PHI_R)
    
    if(PHI_U_norm >=0 and PHI_U_norm<=1):
    
        PHI_F_norm = ( ( (-2 * beta) + 4 ) * ( PHI_U_norm ** 4 ) )\
        + ( ( (4 * beta) - 8 ) * ( PHI_U_norm ** 3 ) )\
        + ( ( ((-5 * beta) + 8) / 2 ) * ( PHI_U_norm ** 2 ) )\
        + ( ( (beta + 2) / 2 ) * PHI_U_norm )
        
        PHI_F = PHI_R + ((PHI_D - PHI_R) * PHI_F_norm)
        
        return PHI_F
    
    else:
        
        return PHI_U
    
#ADBQUICKEST upwind scheme
def ADBQUICKEST(PHI_U, PHI_D, PHI_R, cfl):
    
    PHI_U_norm = norm_var(PHI_U, PHI_D, PHI_R)
    
    a = ( 2 - (3 * abs(cfl)) + (cfl ** 2) )\
    / ( 7 - (6 * cfl) - (3 * abs(cfl)) + (2 * (cfl **2)))
    
    b = ( -4 + (6 *cfl) - (3 * abs(cfl)) + (cfl ** 2))\
    / ( -5 + (6 * cfl) - (3 * abs(cfl)) + (2 * (cfl **2)))
    
    PHI_F_norm = 0
    
    if(PHI_U_norm >= 0 and PHI_U_norm < a):
        
        PHI_F_norm = (2 - cfl) * PHI_U_norm
    
    if(PHI_U_norm >= a and PHI_U_norm <= b):
        
        PHI_F_norm = PHI_U_norm\
        + 0.5 * (1 - abs(cfl)) * (1 - PHI_U_norm)\
        - (1/6) * (1 - (cfl ** 2)) * (1 - (2 * PHI_U_norm))
        
    if(PHI_U_norm > b and PHI_U_norm <= 1):
        
        PHI_F_norm = 1 - cfl + (cfl * PHI_U_norm)
        
    if(PHI_U_norm < 0 or PHI_U_norm > 1):
        
        PHI_F_norm = PHI_U_norm
        
    PHI_F = PHI_R + ((PHI_D - PHI_R) * PHI_F_norm)
        
    return PHI_F

#convection in face 'f'    
def U_F(U, i):

    return ( U[i+1] + U[i] ) / 2

#convection in face 'g'
def U_G(U, i):
    
    return (U[i] + U[i-1]) / 2

#approximation of value 'u' in face 'f' of node 'i'
def u_f(u, i, SCHEME, param, VEL_F):
    
    U, D, R = 0, 0, 0
    
    if( VEL_F > 0 ):
        
        D = u[i+1]
        R = u[i-1]
        U = u[i]
    
    else:
        
        D = u[i]
        
        if i + 2 <= len(u) - 1:
            
            R = u[i+2]
            
        else:
            
            R = u[i+1]
            
        U = u[i+1]
        
    return SCHEME(U, D, R, param)

#approximation of value 'u' in face 'g' of node 'i'
def u_g(u, i, SCHEME, param, VEL_G):
    
    U, D, R = 0, 0, 0
    
    if( VEL_G > 0 ):
        
        D = u[i]
        
        if i - 2 >=0:
        
            R = u[i-2]
        
        else:
            
            R = u[i-1]
        
        U = u[i-1]
    
    else:
        
        D = u[i-1]
        R = u[i+1]
        U = u[i]
        
    return SCHEME(U, D, R, param)

def apply_boundary_condition(u, uL, uR):
    
    u[0] = uL
    u[-1] = uR
        
def apply_initial_condition(u, domx, equation_type =Equation_types.Burges):
    
    n = len(u)
    
    d = (domx[-1] - domx[0]) / (n - 1)
    
    for i in range(n):
        
        x = domx[0] + i * d
        
        if equation_type == Equation_types.Burges:
            
            u[i] = math.sin(2 * math.pi * x)
            
        if equation_type == Equation_types.Boundary_layer:
            
            u[i] = 0
        

def advection_difusion_equation_solver(nx, dx, domx, domt, cfl, v, SCHEME, param, \
                                       SCHEME_LABEL, marker = '.', PATH = 'results/',
                                       equation_type = Equation_types.Burges, a = 1):
        
    #time discretization size
    dt = cfl * dx
    
    #number of points in time domain
    nt = int((domt[-1] - domt[0]) / dt) + 1
    
    #matrix to store solution 'u' at time 'i' and 'i + 1'
    M_u = [[0 for i in range(nx)], [0 for i in range(nx)]]
    
    #index of layer in matrix 'M_u' that stores 'u' at time 'i' and 'i + 1'
    p, q = 0, 1
    
    apply_initial_condition(M_u[p], domx, equation_type)
    apply_initial_condition(M_u[q], domx, equation_type)
    
    uL, uR = 0, 0
    
    if equation_type == Equation_types.Burges:
        
        uL = 0
        uR = 0
        
    if equation_type == Equation_types.Boundary_layer:
        
        uL = 0
        uR = 1
        
    apply_boundary_condition(M_u[p], uL, uR)
    apply_boundary_condition(M_u[q], uL, uR)
        
    params_log(nx, dx, nt, dt, domx, domt, cfl, v, param)
    
    #loop in time
    for t in range(1, nt):
        
        #switching index of layer in matrix 'M_u' that stores 'u'
        p, q = q, p
        
        #loop in space
        for i in range(1, nx - 1):
                
            #advection velocities in face f and g, respectively
            VEL_F, VEL_G = 0, 0
            
            if equation_type == Equation_types.Burges:
            
                VEL_F = UF = U_F(M_u[q], i)
                VEL_G = UG = U_G(M_u[q], i)
            
            if equation_type == Equation_types.Boundary_layer:
                
                VEL_F = a
                VEL_G = a
            
            uf = u_f(M_u[q], i, SCHEME, param, VEL_F)
            ug = u_g(M_u[q], i, SCHEME, param, VEL_G)
            
            if equation_type == Equation_types.Burges:
            
                M_u[p][i] = ( dt * v * ( M_u[q][i-1] - 2 * M_u[q][i] + M_u[q][i+1] )\
                   / (dx ** 2) )\
                 - ( dt * 0.5 * ( (UF * uf) - (UG * ug) ) / dx )\
                 + M_u[q][i]
                 
            if equation_type == Equation_types.Boundary_layer:
                
                M_u[p][i] = ( dt * v * ( M_u[q][i-1] - 2 * M_u[q][i] + M_u[q][i+1] )\
                   / (dx ** 2) )\
                 - ( dt * a * (uf -  ug) / dx )\
                 + M_u[q][i]
            
        
        #save result in last iteration
        if t == nt-1:
            
            fileName = PATH + 'result'\
            + '_time=' + str(domt[-1])\
            + '_v=' + str(v)\
            + '_cfl=' + str(cfl) + '.png'
            save_result(domx, nx, dx, M_u[p], fileName, domt[-1], v, cfl,\
                        SCHEME_LABEL, marker)
        
#plot result and save in figure
def save_result(domx, nx, dx, u, fileName, time, v, cfl, label, 
                marker = '.', ylimD = -2, ylimU = 2):
    
    #create folder if not exist
    if not os.path.exists(fileName.split('/')[0]):
        os.makedirs(fileName.split('/')[0])
    
    title = 'tempo = ' + str(time)\
            + '  v =' + str(v)\
            + '  cfl = ' + str(cfl)
    plt.title(title)
    
    plt.ylim([ylimD, ylimU])
    
    plt.xlabel('x')
    plt.ylabel('u')
    
    plt.plot([domx[0] + i * dx for i in range(nx)], u, marker = marker,\
              label = label)
    plt.legend(loc = "best")
    plt.grid(True)
    plt.tight_layout()
    
    print('saving ', fileName)
    plt.savefig(str(fileName), dpi = 500)
    
#generate log of params values    
def params_log(nx, dx, nt, dt, domx, domt, cfl, v, param):
    
    print("\n\n------------ PARAMS LOG ------------")
    print("nx:", nx)
    print("dx:", dx)
    print("nt:", nt)
    print("dt:", dt)
    print("domx:", domx)
    print("domt:", domt)
    print("cfl(param):", cfl)
    print("cfl(real):", dt / dx)
    print("v:", v)
    print("param:", param)
    print("------------------------------------")

def TRAB1_CFD_MESTRADO():
    
    #domain in x direction
    domx = [0, 1]
    
    #number of points in x direction
    nx = 200
    
    #spacing discretization size
    dx = (domx[-1] - domx[0]) / (nx-1)
    
    #FSFL param
    beta = 0;
    
    #TOPUS param
    alpha = 2;
    
    #domains in time
    domts = [[0, 0.1], [0, 0.3], [0, 0.5]]
    
    #Courant numbers
    cfls = [0.1, 0.5, 0.9]
    
    #viscosity numbers
    vs = [0.001, 0.0005, 0.00025]
    
    #plot and save initial condition
    initial_condition = [0 for i in range(nx)]
    apply_initial_condition(initial_condition, domx)
    name = 'results_TRAB1_CFD_MESTRADO/initial.png'
    save_result(domx, nx, dx,\
                initial_condition,\
                name, 0, '-', '-', 'Condicao inicial')
    
    #clean plot
    plt.cla()
    plt.clf()
    
    #-----------------------EXERCISE 2-----------------------
    
    #path to save results
    PATH = 'results_TRAB1_CFD_MESTRADO/EXE2_'
    
    #viscosity
    v = 0.001
    
    for cfl in cfls:
        
        for domt in domts:
            
            #seting upwind scheme, label and parameter
            SCHEME = FSFL
            SCHEME_LABEL = 'FSFL'
            param = beta
            advection_difusion_equation_solver(nx, dx, domx, domt, cfl, v, \
                                               SCHEME, param, SCHEME_LABEL, marker = '.', \
                                               PATH = PATH)
            
            #seting upwind scheme, label and parameter
            SCHEME = ADBQUICKEST
            SCHEME_LABEL = 'ADBQUICKEST'
            param = cfl
            advection_difusion_equation_solver(nx, dx, domx, domt, cfl, v, \
                                               SCHEME, param, SCHEME_LABEL, marker = None, \
                                               PATH = PATH)
            
            #clean plot
            plt.cla()
            plt.clf()
    
    #--------------------------------------------------------
    
    
    #-----------------------EXERCISE 3-----------------------
    
    #path to save results
    PATH = 'results_TRAB1_CFD_MESTRADO/EXE3_'
    
    #Courant number
    cfl = 0.5
    
    for v in vs:
        
        for domt in domts:
            
            #seting upwind scheme, label and parameter
            SCHEME = FSFL
            SCHEME_LABEL = 'FSFL'
            param = beta
            advection_difusion_equation_solver(nx, dx, domx, domt, cfl, v, \
                                               SCHEME, param, SCHEME_LABEL, marker = '.', \
                                               PATH = PATH)
            
            #seting upwind scheme, label and parameter
            SCHEME = ADBQUICKEST
            SCHEME_LABEL = 'ADBQUICKEST'
            param = cfl
            advection_difusion_equation_solver(nx, dx, domx, domt, cfl, v, \
                                               SCHEME, param, SCHEME_LABEL, marker = None, \
                                               PATH = PATH)
            
            #clean plot
            plt.cla()
            plt.clf()
    
    #--------------------------------------------------------


def boundary_layer_equation_test():
    
    cfl = 0.1
    
    nx = 100
    
    v = 0.25
    
    domx = [0, 1]
    
    domt = [0, 1]
    
    dx = (domx[-1] - domx[0]) / (nx-1)
 
    PATH = 'test/'
    
    n_v = 1
    
    beta = 0
    
    for i in range(n_v):
        
        analitic = [analitic_boundary_layer_equation(i*dx, v) for i in range(nx)]
        
        fileName = PATH + 'boundary_layer'\
            + '_time=' + str(0)\
            + '_v=' + str(v)\
            + '_cfl=' + str(cfl) + '.png'
        
        save_result(domx, nx, dx, analitic, fileName, 0, v, cfl, 
                'analitica', marker = None, ylimD=-0.5, ylimU=1.5)
        
        SCHEME = FSFL
        SCHEME_LABEL = 'FSFL'
        param = beta
        
        advection_difusion_equation_solver(nx, dx, domx, domt, cfl, v, \
                   SCHEME, param, SCHEME_LABEL, marker = '.', \
                   PATH = PATH, equation_type = Equation_types.Boundary_layer)
        
        
        #clean plot
        plt.cla()
        plt.clf()
        
        v = v / 10
        
#calling main function
TRAB1_CFD_MESTRADO()
#boundary_layer_equation_test()
