#MATEUS TEIXEIRA MAGALHÃES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import numpy as np
import enum
import tools
import upwind_schemes as SCHEMES

#class of equations of 'advection_difusion_equation_solver'
class Equation_types(enum.Enum):
    
    Burges = 'burges'
    Boundary_layer = 'boundary_layer'
    Linear_advection = 'linear_advection'

#aproximation of convection velocity in face 'f' for Burges equation
def U_F(U, i):

    return ( U[i+1] + U[i] ) / 2

#aproximation of convection velocity in face 'g' for Burges equation
def U_G(U, i):
    
    return (U[i] + U[i-1]) / 2

#approximation of value 'u' in face 'f' of node 'i'
def u_f(u, i, SCHEME, param, VEL_F, dx):
    
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
    
    #using FOU when normalized variable will give us division by zero
    if (D - R) == 0.0 or (D - R) <= 1e-10:
        
        return SCHEMES.FOU(u, dx, i, VEL_F)
    
    return SCHEME(U, D, R, param)

#approximation of value 'u' in face 'g' of node 'i'
def u_g(u, i, SCHEME, param, VEL_G, dx):
    
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
        
    #using FOU when normalized variable will give us division by zero
    if (D - R) == 0.0 or (D - R) <= 1e-10:
        
        return SCHEMES.FOU(u, dx, i, VEL_G)
        
    return SCHEME(U, D, R, param)
        
def advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func,\
                                       uL, uR, SCHEME, param,\
                                       analitic_sol_func,\
                                       SCHEME_LABEL, marker = '.',\
                                       PATH = 'results/',\
                                       equation_type = Equation_types.Burges,\
                                       a = 1,\
                                       save_step_by_step = False, clean_plot = True,\
                                       ymin = None, ymax = None):
    #spacing discretization size
    dx = (domx[-1] - domx[0]) / (nx-1)  
    
    #time discretization size
    dt = cfl * dx
    
    #number of points in time domain
    nt = int((domt[-1] - domt[0]) / dt) + 1
    
    #matrix to store solution 'u' at time 'i' and 'i + 1'
    M_u = [[0 for i in range(nx)], [0 for i in range(nx)]]
    
    #index of layer in matrix 'M_u' that stores 'u' at time 'i' and 'i + 1'
    p, q = 0, 1
    
    #x axes
    x_axes = []
    
    #initial condition:    
    for i in range(nx):
        
        x_axes.append(domx[0] + i*dx)
                        
        M_u[p][i] = initial_cond_func(x_axes[i])
        M_u[q][i] = initial_cond_func(x_axes[i])
    
    #boundary condition:
    M_u[p][0] = uL
    M_u[q][0] = uL
    M_u[p][-1] = uR
    M_u[q][-1] = uR  
        
    #save initial result   
    if save_step_by_step:
        save_solution(x_axes, M_u[p], analitic_sol_func,\
                      cfl, v, PATH, domt[0], SCHEME_LABEL,\
                      marker = marker, clean_plot=True, ymin = ymin, ymax = ymax)
    
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
            
            if equation_type == Equation_types.Linear_advection:
                
                VEL_F = a
                VEL_G = a
            
            uf = u_f(M_u[q], i, SCHEME, param, VEL_F, dx)
            ug = u_g(M_u[q], i, SCHEME, param, VEL_G, dx)
            
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
                 
            if equation_type == Equation_types.Linear_advection:
                
                M_u[p][i] = - ( dt * a * (uf -  ug) / dx ) + M_u[q][i]
            
        #save result in each iteration
        if save_step_by_step:
            time = t*dt
            save_solution(x_axes, M_u[p], analitic_sol_func,\
                  cfl, v, PATH, time, SCHEME_LABEL,\
                  marker = marker, clean_plot=clean_plot, ymin = ymin, ymax = ymax)
 
    #save result in last iteration
    if (not save_step_by_step):
        save_solution(x_axes, M_u[p], analitic_sol_func,\
                  cfl, v, PATH, domt[-1], SCHEME_LABEL,\
                  marker = marker, clean_plot=clean_plot, ymin = ymin, ymax = ymax)
        


def save_solution(x_axes, u, analitic_sol_func, cfl, v, PATH, time, SCHEME_LABEL,\
                  marker = '.', clean_plot = False, ymin = None, ymax = None):
    
    time_precision = 4
    
    time_string = "{:." + str(time_precision) + "f}"
    
    fileName = PATH + 'result'\
    + '_time=' + time_string.format(time)\
    + '_v=' + str(v)\
    + '_cfl=' + str(cfl) + '.png'
    
    title = 'tempo = ' + str(time)\
            + '  v =' + str(v)\
            + '  cfl = ' + str(cfl)
    
    analitic = None
    
    if analitic_sol_func!=None:
    
        analitic = list(analitic_sol_func(np.array(x_axes), time))        
    
    tools.save_fig(x_axes, u, fileName, title, SCHEME_LABEL,\
                marker = marker,\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = analitic==None and clean_plot, ymin = ymin, ymax = ymax)
    
    if analitic != None:
            
        tools.save_fig(x_axes, analitic, fileName, title, 'analitica',\
                marker = None,\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = True, ymin = ymin, ymax = ymax)

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