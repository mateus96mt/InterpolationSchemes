#MATEUS TEIXEIRA MAGALHÃES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import enum
import tools
import os

#class of equations of 'advection_difusion_equation_solver'
class Equation_types(enum.Enum):
    
    Burges = 'burges'
    Boundary_layer = 'boundary_layer'
    Linear_advection = 'linear_advection'

#normalized variable
def norm_var(PHI, PHI_D, PHI_R):
    
    return (PHI - PHI_R) / (PHI_D - PHI_R)

#aproximation of convection velocity in face 'f' for Burges equation
def U_F(U, i):

    return ( U[i+1] + U[i] ) / 2

#aproximation of convection velocity in face 'g' for Burges equation
def U_G(U, i):
    
    return (U[i] + U[i-1]) / 2

#approximation of value 'u' in face 'f' of node 'i'
def u_f(u, i, SCHEME, param, VEL_F, dx):
    
    phi_u, phi_d, phi_r = 0, 0, 0
    
    if( VEL_F > 0 ):
        
        phi_d = u[i+1]
        phi_r = u[i-1]
        phi_u = u[i]
    
    else:
        
        phi_d = u[i]
        
        if i + 2 <= len(u) - 1:
            
            phi_r = u[i+2]
            
        else:
            
            phi_r = u[i+1]
            
        phi_u = u[i+1]
    
    #using FOU when normalized variable will give us division by zero
    if (phi_d - phi_r) == 0.0 or (phi_d - phi_r) <= 1e-10:
        
        return phi_u
    
    else:
        
        PHI_U_norm = norm_var(phi_u, phi_d, phi_r)
        
        phi_f_norm = SCHEME(PHI_U_norm, param)
        
        phi_f = phi_r + ((phi_d - phi_r) * phi_f_norm)
        
        return phi_f

#approximation of value 'u' in face 'g' of node 'i'
def u_g(u, i, SCHEME, param, VEL_G, dx):
    
    phi_u, phi_d, phi_r = 0, 0, 0
    
    if( VEL_G > 0 ):
        
        phi_d = u[i]
        
        if i - 2 >=0:
        
            phi_r = u[i-2]
        
        else:
            
            phi_r = u[i-1]
        
        phi_u = u[i-1]
    
    else:
        
        phi_d = u[i-1]
        phi_r = u[i+1]
        phi_u = u[i]
        
    #using FOU when normalized variable will give us division by zero
    if (phi_d - phi_r) == 0.0 or (phi_d - phi_r) <= 1e-10:
        
        return phi_u
    
    else:
        
        PHI_U_norm = norm_var(phi_u, phi_d, phi_r)
        
        phi_f_norm = SCHEME(PHI_U_norm, param)
        
        phi_f = phi_r + ((phi_d - phi_r) * phi_f_norm)
        
        return phi_f
        
def advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func,\
                                       uL, uR, SCHEME, param,\
                                       analitic_sol_func,\
                                       folderName,\
                                       equation_type = Equation_types.Burges,\
                                       a = 1,\
                                       step_interval = 0.05,\
                                       errorType = 2):
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
    
    params_log(nx, dx, nt, dt, domx, domt, cfl, v, param)
    
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
    save_solution(folderName, x_axes, M_u[p], analitic_sol_func, cfl, v, 0,
                  nx, nt, dx, dt, domt)
    
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
        save_solution(folderName, x_axes, M_u[p], analitic_sol_func, cfl, v, t,
                      nx, nt, dx, dt, domt)
 
    #save result in last iteration
    save_solution(folderName, x_axes, M_u[p], analitic_sol_func, cfl, v, nt,
                  nx, nt, dx, dt, domt)
        
        
    return tools.calculateError(analitic_sol_func,
                                M_u[p], nx, domx, dx, domt[-1],
                                tipo=errorType)


def save_solution(folderName, x_axes, u, analitic_sol_func, cfl, v, t, 
                  nx, nt, dx, dt, domt):
    
    time = domt[0] + t*dt
    
    #create folder if not exist
    if not os.path.exists(folderName):
        os.makedirs(folderName)
    
    if t == nt:
        
        output = open(folderName + "/" + "FINAL" + ".data","w")
        print(folderName + "/" + "FINAL" + ".data")
        
    else:
        
        output = open(folderName + "/" + str(t) + ".data","w")
        print(folderName + "/" + str(t) + ".data")

    
    output.write("time = " + str(time)+
                 "    cfl = " + str(cfl) + 
                 "    nx = " + str(nx) +
                 "    nt = " + str(nt) +
                 "    dx = " + str(dx) +
                 "    dt = " + str(dt) +"\n")
    
    for i in range(len(x_axes)):
        
        output.write(str(x_axes[i]) +
                     " " + str(analitic_sol_func(x_axes[i], time)) + 
                     " " + str(u[i]) +
                     "\n")

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