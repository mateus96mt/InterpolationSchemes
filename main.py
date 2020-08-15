import matplotlib.pyplot as plt
import math

#normalized variable
def norm_var(PHI, PHI_D, PHI_R):
    
    return (PHI - PHI_R) / (PHI_D - PHI_R)

#FSFL upwind scheme
def FSFL(PHI_U, PHI_D, PHI_R, beta):
    
    PHI_U_norm = norm_var(PHI_U, PHI_D, PHI_R)
    
    if(PHI_U_norm >=0 and PHI_U_norm<=1):
    
        PHI_F_norm = ( ( -2 * beta + 4 ) * ( PHI_U_norm ** 4 ) )\
        + ( ( 4 * beta - 8 ) * ( PHI_U_norm ** 3 ) )\
        + ( ( (-5 * beta + 8) / 2 ) * ( PHI_U_norm ** 2 ) )\
        + ( ( (beta + 2) / 2 ) * PHI_U_norm )
        
        PHI_F = PHI_R + (PHI_D - PHI_R) * PHI_F_norm
        
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
    
    print("[a,b] = [", a, b, "]")

    if(PHI_U_norm >= 0 and PHI_U_norm < a):
        
        return (2 - cfl) * PHI_U_norm
    
    if(PHI_U_norm >= a and PHI_U_norm <= b):
        
        return PHI_U_norm\
        + 0.5 * (1 - abs(cfl)) * (1 - PHI_U_norm)\
        - (1/6) * (1 - (cfl ** 2)) * (1 - (2 * PHI_U_norm))
        
    if(PHI_U_norm > b and PHI_U_norm <= 1):
        
        return 1 - cfl + (cfl * PHI_U_norm)
        
    if(PHI_U_norm < 0 or PHI_U_norm > 1):
        
        return PHI_U_norm

#convection in face 'f'    
def U_F(U, i):

    return ( U[i+1] + U[i] ) / 2

#convection in face 'g'
def U_G(U, i):
    
    return (U[i] + U[i-1]) / 2

#approximation of value 'u' in face 'f' of node 'i'
def u_f(u, i, SCHEME, param):
    
    U, D, R = 0, 0, 0
    
    if( U_F(u, i) > 0 ):
        
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
def u_g(u, i, SCHEME, param):
    
    U, D, R = 0, 0, 0
    
    if( U_G(u, i) > 0 ):
        
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
    
    u[:][0] = uL
    u[:][-1] = uR
    
def apply_initial_condition(u, domx):
    
    n = len(u[0])
    
    d = (domx[-1] - domx[0]) / (n - 1)
    
    for i in range(n):
        
        x = domx[0] + i * d
        
        u[0][i] = math.sin(2 * math.pi * x)
        
        u[1][i] = math.sin(2 * math.pi * x)

def burges_equation_solver(nx, domx, domt, cfl, v, SCHEME, param):
    
    dx = (domx[-1] - domx[0]) / (nx-1)
    
    dt = cfl * dx
    
    u = [[0 for i in range(nx)], [1 for i in range(nx)]]
    
    apply_boundary_condition(u, 0, 0)
    
    apply_initial_condition(u, domx)
    
    p, q = 0, 1
    
    t = domt[0]
    
    params_log(nx, dx, dt, domx, domt, cfl, v, param)
        
    while t < domt[-1]:
        
        p, q = q, p
        
        for i in range(1, len(u[p])-1):
                
            UF = U_F(u[q], i)
            UG = U_G(u[q], i)
            
            uf = u_f(u[q], i, SCHEME, param)
            ug = u_g(u[q], i, SCHEME, param)
            
            u[p][i] = ( dt * v * ( u[q][i-1] + 2 * u[q][i] + u[q][i+1] ) / (dx ** 2) )\
             - ( dt * ( (UF * uf) - (UG * ug) ) / dx )\
             + u[q][i]
        
        t = t + dt
        
    plt.plot([domx[0] + i * dx for i in range(nx)], u[p], '.')
    
#generate log of params values    
def params_log(nx, dx, dt, domx, domt, cfl, v, param):
    
    print("------------ PARAMS LOG ------------")
    print("nx:", nx)
    print("dx:", dx)
    print("dt:", dt)
    print("domx:", domx)
    print("domx:", domt)
    print("cfl(param):", cfl)
    print("cfl(real):", dt / dx)
    print("v:", v)
    print("param:", param)
    print("------------------------------------")
        
def main():
    
    domx = [0, 1]
    
    domt = [0, 0.5]
    
    nx = 200
    
    v = 0.1
    
    cfl = 0.1
    
    beta = 1;
    
    burges_equation_solver(nx, domx, domt, cfl, v, FSFL, beta)
    
main()
        
