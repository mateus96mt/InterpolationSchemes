#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

#TOPUS upwind scheme
def TOPUS(PHI_U_norm, alpha):
    
    if(PHI_U_norm >=0.0 and PHI_U_norm<=1.0):
    
        return ( alpha * (PHI_U_norm ** 4) )\
        + ( ((-2 * alpha) + 1) * (PHI_U_norm ** 3) )\
        + ( (( (5 * alpha) - 10) / 4 ) * (PHI_U_norm ** 2) )\
        + ( ( (-alpha + 10) / 4) * PHI_U_norm)

    else:
        
        return PHI_U_norm

#FSFL upwind scheme
def FSFL(PHI_U_norm, beta):
        
    if(PHI_U_norm >=0.0 and PHI_U_norm<=1.0):
    
        return ( ( (-2 * beta) + 4 ) * ( PHI_U_norm ** 4 ) )\
        + ( ( (4 * beta) - 8 ) * ( PHI_U_norm ** 3 ) )\
        + ( ( ((-5 * beta) + 8) / 2 ) * ( PHI_U_norm ** 2 ) )\
        + ( ( (beta + 2) / 2 ) * PHI_U_norm )
    
    else:
        
        return PHI_U_norm
    
#ADBQUICKEST upwind scheme
def ADBQUICKEST(PHI_U_norm, cfl):
        
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
        
    return PHI_F_norm

#FIRST ORDER upwind scheme, FOU
def FOU(u, dx, i, VEL):
    
    if VEL > 0.0:
            
        return (u[i] - u[i-1]) / dx
    
    else:
        
        return (u[i+1] - u[i]) / dx
    
    
    
    
    
    












