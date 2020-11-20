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

#FIRST ORDER upwind scheme, FOU
def FOU(u, dx, i, VEL):
    
    if VEL > 0.0:
            
        return (u[i] - u[i-1]) / dx
    
    else:
        
        return (u[i+1] - u[i]) / dx