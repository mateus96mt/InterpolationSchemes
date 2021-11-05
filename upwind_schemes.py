#MATEUS TEIXEIRA MAGALHAES - UFJF - PGMC
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

#EPUS upwind scheme  
def EPUS(PHI_U_norm, lam):
    
    if(PHI_U_norm >=0.0 and PHI_U_norm<=1.0):

        return (-4 * (lam - 24.0) * (PHI_U_norm ** 8))\
            + (16.0 * (lam - 23.0) * (PHI_U_norm ** 7))\
            + ((528.0 - (25 * lam)) * (PHI_U_norm ** 6))\
            + ( ((19.0 * lam) - 336.0) * (PHI_U_norm ** 5))\
            + ((80.0 - (7.0 * lam)) * (PHI_U_norm ** 4))\
            + (lam * (PHI_U_norm ** 3))\
            + PHI_U_norm
    
    else:
    
        return PHI_U_norm


def SDPUS_C1(PHI_U_norm, gamma):
    
    if(PHI_U_norm >=0.0 and PHI_U_norm<=1.0):

        return ((-24.0 + (4 * gamma)) * (PHI_U_norm ** 6))\
                + ((68.0 - (12.0 * gamma)) * (PHI_U_norm ** 5))\
                + ((-64.0 + (13 * gamma)) * (PHI_U_norm ** 4))\
                + ((20.0 - (6 * gamma)) * (PHI_U_norm ** 3))\
                + (gamma * (PHI_U_norm ** 2))\
                + PHI_U_norm
    
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
    
def CUBISTA(PHI_U_norm, param):
    
    if(PHI_U_norm >= 0.0 and PHI_U_norm<0.375):
        
        return 1.75 * PHI_U_norm
    
    if(PHI_U_norm >= 0.375 and PHI_U_norm <= 0.75):
        
        return (0.75 * PHI_U_norm) + 0.375
    
    if(PHI_U_norm > 0.75 and PHI_U_norm <= 1.0):
        
        return (0.25 * PHI_U_norm) + 0.75
    
    return PHI_U_norm

def VONOS(PHI_U_norm, param):
    
    if(PHI_U_norm >= 0.0 and PHI_U_norm<(3.0/74.0)):
        
        return 10.0 * PHI_U_norm
    
    if(PHI_U_norm >= (3.0/74.0) and PHI_U_norm < 0.5):
        
        return (0.375 * (1.0 + (2.0*PHI_U_norm)))
    
    if(PHI_U_norm >= 0.5 and PHI_U_norm < (2.0/3.0)):
        
        return (1.5 * PHI_U_norm)
    
    if(PHI_U_norm >= (2.0/3.0) and PHI_U_norm <= 1.0):
        
        return 1.0
    
    return PHI_U_norm
    
def WACEB(PHI_U_norm, param):
    
    if(PHI_U_norm >= 0.0 and PHI_U_norm<0.3):
        
        return 2.0 * PHI_U_norm
    
    if(PHI_U_norm >= 0.3 and PHI_U_norm <= (5.0/6.0)):
        
        return (0.75 * PHI_U_norm) + 0.375
    
    if(PHI_U_norm > (5.0/6.0) and PHI_U_norm <= 1.0):
        
        return 1.0
    
    return PHI_U_norm










