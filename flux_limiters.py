#MATEUS TEIXEIRA MAGALHAES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git


#TOPUS flux limiter
def TOPUS(rf, alpha):
    
    if rf < 0:
        
        return 0;
    
    else:
    
        return ( ( (-0.5*alpha + 1.0)*(rf**3) ) +\
               ( (alpha + 4.0)*(rf**2) ) +\
               ((-0.5*alpha+3.0)*rf) ) / ((1 + rf)**3)

#FSFL flux limiter
def FSFL(rf, beta):
        
    if rf < 0:
        
        return 0;
    
    else:
    
        return ( (beta*(rf**3)) +\
                 ((8.0-2.0*beta)*(rf**2)) +\
                 (beta*rf)
                ) / ((1 + rf)**3)

#EPUS flux limiter 
def EPUS(rf, lam):
    
    if rf < 0:
        
        return 0;
    
    else:
    
        return ( ((2*lam-32.0)*(rf**5)) +\
                 ((160.0 - 4.0*lam)*(rf**4)) +\
                 (2*lam*(rf**3))
                ) / ((1 + rf)**7)

#SDPUS-C1 flux limiter
def SDPUS_C1(rf, gamma):
    
    if rf < 0:
        
        return 0;
    
    else:
    
        return ( (4.0*(rf**3)) +\
                 (12*(rf**2))
                ) / ((1 + rf)**4)
 
