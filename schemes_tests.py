#MATEUS TEIXEIRA MAGALHAES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import numpy as np
import matplotlib.pyplot as plt
import upwind_schemes as SCHEMES
import flux_limiters as FLUX_LIMITERS
import os

def generate_TVD_Harten(color='b'):
    
    x, y, = [0.0, 0.5, 1.0], [0.0, 1.0, 1.0]

    plt.figure(figsize=(8, 8))
    plt.axis('equal')
    plt.xlabel(r'$\hat{\phi}_U$', fontsize=18)
    plt.ylabel(r'$\hat{\phi}_f$', fontsize=18, rotation=0)
    plt.title('Região TVD de Harten')
    plt.fill(x, y, color=color, label = 'Região TVD de Harten')
    plt.plot([x[0], x[-1]], [y[0], y[-1]], color='black')
    plt.plot([x[0], x[1]], [y[0], y[1]], color='black')
    plt.plot([x[1], x[-1]], [y[1], y[-1]], color='black')
    plt.tight_layout()
    plt.savefig('SCHEMES_TESTS/TVD_Harten.png', dpi=100)

def generate_TVD_Sweb(color='b', xmax=4.0):
    
    x, y, = [0.0, 1.0, xmax, xmax], [0.0, 2.0, 2.0, 0.0]

    plt.figure(figsize=(8, 8))
    plt.ylim(0.0, 2.5)
    plt.xlim(0.0, xmax)
    plt.xlabel(r'$r_f$', fontsize=18)
    plt.ylabel(r'$\psi$', fontsize=18, rotation=0)
    plt.title('Região TVD de Sweby')
    plt.fill(x, y, color=color, label = 'Região TVD de Sweby')
    plt.plot([x[0], x[1]], [y[0], y[1]], color='black')
    plt.plot([x[1], x[2]], [y[1], y[2]], color='black')
    
    plt.tight_layout()
    plt.savefig('SCHEMES_TESTS/TVD_Sweb.png', dpi=100)
  
    
def generate_CBC(color='b'):
    
    x, y, = [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]

    plt.figure(figsize=(8, 8))
    plt.axis('equal')
    plt.xlabel(r'$\hat{\phi}_U$', fontsize=18)
    plt.ylabel(r'$\hat{\phi}_f$', fontsize=18, rotation=0)
    plt.title('Região CBC')
    plt.fill(x, y, color=color, label = 'Região CBC')
    plt.plot([x[0], x[-1]], [y[0], y[-1]], color='black')
    plt.plot([x[0], x[1]], [y[0], y[1]], color='black')
    plt.plot([x[1], x[-1]], [y[1], y[-1]], color='black')
    plt.tight_layout()
    plt.savefig('SCHEMES_TESTS/CBC.png', dpi=100)

def generate_OPQ_LEONARD(color='lightgray', schemes = [], params = [], schemesNames = [], markers = [], n = 100):
    
    x, y, = [0.0, 0.5, 1.0], [0.0, 1.0, 1.0]
    o = [0,0]
    p = [1,1]
    q = [0.5, 0.75]

    plt.figure(figsize=(8, 8))
    plt.axis('equal')
    plt.xlabel(r'$\hat{\phi}_U$', fontsize=18)
    plt.ylabel(r'$\hat{\phi}_f$', fontsize=18, rotation=0)
#    plt.title('Pontos O, P e Q de Leonard')
    
    plt.plot([x[0], x[-1]], [y[0], y[-1]], color=color)
    plt.plot([x[0], x[1]], [y[0], y[1]], color=color)
    plt.plot([x[1], x[-1]], [y[1], y[-1]], color=color)
    plt.scatter(o[0], o[1], label='O', marker='D', s=70, color = 'black')
    plt.scatter(p[0], p[1], label='P', marker='X', s=70, color = 'black')
    plt.scatter(q[0], q[1], label='Q', marker='o', s=70, color = 'black')
#    plt.legend(loc='best')
#    plt.tight_layout()
#    
#    xmax = 1.0
#    markerpercent = 0.1
##    #create folder if not exist
##    if not os.path.exists(fileName.split('/')[0]):
##        os.makedirs(fileName.split('/')[0])
#    
    
    xmax = 1.0
    markerpercent = 0.1
    
#    plt.ylim(0.0, 2.5)
    
#    generate_TVD_Sweb(color=TVD_color, xmax = xmax)
    
    domx = [0.0, xmax]
    
    x_axes = np.linspace(domx[0], domx[1], n)
    
#    plt.title(SCHEME_NAME)
        

    for i in range(len(schemes)):
        
        SCHEME = schemes[i]
        
        param = params[i]
        
        schemeName = schemesNames[i]
        
#        color = colors[i]
        
        marker = markers[i]
        
        y_axes_FSFL = [SCHEME(x, param) for x in x_axes]
        
        
        plt.plot(x_axes, y_axes_FSFL,\
                 label=schemeName,\
                 marker = marker, markevery = int(markerpercent*n))
#            
#        
#        
#        
#    #    plt.savefig(fileName, dpi=100)
    plt.legend(loc="best")
    plt.savefig('SCHEMES_TESTS/opq_leonard.png', dpi=100)

def generate_scheme_curve(TVD_color, colors, markers, params, n,\
                          param_name, SCHEME, SCHEME_NAME,
                          fileName, markerpercent = 0.1):
    
    #create folder if not exist
    if not os.path.exists(fileName.split('/')[0]):
        os.makedirs(fileName.split('/')[0])
    
    generate_TVD_Harten(color=TVD_color)
    
    domx = [0.0, 1.0]
    
    x_axes = np.linspace(domx[0], domx[1], n)
    
    plt.title(SCHEME_NAME)
    
    for i in range(len(params)):
            
        param = params[i]
        
        color = colors[i]
        
        marker = markers[i]
        
        y_axes_FSFL = [SCHEME(x, param) for x in x_axes]
        
        
        plt.plot(x_axes, y_axes_FSFL,\
                 label=str(param_name)+'='+str(param), color = color,\
                 marker = marker, markevery = int(markerpercent*n))
        
    plt.legend(loc="best")
    
    plt.savefig(fileName, dpi=100)

def generate_scheme_curve_LIM_FLUX(TVD_color, colors, markers, params, n,\
                          param_name, SCHEME, SCHEME_NAME,
                          fileName, xmax = 4.0, markerpercent = 0.1):
    
    #create folder if not exist
    if not os.path.exists(fileName.split('/')[0]):
        os.makedirs(fileName.split('/')[0])
    
    plt.ylim(0.0, 2.5)
    
    generate_TVD_Sweb(color=TVD_color, xmax = xmax)
    
    domx = [0.0, xmax]
    
    x_axes = np.linspace(domx[0], domx[1], n)
    
    plt.title(SCHEME_NAME)
    
    for i in range(len(params)):
            
        param = params[i]
        
        color = colors[i]
        
        marker = markers[i]
        
        y_axes_FSFL = [SCHEME(x, param) for x in x_axes]
        
        
        plt.plot(x_axes, y_axes_FSFL,\
                 label=str(param_name)+'='+str(param), color = color,\
                 marker = marker, markevery = int(markerpercent*n))
        
    
    plt.legend(loc="best")
    
    plt.savefig(fileName, dpi=100)

def generate_schemeS_curve_LIM_FLUX_EXAMPLE(TVD_color, colors, markers, params, n,\
                          param_name, _SCHEMES_, SCHEME_NAMES,
                          fileName, xmax = 4.0, markerpercent = 0.1):
    
    #create folder if not exist
    if not os.path.exists(fileName.split('/')[0]):
        os.makedirs(fileName.split('/')[0])
    
    plt.ylim(0.0, 2.5)
    
    generate_TVD_Sweb(color=TVD_color, xmax = xmax)
    
    domx = [0.0, xmax]
    
    x_axes = np.linspace(domx[0], domx[1], n)
    
    for j in range(len(_SCHEMES_)):
        
        SCHEME = _SCHEMES_[j]
        
        SCHEME_NAME = SCHEME_NAMES[j]
    
        plt.title('')
        
#        for i in range(len(params)):
                
        param = params[j]
        
#        color = colors[i]
        
        marker = markers[j]
        
        y_axes_FSFL = [SCHEME(x, param) for x in x_axes]
        
        plt.plot(x_axes, y_axes_FSFL,\
                 label=SCHEME_NAME,\
                 marker = marker, markevery = int(markerpercent*n))
        
    
    plt.legend(loc="best")
    
    plt.savefig(fileName, dpi=100)


alphas = [-2.0, 0.0, 2.0]

betas = [0.0, 1.0, 2.0]

gammas = [4.0, 8.0, 12.0]

lams = [16.0, 48.0, 95.0]

TVD_color = 'c'

colors = ['b', 'r', 'g']

markers = ['+', '.', '^']


n = 100

#PATH = 'schemes_graphs_norm_var/'
##
##generate_TVD_Harten(color='c')
##generate_CBC(color='c')
##generate_OPQ_LEONARD()
#
#schemes = [SCHEMES.ADBQUICKEST, SCHEMES.CUBISTA, SCHEMES.VONOS, SCHEMES.WACEB]
##schemes = [SCHEMES.ADBQUICKEST]
#params = [0.1, 0.1, 0.1, 0.1]
##params = [0.1]
#schemeNames = ['ADBQUICKEST', 'CUBISTA', 'VONOS', 'WACEB']
##schemeNames = ['ADBQUICKEST']
#markers = ['+', 'd', '^', '.']
##markers = ['+']
#
#generate_OPQ_LEONARD(color='lightgray', 
#                     schemes = schemes, 
#                     params = params, 
#                     schemesNames = schemeNames, 
#                     markers = markers,
#                     n = n)

#param_name = r'$\beta$'
#SCHEME = SCHEMES.FSFL
#SCHEME_NAME = 'FSFL'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, markers, betas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#param_name = r'$\alpha$'
#SCHEME = SCHEMES.TOPUS
#SCHEME_NAME = 'TOPUS'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, markers, alphas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#param_name = r'$\gamma$'
#SCHEME = SCHEMES.SDPUS_C1
#SCHEME_NAME = 'SDPUS_C1'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, markers, gammas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#param_name = r'$\lambda$'
#SCHEME = SCHEMES.EPUS
#SCHEME_NAME = 'EPUS'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, markers, lams, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#PATH = 'schemes_graphs_lim_flux/'
##
#TVD_color = 'khaki'
##
#param_name = r'$\alpha$'
#SCHEME = FLUX_LIMITERS.TOPUS
#SCHEME_NAME = 'TOPUS'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve_LIM_FLUX(TVD_color, colors, markers, alphas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName, xmax = 4.0, markerpercent=0.05)
#
#param_name = r'$\beta$'
#SCHEME = FLUX_LIMITERS.FSFL
#SCHEME_NAME = 'FSFL'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve_LIM_FLUX(TVD_color, colors, markers, betas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName, xmax = 20.0, markerpercent=0.05)
#
#
#param_name = r'$\gamma$'
#SCHEME = FLUX_LIMITERS.SDPUS_C1
#SCHEME_NAME = 'SDPUS_C1'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve_LIM_FLUX(TVD_color, colors, markers, gammas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName, xmax = 3.0, markerpercent=0.05)
#
#param_name = r'$\lambda$'
#SCHEME = FLUX_LIMITERS.EPUS
#SCHEME_NAME = 'EPUS'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve_LIM_FLUX(TVD_color, colors, markers, lams, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName, xmax = 3.0, markerpercent=0.05)
#
#
#generate_TVD_Sweb(color='khaki', xmax=20.0)




PATH = 'schemes_graphs_lim_flux/'
#
TVD_color = 'khaki'
#
#param_name = r'$\cfl$'
#param_name = r'$\theta$'
#SCHEME = FLUX_LIMITERS.CUBISTA
#SCHEME_NAME = 'ADBQUICKEST'
param_names = []
_SCHEMES_ = [FLUX_LIMITERS.ADBQUICKEST, 
             FLUX_LIMITERS.CUBISTA, 
             FLUX_LIMITERS.VONOS, 
             FLUX_LIMITERS.WACEB]
#_SCHEMES_ = [FLUX_LIMITERS.TOPUS]
SCHEME_NAMES = ['ADBQUICKEST',
                'CUBISTA',
                'VONOS',
                'WACEB']
#SCHEME_NAMES = ['TOPUS']
params = [1.0, 0.25, 0.0, 0.0]
#params = [2.0]

markers = ['+', 'd', '^', '.']
fileName = PATH + 'examples_flux_lim' + '.png'
generate_schemeS_curve_LIM_FLUX_EXAMPLE(TVD_color, colors, markers, params, n, param_names,\
                      _SCHEMES_, SCHEME_NAMES, fileName, xmax = 4.0, markerpercent=0.05)

