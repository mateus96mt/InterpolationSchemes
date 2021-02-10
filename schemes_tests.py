#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import numpy as np
import matplotlib.pyplot as plt
import upwind_schemes as SCHEMES
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

def generate_OPQ_LEONARD(color='lightgray'):
    
    x, y, = [0.0, 0.5, 1.0], [0.0, 1.0, 1.0]
    o = [0,0]
    p = [1,1]
    q = [0.5, 0.75]

    plt.figure(figsize=(8, 8))
    plt.axis('equal')
    plt.xlabel(r'$\hat{\phi}_U$', fontsize=18)
    plt.ylabel(r'$\hat{\phi}_f$', fontsize=18, rotation=0)
    plt.title('Pontos O, P e Q de Leonard')
    
    plt.plot([x[0], x[-1]], [y[0], y[-1]], color=color)
    plt.plot([x[0], x[1]], [y[0], y[1]], color=color)
    plt.plot([x[1], x[-1]], [y[1], y[-1]], color=color)
    plt.scatter(o[0], o[1], label='O', marker='D', s=70, color = 'black')
    plt.scatter(p[0], p[1], label='P', marker='X', s=70, color = 'black')
    plt.scatter(q[0], q[1], label='Q', marker='o', s=70, color = 'black')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('SCHEMES_TESTS/opq_leonard.png', dpi=100)

def generate_scheme_curve(TVD_color, colors, params, n,\
                          param_name, SCHEME, SCHEME_NAME,
                          fileName):
    
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
        
        y_axes_FSFL = [SCHEME(x, param) for x in x_axes]
        
        
        plt.plot(x_axes, y_axes_FSFL,\
                 label=str(param_name)+'='+str(param), color = color)
        
    plt.legend(loc="best")
    
    plt.savefig(fileName, dpi=100)

betas = [0.0, 1.0, 2.0]

alphas = [2.0, 2.8, 3.2]

gammas = [4.0, 8.0, 12.0]

lams = [16.0, 48.0, 95.0]

TVD_color = 'c'

colors = ['b', 'r', 'g']

n = 100

PATH = 'schemes_graphs/'

#generate_TVD_Harten()
#generate_CBC()
generate_OPQ_LEONARD()

#param_name = r'$\beta$'
#SCHEME = SCHEMES.FSFL
#SCHEME_NAME = 'FSFL'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, betas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#param_name = r'$\alpha$'
#SCHEME = SCHEMES.TOPUS
#SCHEME_NAME = 'TOPUS'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, alphas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#param_name = r'$\gamma$'
#SCHEME = SCHEMES.SDPUS_C1
#SCHEME_NAME = 'SDPUS_C1'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, gammas, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)
#
#param_name = r'$\lambda$'
#SCHEME = SCHEMES.EPUS
#SCHEME_NAME = 'EPUS'
#fileName = PATH + SCHEME_NAME + '.png'
#generate_scheme_curve(TVD_color, colors, lams, n, param_name,\
#                      SCHEME, SCHEME_NAME, fileName)













