#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import numpy as np
import math
import tools
import upwind_schemes as SCHEMES
import advection_diffusion_solver_1D as solver
import matplotlib.pyplot as plt
import os

#-----------------------------------INITIAL CONDITIONS------------------------

#INICIAL CONDITION FOR LINEAR-ADVECTION, (DISSERTACAO RAFAEL, CASOS 1, 2 e 3)
def u_0(x, case = 1):
    
    if case == 1:
    
        if x>=0 and x<0.2:
            
            return np.exp(-math.log(((x-0.15)/0.05)**2, 50))
        
        if x>0.3 and x<0.4:
            
            return 1
        
        if x>0.5 and x<0.55:
            
            return 20.0*x-10.0
        
        if x>=0.55 and x<0.66:
            
            return -20*x + 12
        
        if x>0.7 and x<0.8:
            
            return math.sqrt(1 - (((x-0.75)/0.05)**2))
        
        return 0

    if case == 2:
        
        if x>=0 and x<=0.2:
            
            return 1
        
        if x>0.2 and x<=0.4:
            
            return 4*x-0.6
        
        if x>0.4 and x<=0.6:
            
            return -4*x + 2.6
        
        if x>0.6 and x<=0.8:
            
            return 1
        
        return 0
    
    if case == 3:
        
        if x>=-1 and x<=-1/3:
            
            return -x*np.sin(3*np.pi*(x**2)/2.0)
            
        if x>-1/3 and x<1/3:
            
            return abs(np.sin(2*np.pi*x))
        
        if x>=1/3 and x<=1:
            
            return 2*x - 1.0 - (1/6)*np.sin(3*np.pi*x)
        
        return 0
        
def initial_cond_func_Burges(x):
    
    return np.sin(2 * math.pi * x)
            

def initial_cond_func_Boundary_Layer(x):
    
    return 0.0
            
def initial_cond_func_Linear_Advection(x, case = 1):
    
    return  u_0(x, case = case)

#-----------------------------------INITIAL CONDITIONS------------------------



#-----------------------------------ANALITIC SOLUTIONS------------------------

def analitic_boundary_layer_equation(x, v):
    
    return (1 - np.exp(x / v)) / (1 - np.exp(1 / v))


def analitic_linear_advection(x, u_0, t, a, case = 1):
    
    return u_0(x - a*t, case)

#-----------------------------------------------------------------------------
    
def linear_advection_parameter_comparation(cases = [2, 3]):
    
    folderName = "RESULTADOS_DATA"
    
    number_of_params = 3
    
    alphas = np.linspace(-2.0, 2.0, number_of_params)
    alphas = [float('%.2f' % el) for el in alphas]
    
    betas = np.linspace(0.0, 2.0, number_of_params)
    betas = [float('%.2f' % el) for el in betas]
    
    gammas = np.linspace(4.0, 12.0, number_of_params)
    gammas = [float('%.2f' % el) for el in gammas]
    
    lams = np.linspace(16.0, 95.0, number_of_params)
    lams = [float('%.2f' % el) for el in lams]
    
#    params_list = [alphas]
    params_list = [alphas, betas, gammas, lams]
    
    
#    schemes = [[SCHEMES.TOPUS, 'alpha']]
    schemes = [[SCHEMES.TOPUS, 'alpha'],
               [SCHEMES.FSFL, 'beta'],
               [SCHEMES.SDPUS_C1, 'gamma'],
               [SCHEMES.EPUS, 'lambda']]
    
    #advection velocity
    a = 1
    
    cfl = 0.05
    
    nx = 400
    
    v = 0.25
    
    domx = [-1, 1]
    
    domts = [[0.0, 0.25], [0, 0.125]]

    for case_index in range(len(cases)):
    
        case = cases[case_index]
        
        domt = domts[case_index]
                
        for scheme_index in range(len(schemes)):
            
            tools.clear()
            
            scheme = schemes[scheme_index]
            
            SCHEME = scheme[0]
            
            params = params_list[scheme_index]
            
            setter_index = -1
            for param in params:
                setter_index = setter_index + 1
                
                initial_cond_func = lambda x: initial_cond_func_Linear_Advection(x, case = case)
                                    
                analitic_sol_func = lambda x, t: analitic_linear_advection(x, u_0, t, a, case = case)
            
                SCHEME_LABEL = SCHEME.__name__ + "_" + scheme[-1] + '=' + str(param)
                folderName = "RESULTADOS_DATA_" + SCHEME_LABEL
            
                solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func,\
                                       initial_cond_func(domx[0]), initial_cond_func(domx[-1]), SCHEME, param,\
                                       analitic_sol_func,\
                                       folderName,\
                                       equation_type = solver.Equation_types.Linear_advection,\
                                       a = 1,\
                                       step_interval = 0.05,\
                                       errorType = 2)

def param_compar_graficos_TOPUS():
             
    name = "RESULTADOS_DATA_TOPUS_alpha="
    
    params = np.linspace(-2.0, 2.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\alpha=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\alpha=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\alpha=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema TOPUS')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(_input.split('/')[0]+'/'+_input.split('/')[0] + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def param_compar_graficos_FSFL():
              
    name = "RESULTADOS_DATA_FSFL_beta="
    
    params = np.linspace(0.0, 2.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\beta=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\beta=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\beta=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema FSFL')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(_input.split('/')[0]+'/'+_input.split('/')[0] + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def param_compar_graficos_EPUS():
              
    name = "RESULTADOS_DATA_EPUS_lambda="
    
    params = np.linspace(16.0, 95.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\lambda=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\lambda=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\lambda=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema EPUS')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(_input.split('/')[0]+'/'+_input.split('/')[0] + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def param_compar_graficos_SDPUS():
              
    name = "RESULTADOS_DATA_SDPUS_C1_gamma="
    
    params = np.linspace(4.0, 12.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\gamma=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\gamma=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\gamma=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema SDPUS-C1')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(_input.split('/')[0]+'/'+_input.split('/')[0] + '.png', dpi = 200)
    plt.cla()
    plt.clf()

def init_cond():
              
    name = "RESULTADOS_DATA_SDPUS_C1_gamma="
    
    params = np.linspace(4.0, 12.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/0.data"
    x, y_exata, y_num = tools.readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('Condição inicial')
    
    plt.savefig('caso2.png', dpi = 200)
    plt.cla()
    plt.clf()

def linear_advection_parameter_comparation_with_error(errorType=2):
    
    number_of_params = 8
    
    alphas = np.linspace(-2.0, 2.0, number_of_params)
    
    betas = np.linspace(0.0, 2.0, number_of_params)
    
    gammas = np.linspace(4.0, 12.0, number_of_params)
    
    lams = np.linspace(16.0, 95.0, number_of_params)
    
#    params_list = [alphas]
    params_list = [alphas, betas, gammas, lams]
    
    
#    schemes = [[SCHEMES.TOPUS, 'alpha']]
    schemes = [[SCHEMES.TOPUS, 'alpha'],
               [SCHEMES.FSFL, 'beta'],
               [SCHEMES.SDPUS_C1, 'gamma'],
               [SCHEMES.EPUS, 'lambda']]
    
    #advection velocity
    a = 1
    
    cases = [2]    
    
    cfl = 0.05
    
    nx = 400
    
    v = 0.25
    
    domx = [-1, 1]
    
    domts = [[0.0, 0.25], [0, 0.125]]
    
    for case_index in range(len(cases)):
    
        case = cases[case_index]
        
        domt = domts[case_index]
                
        folderName = 'ERROR_param_compar_lin_ad_CASE=' + str(case)
        
        for scheme_index in range(len(schemes)):
            
            tools.clear()
            
            scheme = schemes[scheme_index]
            
            SCHEME = scheme[0]
            
            params = params_list[scheme_index]
                        
            analitic_sol = lambda x, t: analitic_linear_advection(x, u_0, t, a, case = case)
            
            ERRORS = []
            
            for param in params:
                                                                
                initial_cond_func = lambda x: initial_cond_func_Linear_Advection(x, case = case)
            
                error = solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func,\
                                       initial_cond_func(domx[0]), 
                                       initial_cond_func(domx[-1]),
                                       SCHEME, param,\
                                       analitic_sol,\
                                       folderName,\
                                       equation_type = solver.Equation_types.Linear_advection,\
                                       a = 1,\
                                       step_interval = 1,\
                                       errorType = 2)
                
                ERRORS.append(error)
                
            plt.grid(True, linestyle='--')
            plt.title('Esquema ' + SCHEME.__name__)
            plt.xlabel(r'$' + '\\'+ scheme[1] + '$')
            plt.ylabel('Erro')
            plt.ylim([0.1575, 0.1750])
            plt.tight_layout()
            plt.plot(params, ERRORS, marker='*')
            plt.savefig(folderName + '/' + SCHEME.__name__ + 'cfl=' + str(cfl) + '.png', dpi=200)
            plt.cla()
            plt.clf()

linear_advection_parameter_comparation(cases = [2])
param_compar_graficos_TOPUS()
param_compar_graficos_FSFL()
param_compar_graficos_SDPUS()
param_compar_graficos_EPUS()
#init_cond()
#linear_advection_parameter_comparation_with_error(errorType=2)




















