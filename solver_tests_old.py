#MATEUS TEIXEIRA MAGALHAES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import numpy as np
import math
import tools
import upwind_schemes as SCHEMES
import advection_diffusion_solver_1D as solver

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
    
def TRAB1_CFD_MESTRADO():
    
    #max and minimum values for plot
    ymin = -2.0
    ymax = 2.0
    
    #domain in x direction
    domx = [0, 1]
    
    #number of points in x direction
    nx = 200
    
    #spacing discretization size
    dx = (domx[-1] - domx[0]) / (nx-1)
    
    #FSFL param
    beta = 0;
    
    #TOPUS param
    alpha = 2;
    
    #domains in time
    domts = [[0, 0.1], [0, 0.3], [0, 0.5]]
    
    #Courant numbers
    cfls = [0.1, 0.5, 0.9]
    
    #viscosity numbers
    vs = [0.001, 0.0005, 0.00025]
        
    x_axes = [domx[0] + i*dx for i in range(nx)]
    
    initial_condition = list(initial_cond_func_Burges(np.array(x_axes)))
    
    name = 'results_TRAB1_CFD_MESTRADO/initial.png'
    
    tools.save_fig(x_axes, initial_condition, name, 'Condicao inicial', 'Condicao inicial',\
                marker = None,\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = True, margin = 0.1, ymin = ymin, ymax = ymax)
    
    #-----------------------EXERCISE 2-----------------------
    
    #path to save results
    PATH = 'results_TRAB1_CFD_MESTRADO/EXE2_'
    
    #viscosity
    v = 0.001
    
    for cfl in cfls:
        
        for domt in domts:
            
            #seting upwind scheme, label and parameter
            SCHEME = SCHEMES.FSFL
            SCHEME_LABEL = 'FSFL'
            param = beta
            solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func_Burges,\
                                       initial_condition[0], initial_condition[-1], SCHEME, param,\
                                       None,\
                                       SCHEME_LABEL, marker = '.',\
                                       PATH = PATH,\
                                       equation_type = solver.Equation_types.Burges,\
                                       a = 1,\
                                       save_step_by_step = False, clean_plot = False,
                                       ymin = ymin, ymax = ymax)
            
            #seting upwind scheme, label and parameter
            SCHEME = SCHEMES.ADBQUICKEST
            SCHEME_LABEL = 'ADBQUICKEST'
            param = cfl
            solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func_Burges,\
                                       initial_condition[0], initial_condition[-1], SCHEME, param,\
                                       None,\
                                       SCHEME_LABEL, marker = None,\
                                       PATH = PATH,\
                                       equation_type = solver.Equation_types.Burges,\
                                       a = 1,\
                                       save_step_by_step = False, clean_plot = True,
                                       ymin = ymin, ymax = ymax)
    
    #--------------------------------------------------------
    
    
    #-----------------------EXERCISE 3-----------------------
    
    #path to save results
    PATH = 'results_TRAB1_CFD_MESTRADO/EXE3_'
    
    #Courant number
    cfl = 0.5
    
    for v in vs:
        
        for domt in domts:
            
            #seting upwind scheme, label and parameter
            SCHEME = SCHEMES.FSFL
            SCHEME_LABEL = 'FSFL'
            param = beta
            solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func_Burges,\
                                       initial_condition[0], initial_condition[-1], SCHEME, param,\
                                       None,\
                                       SCHEME_LABEL, marker = '.',\
                                       PATH = PATH,\
                                       equation_type = solver.Equation_types.Burges,\
                                       a = 1,\
                                       save_step_by_step = False, clean_plot = False,
                                       ymin = ymin, ymax = ymax)
            
            #seting upwind scheme, label and parameter
            SCHEME = SCHEMES.ADBQUICKEST
            SCHEME_LABEL = 'ADBQUICKEST'
            param = cfl
            solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func_Burges,\
                                       initial_condition[0], initial_condition[-1], SCHEME, param,\
                                       None,\
                                       SCHEME_LABEL, marker = None,\
                                       PATH = PATH,\
                                       equation_type = solver.Equation_types.Burges,\
                                       a = 1,\
                                       save_step_by_step = False, clean_plot = True,
                                       ymin = ymin, ymax = ymax)
    
    #--------------------------------------------------------


def parameter_comparation_boundary_layer_equation_test():
    
    ymin = -0.1
    ymax = 1.1
    
    alphas = [[-2.0, '.'], [0.0, 'x'], [2.0, '*']]
    
    betas = [[0.0, '.'], [1.0, 'x'], [2.0, '*']]
    
    gammas = [[4.0, '.'], [8.0, 'x'], [12.0, '*']]
    
    lams = [[16.0, '.'], [48.0, 'x'], [95.0, '*']]
    
    schemes = [[SCHEMES.TOPUS, 'alpha'],
               [SCHEMES.FSFL, 'beta'],
               [SCHEMES.SDPUS_C1, 'gamma'],
               [SCHEMES.EPUS, 'lambda']] 
#    schemes = [[SCHEMES.TOPUS, 'alpha']]   

    params_list = [alphas, betas, gammas, lams]
#    params_list = [alphas]
    
        
    domx = [0, 1]
    
    domt = [0, 0.5]
 
    PATH = 'camada_limite_param_comp_'
    
    vs_nxs = [[1, 10], [0.1, 80], [0.01, 640]]
#    vs_nxs = [[1, 10]]
    
    
    save_analitic = True
        
    t = domt[-1]
    
    for scheme_index in range(len(schemes)):
    
        tools.clear()
        
        scheme = schemes[scheme_index]
        
        params = params_list[scheme_index]
        
        custom_path = PATH + str(scheme[0].__name__) + '/'
        
        time_precision = 4
            
        time_string = "{:." + str(time_precision) + "f}"
        
        fileName = custom_path + 'result'\
                    + '_time=' + time_string.format(t)\
                    + '.png'
                    
        title = 'tempo = ' + str(t)
        
        for v_nx in vs_nxs:
            
            initial_cond_func = lambda x: initial_cond_func_Boundary_Layer(x)
            
            v = v_nx[0]
            
            if save_analitic:
                
                x_axes = np.linspace(domx[0], domx[-1], 640)
                
                analitic = [analitic_boundary_layer_equation(x_axes[i], v)\
                            for i in range(640)]
                
                tools.save_fig(x_axes, analitic, fileName, 'analitica', 'analitica_v=' + str(v),\
                        marker = None,\
                        xlabel = 'x', ylabel = 'y',\
                        clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax)
                        
            for _param in params:
            
                param = _param[0]
                
                marker = '.'#_param[1]
                
                SCHEME = scheme[0]
                
                SCHEME_LABEL = 'v=' + str(v) + '_' + scheme[-1] + '=' + str(param)
                
                nx = v_nx[-1]
                
                dx = (domx[-1] - domx[0]) / (nx-1)
    
                dt = 0.01 / float(nx)
                
                cfl = dt / dx
                
                solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                               initial_cond_func,\
                                               initial_cond_func(domx[0]),\
                                               1.0,\
                                               SCHEME, param,\
                                               None,\
                                               SCHEME_LABEL, marker = marker,\
                                               PATH = custom_path,\
                                               equation_type = solver.Equation_types.Boundary_layer,\
                                               a = 1,\
                                               save_step_by_step = False, clean_plot = False,\
                                               ymin = ymin, ymax = ymax,
                                               customFileName = fileName,
                                               customTitle = title,
                                               outsideLegend=True)

def schemes_comparation_boundary_layer_equation_test():
    
    tools.clear()
    
    ymin = -0.1
    ymax = 1.1
    
    alpha = 2.0
    
    beta = 2.0
    
    gamma = 12.0
    
    lam = 95.0
    
    params = [alpha, beta, gamma, lam]
    
    schemes = [[SCHEMES.TOPUS, 'TOPUS alpha'],
               [SCHEMES.FSFL, 'FSFL beta'],
               [SCHEMES.SDPUS_C1, 'SDPUS-C1 gamma'],
               [SCHEMES.EPUS, 'EPUS lambda']]
        
    domx = [0, 1]
    
    domt = [0, 0.5]
 
    PATH = 'camada_limite_schemes_comp_'
    
    vs_nxs = [[1, 10], [0.1, 80], [0.01, 640]]
#    vs_nxs = [[1, 10]]
    
    save_analitic = True
    
    for v_nx in vs_nxs:
        
        v = v_nx[0]
    
        if save_analitic:
            
            time_precision = 4
                
            time_string = "{:." + str(time_precision) + "f}"
            
            fileName = PATH + '/' + 'result'\
                        + '_time=' + time_string.format(domt[-1])\
                        + '.png'
            
            x_axes = np.linspace(domx[0], domx[-1], 640)
            
            analitic = [analitic_boundary_layer_equation(x_axes[i], v)\
                        for i in range(640)]
            
            tools.save_fig(x_axes, analitic, fileName, 'analitica', 'analitica_v=' + str(v),\
                    marker = None,\
                    xlabel = 'x', ylabel = 'y',\
                    clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax)
    
    t = domt[-1]
    
    for scheme_index in range(len(schemes)):
        
        scheme = schemes[scheme_index]
                
        custom_path = PATH + '/'
        
        fileName = custom_path + 'result'\
                    + '_time=' + time_string.format(t)\
                    + '.png'
                    
        title = 'tempo = ' + str(t)
        
        for v_nx in vs_nxs:
            
            initial_cond_func = lambda x: initial_cond_func_Boundary_Layer(x)
            
            v = v_nx[0]
                           
            param = params[scheme_index]
            
            marker = '.'#_param[1]
            
            SCHEME = scheme[0]
            
            SCHEME_LABEL = 'v=' + str(v) + ' ' + scheme[-1] + '=' + str(param)
            
            nx = v_nx[-1]
            
            dx = (domx[-1] - domx[0]) / (nx-1)

            dt = 0.01 / float(nx)
            
            cfl = dt / dx
            
            solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                           initial_cond_func,\
                                           initial_cond_func(domx[0]),\
                                           1.0,\
                                           SCHEME, param,\
                                           None,\
                                           SCHEME_LABEL, marker = marker,\
                                           PATH = custom_path,\
                                           equation_type = solver.Equation_types.Boundary_layer,\
                                           a = 1,\
                                           save_step_by_step = False, clean_plot = False,\
                                           ymin = ymin, ymax = ymax,
                                           customFileName = fileName,
                                           customTitle = title,
                                           outsideLegend=True)

        
def linear_advection_equation_test():
    
    #max and minimum values for plot
    ymin = -1.1
    ymax = 1.1
    
    #advection velocity
    a = 1
    
    case = 3
    
    cfl = 0.05
    
    nx = 400
    
    v = 0.25
    
    domx = [-1, 1]
    
    domt = [0, 0.125]
    
    dx = (domx[-1] - domx[0]) / (nx-1)
 
    PATH = 'linear_advection_CASE=' + str(case) + '/'
    
    n_v = 1
    
    #FSFL param
    beta = 0
    
    #TOPUS param
    alpha = 2;
    
    t = domt[-1]

    x_axes = [domx[0] + i*dx for i in range(nx)]

    save_analitic = False

    if save_analitic:

        analitic = [analitic_linear_advection(domx[0] + i*dx, u_0, 0, a, case = case)\
                    for i in range(nx)]
        
        fileName = 'analitic_' + PATH + 'INITIAL.png'
        tools.save_fig(x_axes, analitic, fileName, 'analitica', 'analitica',\
                marker = None,\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = True, margin = 0.1, ymin = ymin, ymax = ymax)

    
        analitic = [analitic_linear_advection(domx[0] + i*dx, u_0, t, a, case = case)\
                    for i in range(nx)]
        
        fileName = 'analitic_' + PATH + 'FINAL.png'
        tools.save_fig(x_axes, analitic, fileName, 'analitica', 'analitica',\
                marker = None,\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = True, margin = 0.1, ymin = ymin, ymax = ymax)
    
    for i in range(n_v):
        
        SCHEME = SCHEMES.TOPUS
        SCHEME_LABEL = 'TOPUS'
        param = alpha
        
        initial_cond_func = lambda x: initial_cond_func_Linear_Advection(x, case = case)
        
        analitic_sol = lambda x, t: analitic_linear_advection(x, u_0, t, a, case = case)
        
        solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                       initial_cond_func,\
                                       initial_cond_func(domx[0]),\
                                       initial_cond_func(domx[-1]),\
                                       SCHEME, param,\
                                       analitic_sol,\
                                       SCHEME_LABEL, marker = '.',\
                                       PATH = SCHEME_LABEL + '=' + str(param) + '_' + PATH,\
                                       equation_type = solver.Equation_types.Linear_advection,\
                                       a = 1,\
                                       save_step_by_step = True, clean_plot = True,\
                                       ymin = ymin, ymax = ymax)
        
        v = v / 10
    
    
def linear_advection_parameter_comparation(folderPathSubName = "", cases = [2, 3]):
    
    number_of_params = 3
    
    alphas = np.linspace(-2.0, 2.0, number_of_params)
    alphas = [float('%.2f' % el) for el in alphas]
    
    betas = np.linspace(0.0, 2.0, number_of_params)
    betas = [float('%.2f' % el) for el in betas]
    
    gammas = np.linspace(4.0, 12.0, number_of_params)
    gammas = [float('%.2f' % el) for el in gammas]
    
    lams = np.linspace(16.0, 95.0, number_of_params)
    lams = [float('%.2f' % el) for el in lams]
    
    params_list = [alphas]
#    params_list = [alphas, betas, gammas, lams]
    
    
    schemes = [[SCHEMES.TOPUS, 'alpha']]
#    schemes = [[SCHEMES.TOPUS, 'alpha'],
#               [SCHEMES.FSFL, 'beta'],
#               [SCHEMES.SDPUS_C1, 'gamma'],
#               [SCHEMES.EPUS, 'lambda']]
    
    #advection velocity
    a = 1
    
    ymins = [-1, -1.1]
    
    ymaxs = [2, 1.1]
    
    cfl = 0.9
    
    nx = 400
    
    v = 0.25
    
    domx = [-1, 1]
    
    domts = [[0.0, 0.25], [0, 0.125]]
    
    dx = (domx[-1] - domx[0]) / (nx-1)

    x_axes = [domx[0] + i*dx for i in range(nx)]

    save_analitic = True
    
    setters = ['.', '+', '3']

    for case_index in range(len(cases)):
    
        case = cases[case_index]
        
        ymin = ymins[case_index]
        
        ymax = ymaxs[case_index]
        
        domt = domts[case_index]
        
        t = domt[-1]
        
        PATH = 'param_compar_lin_ad_CASE=' + str(case) + '/'
        
        
        for scheme_index in range(len(schemes)):
            
            tools.clear()
            
            scheme = schemes[scheme_index]
            
            SCHEME = scheme[0]
            
            params = params_list[scheme_index]
            
            setter_index = -1
            for param in params:
                setter_index = setter_index + 1
                
                setter = setters[setter_index]
#                param = param_[0]
                
                marker = None#param_[-1]
                
                SCHEME_LABEL = scheme[-1] + '=' + str(param)
                
                initial_cond_func = lambda x: initial_cond_func_Linear_Advection(x, case = case)
            
#                omitPlot = not (param == params[0] or param == params[-1])
            
                custom_path = folderPathSubName + "_" + str(SCHEME.__name__) + '=' + str(param) + '_' + PATH
            
                if save_analitic:
            
                    time_precision = 4
                
                    time_string = "{:." + str(time_precision) + "f}"
                    
                    fileName = custom_path + 'result'\
                    + '_time=' + time_string.format(t)\
                    + '_v=' + str(v)\
                    + '_cfl=' + str(cfl) + '.png'
                    
                    analitic = [analitic_linear_advection(domx[0] + i*dx, u_0, t, a, case = case)\
                                for i in range(nx)]
                    
                    tools.save_fig(x_axes, analitic, fileName, 'analitica', 'analitica',\
                            marker = None,\
                            xlabel = 'x', ylabel = 'y',\
                            clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax)
            
                solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                               initial_cond_func,\
                                               initial_cond_func(domx[0]),\
                                               initial_cond_func(domx[-1]),\
                                               SCHEME, param,\
                                               None,\
                                               SCHEME_LABEL, marker = marker,\
                                               PATH = custom_path,\
                                               equation_type = solver.Equation_types.Linear_advection,\
                                               a = 1,\
                                               save_step_by_step = False, clean_plot = True,\
                                               ymin = ymin, ymax = ymax,
                                               traced = True)
        
def linear_advection_schemes_comparation():
    
    alpha = 2.0
    
    beta = 2.0
    
    gamma = 12.0
    
    lam = 95.0
    
    params = [alpha, beta, gamma, lam]
    
    schemes = [[SCHEMES.TOPUS, 'TOPUS alpha'],
               [SCHEMES.FSFL, 'FSFL beta'],
               [SCHEMES.SDPUS_C1, 'SDPUS-C1 gamma'],
               [SCHEMES.EPUS, 'EPUS lambda']]
    
    #advection velocity
    a = 1
    
    cases = [2, 3]
    
    ymins = [-0.1, -1.1]
    
    ymaxs = [1.1, 1.1]
    
    cfl = 0.05
    
    nx = 400
    
    v = 0.25
    
    domx = [-1, 1]
    
    domts = [[0.0, 0.25], [0, 0.125]]
    
    dx = (domx[-1] - domx[0]) / (nx-1)

    x_axes = [domx[0] + i*dx for i in range(nx)]

    save_analitic = True

    for case_index in range(len(cases)):
    
        tools.clear()
        
        case = cases[case_index]
        
        ymin = ymins[case_index]
        
        ymax = ymaxs[case_index]
        
        domt = domts[case_index]
        
        t = domt[-1]
        
        PATH = 'schemes_compar_lin_ad_CASE=' + str(case) + '/'
        
        for scheme_index in range(len(schemes)):
                        
            scheme = schemes[scheme_index]
            
            SCHEME = scheme[0]
            
            param = params[scheme_index]
            
            custom_path = PATH
            
            if scheme_index == 0 and save_analitic:
        
                time_precision = 4
            
                time_string = "{:." + str(time_precision) + "f}"
                
                fileName = custom_path + 'result'\
                + '_time=' + time_string.format(t)\
                + '_v=' + str(v)\
                + '_cfl=' + str(cfl) + '.png'
                
                analitic = [analitic_linear_advection(domx[0] + i*dx, u_0, t, a, case = case)\
                            for i in range(nx)]
                
                tools.save_fig(x_axes, analitic, fileName, 'analitica', 'analitica',\
                        marker = None,\
                        xlabel = 'x', ylabel = 'y',\
                        clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax)
                
            marker = None#param_[-1]
            
            SCHEME_LABEL = scheme[-1] + '=' + str(param)
            
            initial_cond_func = lambda x: initial_cond_func_Linear_Advection(x, case = case)
        
            solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                           initial_cond_func,\
                                           initial_cond_func(domx[0]),\
                                           initial_cond_func(domx[-1]),\
                                           SCHEME, param,\
                                           None,\
                                           SCHEME_LABEL, marker = marker,\
                                           PATH = custom_path,\
                                           equation_type = solver.Equation_types.Linear_advection,\
                                           a = 1,\
                                           save_step_by_step = False, clean_plot = False,\
                                           ymin = ymin, ymax = ymax)

def initial_and_final_linear_advection(ymax, ymin, t_final, case = 2, a = 1):
    
    nx = 100
    
    domx = [-1, 1]
    
    dx = (domx[-1] - domx[0]) / (nx-1)
    
    x_axes = [domx[0] + i*dx for i in range(nx)]
    
    analitic = [analitic_linear_advection(domx[0] + i*dx, u_0, 0, a, case = case)\
                            for i in range(nx)]
    
    tools.clear()
    tools.save_fig(x_axes, analitic, 'initials/START_case=' + str(case) + '.png', '', 't=' + str(0),\
            marker = None,\
            xlabel = 'x', ylabel = 'y',\
            clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax, is_traced = False)
            
    tools.clear()
    tools.save_fig(x_axes, analitic, 'initials/initial_case=' + str(case) + '.png', '', 't=' + str(0),\
            marker = None,\
            xlabel = 'x', ylabel = 'y',\
            clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax, is_traced = True)
    
    analitic = [analitic_linear_advection(domx[0] + i*dx, u_0, t_final, a, case = case)\
                            for i in range(nx)]
    
    tools.save_fig(x_axes, analitic, 'initials/initial_case=' + str(case) + '.png', '', 't=' + str(t_final),\
            marker = None,\
            xlabel = 'x', ylabel = 'y',\
            clean_plot = False, margin = 0.1, ymin = ymin, ymax = ymax, is_traced = False)

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
    
    ymins = [-0.1, -1.1]
    
    ymaxs = [1.1, 1.1]
    
    cfl = 0.9
    
    nx = 400
    
    v = 0.25
    
    domx = [-1, 1]
    
    domts = [[0.0, 0.25], [0, 0.125]]
    
    for case_index in range(len(cases)):
    
        case = cases[case_index]
        
        ymin = ymins[case_index]
        
        ymax = ymaxs[case_index]
        
        domt = domts[case_index]
                
        PATH = 'ERROR_param_compar_lin_ad_CASE=' + str(case) + '/'
        
        for scheme_index in range(len(schemes)):
            
            tools.clear()
            
            scheme = schemes[scheme_index]
            
            SCHEME = scheme[0]
            
            params = params_list[scheme_index]
            
            custom_path = str(SCHEME.__name__) + '_' + PATH
            
            analitic_sol = lambda x, t: analitic_linear_advection(x, u_0, t, a, case = case)
            
            ERRORS = []
            
            for param in params:
                                
                marker = None#param_[-1]
                
                SCHEME_LABEL = scheme[-1] + '=' + str(param)
                
                initial_cond_func = lambda x: initial_cond_func_Linear_Advection(x, case = case)
            
                error = solver.advection_difusion_equation_solver(nx, domx, domt, cfl, v,\
                                               initial_cond_func,\
                                               initial_cond_func(domx[0]),\
                                               initial_cond_func(domx[-1]),\
                                               SCHEME, param,\
                                               analitic_sol,\
                                               SCHEME_LABEL, marker = marker,\
                                               PATH = custom_path,\
                                               equation_type = solver.Equation_types.Linear_advection,\
                                               a = 1,\
                                               save_step_by_step = False, clean_plot = False,\
                                               ymin = ymin, ymax = ymax,
                                               generateError = True,
                                               errorType=errorType)
                
                ERRORS.append(error)
                
            tools.save_fig(params, ERRORS, PATH + scheme[0].__name__ + ".png", 
                           scheme[0].__name__ + " cfl = " + str(cfl),
                           None, xlabel = "parâmetro " + scheme[-1],
                           ylabel = 'Erro')
            
            



#TRAB1_CFD_MESTRADO()
#linear_advection_equation_test()
#initial_and_final_linear_advection(1.1, -0.1, 0.25, case = 2)
#initial_and_final_linear_advection(1.1, -1.1, 0.125, case = 3)
#linear_advection_parameter_comparation(folderPathSubName = 'REFINED',
#                                       cases = [2])
#linear_advection_schemes_comparation()
#parameter_comparation_boundary_layer_equation_test()
#schemes_comparation_boundary_layer_equation_test()
#linear_advection_parameter_comparation_with_error(errorType=2)
























