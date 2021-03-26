#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import matplotlib.pyplot as plt
import numpy as np

def readOutPut(fileName):
    
    x = []
    y_exata = []
    y_num = []
    
    file = open(fileName, 'r')
    
    str_file = file.readlines()
    
    for i in range(2, len(str_file)):
    
        str_lines = str_file[i].split('\n')
        str_lines = str_lines[0].split(' ')
        
        x.append(float(str_lines[0]))
        y_exata.append(float(str_lines[1]))
        y_num.append(float(str_lines[2]))
        
    return x, y_exata, y_num             

def solutionOverTime():
    
    n = 56

    schemes = ["TOPUS", "FSFL", "SDPUS_C1", "EPUS"]
    params_name = ["\\alpha", "\\beta", "\\gamma", "\\lambda"]
    
    alphas = [-2.0, 0.0, 2.0]
    betas = [0.0, 1.0, 2.0]
    gammas = [4.0, 8.0, 12.0]
    lambdas = [16.0, 55.5, 95.0]
    
    params = [alphas, betas, gammas, lambdas]
    
    for scheme_index in range(len(schemes)):
        
        param = params[scheme_index]
        
        for param_value in param:
            
            name = schemes[scheme_index] + "_" + str(param_value)
            
            print(name)
            
            for i in range(n):
                                
                sub_path = name
                _input = sub_path + "/" + str(i) + ".data"
                x, y_exata, y_num = readOutPut(_input)
                
                plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
                sub_path = name
                _input = sub_path + "/" + str(i) + ".data"
                x, y_exata, y_num = readOutPut(_input)
                plt.plot(x, y_num, 'x', color='green', ms=2.5, 
                         label=r'$' + params_name[scheme_index] +  '=' + str(param_value) + '$')
                
                plt.grid(True, linestyle='--')
                plt.title('Esquema ' + schemes[scheme_index])
                plt.legend(loc="best")
                plt.xlabel(r'$\mathrm{x}$')
                plt.ylabel(r'$\mathrm{u}$')
                plt.ylim([-1.2, 1.2])
                plt.tight_layout()
                
                plt.savefig(_input.split('/')[0]+'/'+_input.split('/')[1] + '.png', dpi = 200)
                plt.cla()
                plt.clf()
                plt.show()


def param_compar_graficos_TOPUS():
             
    name = "TOPUS_"
    
    params = np.linspace(-2.0, 2.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num =readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\alpha=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\alpha=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\alpha=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema TOPUS')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(name + '_param_comp' + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def param_compar_graficos_FSFL():
              
    name = "FSFL_"
    
    params = np.linspace(0.0, 2.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\beta=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\beta=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\beta=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema FSFL')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(name + '_param_comp' + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def param_compar_graficos_EPUS():
              
    name = "EPUS_"
    
    params = np.linspace(16.0, 95.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\lambda=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\lambda=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\lambda=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema EPUS')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(name + '_param_comp' + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def param_compar_graficos_SDPUS():
              
    name = "SDPUS_C1_"
    
    params = np.linspace(4.0, 12.0, 3)
    param = params[0]
    
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    
    plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
    
    param = params[0]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'x', color='green', ms=2.5, 
             label=r'$\gamma=' + str(param) + '$')
    
    param = params[1]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, '.', color='blue', ms=2.5, 
             label=r'$\gamma=' + str(param) + '$')
    
    param = params[2]
    sub_path = name + str(param)
    _input = sub_path + "/FINAL.data"
    x, y_exata, y_num = readOutPut(_input)
    plt.plot(x, y_num, 'v', color='tab:blue', ms=2.5, 
             label=r'$\gamma=' + str(param) + '$')
    
    
    plt.grid(True, linestyle='--')
    plt.title('Esquema SDPUS-C1')
    plt.legend(loc="best")
    plt.xlabel(r'$\mathrm{x}$')
    plt.ylabel(r'$\mathrm{u}$')
    plt.ylim([-1.2, 1.2])
    plt.tight_layout()
    
    plt.savefig(name + '_param_comp' + '.png', dpi = 200)
    plt.cla()
    plt.clf()
    
def readErrors(fileName, folderName, cfl, n = 8):
    
    file = open(fileName, 'r')
    
    str_file_ = file.readlines()
    
    str_file = []
    
    for i in range(len(str_file_)):
        
        str_file.append(float(str_file_[i].split('\n')[0]))
    
    alphas = np.linspace(-2.0, 2.0, n)
    betas = np.linspace(0.0, 2.0, n)
    gammas = np.linspace(4.0, 12.0, n)
    lams = np.linspace(16.0, 95.0, n)
    
    params = [alphas, betas, gammas, lams]
    
    erro_topus = str_file[0:n]
    erro_fsfl = str_file[n:2*n]
    erro_sdpus = str_file[2*n:3*n]
    erro_epus = str_file[3*n:4*n]
    
    erros = [erro_topus, erro_fsfl, erro_sdpus, erro_epus]
    
    schemes = ["TOPUS", "FSFL", "SDPUS-C1", "EPUS"]
    params_names = ["alpha", "beta", "gamma", "lambda"]
    
    ylims = [min(str_file)-0.1*abs(max(str_file)), max(str_file)+0.1*abs(max(str_file))]
    
    for scheme_index in range(len(schemes)):
    
        print(erros[scheme_index])
        
        param = params[scheme_index]
        
        plt.grid(True, linestyle='--')
        plt.title('Esquema ' + schemes[scheme_index])
        plt.xlabel(r'$' + '\\'+ params_names[scheme_index] + '$')
        plt.ylabel(r'$\Vert E  \Vert_2$', rotation=0, labelpad=20)
        plt.ylim([ylims[0], ylims[1]])
        plt.tight_layout()
        plt.plot(param, erros[scheme_index], marker='*')
        plt.savefig(folderName + '/' + schemes[scheme_index] + 'cfl=' + str(cfl) + '.png', dpi=200)
        plt.cla()
        plt.clf()
    
    
    
#solutionOverTime()
    
#param_compar_graficos_TOPUS()
#param_compar_graficos_FSFL()
#param_compar_graficos_SDPUS()
#param_compar_graficos_EPUS()
    
#readErrors("erros_cfl_0.05.txt", "ERROS_cfl=0.05", 0.9, n=8)















 
