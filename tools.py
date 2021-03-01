#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import os
import matplotlib.pyplot as plt
import numpy as np

#plot result and save in figure
def save_fig(x_axes, y_axes, fileName, title, label,\
                marker = '.',\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = True, margin = 0.1, ymin = None, ymax = None, 
                is_traced = False,
                outsideLegend = False, 
                setter=None):
    
    min_y_axes = min(y_axes)
    max_y_axes = max(y_axes)
    
    ysize =  max_y_axes - min_y_axes
    
    if ymin == None:
        
        ymin = min_y_axes - margin*ysize
    
    if ymax == None:
    
        ymax = max_y_axes + margin*ysize
    
    #create folder if not exist
    if not os.path.exists(fileName.split('/')[0]):
        os.makedirs(fileName.split('/')[0])

    plt.title(title)
    
    plt.ylim([ymin, ymax])
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    if is_traced:  
        
        plt.plot(x_axes, y_axes, marker = marker,\
                  label = label, dashes = [int(len(x_axes)/10)])
    else:
        
        if setter == None:
        
            plt.plot(x_axes, y_axes, marker = marker,\
                      label = label)
            
        else:
            plt.plot(x_axes, y_axes, setter,\
                      label = label)
            
        
    if outsideLegend:
        plt.legend(loc = "best")
    else:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.grid(True)
    plt.tight_layout()
    
    print('saving ', fileName)
    plt.savefig(str(fileName), dpi = 500)
    
    if clean_plot:
        
        #clean plot
        plt.cla()
        plt.clf()
        
def calculateError(exata, y_numerica, nx, domx, dx, time, tipo = 2):
    
    if exata == None:
    
        return 0.0
        
    else:
      
        x = np.linspace(domx[0], domx[1], nx)
        
        y_exata = [exata(x[i], time) for i in range(nx)]
            
        sum_num = 0.0
        sum_dem = 0.0
        max_error = max(np.array(y_exata)[:]-np.array(y_numerica)[:])
        max_exata = max(np.array(y_exata))
        
        for i in range(nx):
        
            if tipo == 1:
            
                sum_num = sum_num + (y_exata[i] - y_numerica[i])
                sum_dem = sum_dem + y_exata[i]
            
            if tipo == 2:
            
                sum_num = sum_num + ((y_exata[i] - y_numerica[i])**2)
                sum_dem = sum_dem + (y_exata[i]**2)
        
        if tipo == 1:
            
            return (sum_num/sum_dem)
        
        if tipo == 2:
            
            return np.sqrt(sum_num/sum_dem)
        
        if tipo == 3:
            
            return max_error/max_exata
    
def readOutPut(fileName):
    
    x = []
    y_exata = []
    y_num = []
    
    file = open(fileName, 'r')
    
    str_file = file.readlines()
    
    for i in range(1, len(str_file)):
    
        str_lines = str_file[i].split('\n')
        str_lines = str_lines[0].split(' ')
        
        x.append(float(str_lines[0]))
        y_exata.append(float(str_lines[1]))
        y_num.append(float(str_lines[2]))
        
    return x, y_exata, y_num
    
    
def clear():
    
    plt.cla()
    plt.clf()
