#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import os
import matplotlib.pyplot as plt

#plot result and save in figure
def save_fig(x_axes, y_axes, fileName, title, label,\
                marker = '.',\
                xlabel = 'x', ylabel = 'y',\
                clean_plot = True, margin = 0.1, ymin = None, ymax = None):
    
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
        
    plt.plot(x_axes, y_axes, marker = marker,\
              label = label)
    
    plt.legend(loc = "best")
    plt.grid(True)
    plt.tight_layout()
    
    print('saving ', fileName)
    plt.savefig(str(fileName), dpi = 500)
    
    if clean_plot:
        
        #clean plot
        plt.cla()
        plt.clf()