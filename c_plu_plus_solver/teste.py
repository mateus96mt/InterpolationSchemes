#MATEUS TEIXEIRA MAGALHÃƒES - UFJF - PGMC
#github repository: https://github.com/mateus96mt/InterpolationSchemes.git

import matplotlib.pyplot as plt

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


name = "OUTPUT"

param = 2.0

sub_path = name
_input = sub_path + "/FINAL.data"
x, y_exata, y_num = readOutPut(_input)

plt.plot(x, y_exata, linewidth=1.0, color='red', label="Solução exata")
sub_path = name
_input = sub_path + "/FINAL.data"
x, y_exata, y_num = readOutPut(_input)
plt.plot(x, y_num, 'x', color='green', ms=2.5, 
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
plt.show()



















 
