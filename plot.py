'''
 ____________________________________________________________________________________________________________________
|    __________________________________________________________________________________________________________     |
|    |                                   @author: Sadegh Shokatian                                             |    |
|    |                                Supervisor: Stephen G. Urquhart                                          |    |
|    |                                    department of Chemistry                                              |    |
|    |                                   University of Saskatchewan                                            |    |
|    |_________________________________________________________________________________________________________|    |
|___________________________________________________________________________________________________________________|
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import re
import os
import sys
import shutil


###
steps = 40.0
ff = 'Separated_files'
Inp = 'Inputs'

n_carbon = 0
c_counter = 1
xasfiles = 'xas'

###


#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |**************************To calculate numbers of carbons in the system**********************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
#

if n_carbon == 0:
    with open('Sources//molecule.xyz','r') as file:
        for line in file:
            if line[0] == 'C':
                n_carbon = n_carbon + 1
    n_carbon = n_carbon + 1
else:
    n_carbon = n_carbon
###

###

data2 = []        #       Empty array to put all the XAD data in it
n1 = int(10)
###


#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |********This class is defined for XAS files that are not separated from the others and*******************|    |
#                                             |    |**************************** still stored in the Inputs directory****************************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
#


class Load_Inp:
    def __init__(self,i):
        plt.subplot(2,1,2)
        plt.ylabel('a. u.')
        plt.xlabel('eV')
        plt.title('Individual spectra (C1 is the spectrum at the bottom)')
        
        for i in range(1,n_carbon):
            data = np.loadtxt(Inp+'//c'+str(i)+'//c'+str(i)+'.xas')
            plt.plot(data[:,0],(data[:,1]+ 1 +  (i*steps)))
        
                
        for i in range(1,n_carbon):
            
            if i == 1:
                with open(Inp+'//c'+str(1)+'//c'+str(1)+'.xas','r') as f:
                    for line in f:
                        parts = re.split("[\s]+", line);
                        data2.append([parts[1],float(parts[2])])
            else:
                j = 0
                with open(Inp+'//c'+str(i)+'//c'+str(i)+'.xas','r') as f:
                    for line in f:
                        parts = re.split("[\s]+", line);
                        data2[j][1]+=float(parts[2])
                        j+=1            
        
        with open(Inp+'//Final.xas','w') as fout: 
            for line in data2:
                fout.write('       '+line[0]+'       '+str(line[1])+'\n')
            
        data = np.loadtxt(Inp+'//Final.xas')
        plt.subplot(2,1,1)
        plt.title('Final Spectrum')
        plt.plot(data[:,0],(data[:,1]+ 10 + n_carbon*steps))
        plt.savefig(Inp+'//sadegh.png')
        plt.ylabel('a. u.')
        #plt.xlabel('eV')
        plt.axis()
        return plt.show()
###



#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |************This class is defined for XAS files that are separated from the others and*******************|    |
#                                             |    |**************************** stored in the Separated_files directory*************************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
class Load_Sep:
    def __init__(self,i):
        #plt.subplot(2,1,2)
        plt.ylabel('a. u.')
        plt.xlabel('eV')
        
        for i in range(1,n_carbon):
            plt.title('Individual spectra (C1 is the spectrum at the bottom)')
            data = np.loadtxt(ff+'//XAS_files'+'//c'+str(i)+'.xas')
            plt.xlim(287,292)
            plt.subplot(2,1,2)
            plt.plot(data[:,0],(data[:,1]+ 1 + (i*steps)))
        
                
        for i in range(1,n_carbon):
            
            if i == 1:
                with open(ff+'//XAS_files'+'//c'+str(1)+'.xas','r') as f:
                    for line in f:
                        parts = re.split("[\s]+", line);
                        data2.append([parts[1],float(parts[2])])
            else:
                j = 0
                with open(ff+'//XAS_files'+'//c'+str(i)+'.xas','r') as f:
                    for line in f:
                        parts = re.split("[\s]+", line);
                        data2[j][1]+=float(parts[2])
                        j+=1            
        
        with open(ff+'//XAS_files'+'//Final.xas','w') as fout: 
            for line in data2:
                fout.write('       '+line[0]+'       '+str(line[1])+'\n')
            
        data = np.loadtxt(ff+'//XAS_files'+'//Final.xas')
        plt.subplot(2,1,1)
        plt.xlim(287,292)
        plt.title('Final Spectrum')
        plt.plot(data[:,0],((data[:,1]+ 10 + n_carbon*steps)))
        plt.savefig(ff+'//XAS_files'+'//sadegh.png')
        plt.ylabel('a. u.')
        #plt.xlabel('eV')
        plt.axis()
        return plt.show()    
###



#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |*****To make some changes in the XAS files which makes them readable for python(By changing D to e)******|    |
#                                             |    |************** This part is just for XAS files which are stored in the Inputs directory******************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
#               
if os.path.exists(Inp+'//c'+str(1)+'//c'+str(1)+'.xas'):
    for i in range(1,n_carbon):
        with open(Inp+'//c'+str(i)+'//c'+str(i)+'.xas','r+') as file1:
            with open(Inp+'//c'+str(i)+'//c'+str(i)+'-1.xas','w') as file2:
                for line in file1:
                    file2.write(line.replace('D','e'))
    file1.close()
    file2.close()
    for i in range(1,n_carbon):
        with open(Inp+'//c'+str(i)+'//c'+str(i)+'-1.xas','r+') as file1:
            with open(Inp+'//c'+str(i)+'//c'+str(i)+'.xas','w') as file2:
                for line in file1:
                    file2.write(line)
    file1.close()
    file2.close()
    for i in range(1,n_carbon):
        os.remove(Inp+'//c'+str(i)+'//c'+str(i)+'-1.xas')

###

#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |*****To make some changes in the XAS files which makes them readable for python(By changing D to e)******|    |
#                                             |    |********** This part is just for XAS files which are stored in the Separated_files directory*************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
if os.path.exists(ff+'//XAS_files'+'//c'+str(1)+'.xas'):
    for i in range(1,n_carbon):
        with open(ff+'//XAS_files'+'//c'+str(i)+'.xas','r+') as file1:
            with open(ff+'//XAS_files'+'//c'+str(i)+'1.xas','w') as file2:
                for line in file1:
                    file2.write(line.replace('D','e'))
    file1.close()
    file2.close()
    for i in range(1,n_carbon):
        with open(ff+'//XAS_files'+'//c'+str(i)+'1.xas','r+') as file1:
            with open(ff+'//XAS_files'+'//c'+str(i)+'.xas','w') as file2:
                for line in file1:
                    file2.write(line)
    file1.close()
    file2.close()

    for i in range(1,n_carbon):
        os.remove(ff+'//XAS_files'+'//c'+str(i)+'1.xas')

    



#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |*******************To Load one of the classes base on the adress of XAS files****************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
#

if os.path.exists(Inp+'//c'+str(1)+'//c'+str(1)+'.xas'):
    Load_Inp(i)
else:
    Load_Sep(i)


