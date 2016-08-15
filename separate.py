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

import os
import sys
import shutil


xas = 'yes'
mol = 'yes'
out = 'yes'
rst = 'keep'
xry = 'keep'
xryScripts = 'delete'
batch = 'delete'
inputfiles = 'yes'
reports = 'delete'

xasfolder = 'XAS_files'
xasfiles ='XAS_files//c1.xas' 
molFolder = 'Molden_files'
molFiles = 'Molden_files//c1.mol'
outputfolder = 'Output_files'
outputfiles = 'Output_files//c1.out'
rstfolder = 'rst_files'
rstfiles = 'rst_files//c1.rst'
xryfolder = 'XRY_files'
xryfiles = 'XRY_files//c1.xry'
inputfolder = 'Input_files'
inputfiles = 'Input_files//c1.inp'
crystal = 'Crystal_file'


ff = 'Inputs'
Sep= 'Separated_files'




os.mkdir(Sep)
os.mkdir (xasfolder)
os.mkdir(molFolder)
os.mkdir (outputfolder)
os.mkdir(rstfolder)
os.mkdir (xryfolder)
os.mkdir(inputfolder)


n_carbon = 0
c_counter = 1


with open('Sources//molecule.xyz','r') as file:
    for line in file:
        if line[0] == 'C':
            n_carbon = n_carbon + 1
n_carbon = n_carbon + 1



#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |******************************To make Xas folder including xas files*************************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|
#
for i in range(1,n_carbon):
    if os.path.exists(ff+'//c'+str(i)+'//c'+str(i)+'.xas'):
        shutil.move(ff+'//c'+str(i)+'//c'+str(i)+'.xas', xasfolder)
    else:
        print 'xas fils are copied'
        
if os.path.exists(xasfolder):
    shutil.move(xasfolder, Sep)
    
    
   
    
    

    

#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |******************************To make mol folder including mol files*************************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|

 

for i in range(1,n_carbon):
    if os.path.exists(ff+'//c'+str(i)+'//c'+str(i)+'.mol'):
        shutil.move(ff+'//c'+str(i)+'//c'+str(i)+'.mol', molFolder)
    else:
        print 'mol fils are copied'
        
if os.path.exists(molFolder):
    shutil.move(molFolder, Sep)

 

#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |******************************To make output folder including output files*******************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|



for i in range(1,n_carbon):
    if os.path.exists(ff+'//c'+str(i)+'//c'+str(i)+'.out'):
        shutil.move(ff+'//c'+str(i)+'//c'+str(i)+'.out', outputfolder)
    else:
        print 'out fils are copied'
        
if os.path.exists(outputfolder): 
    shutil.move(outputfolder, Sep)

#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |******************************To make input folder including input files*********************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|

 



for i in range(1,n_carbon):
    if os.path.exists(ff+'//c'+str(i)+'//c'+str(i)+'.inp'):
        shutil.move(ff+'//c'+str(i)+'//c'+str(i)+'.inp', inputfolder)
    else:
        print 'inp fils are copied'
        
if os.path.exists(inputfolder):
    shutil.move(inputfolder, Sep)

#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |******************************************To decide about rst files**************************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|

 


for i in range(1,n_carbon):
    if os.path.exists(ff+'//c'+str(i)+'//c'+str(i)+'.rst'):
        shutil.move(ff+'//c'+str(i)+'//c'+str(i)+'.rst', rstfolder)
    else:
        print 'rst fils are copied'
            
if os.path.exists(rstfolder):
    shutil.move(rstfolder, Sep)
    
#
#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |******************************************To decide about xry files**************************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|




for i in range(1,n_carbon):
    if os.path.exists(ff+'//c'+str(i)+'//c'+str(i)+'.xry'):
        shutil.move(ff+'//c'+str(i)+'//c'+str(i)+'.xry', xryfolder)
    else:
        print 'xry fils are copied'
if os.path.exists(xryfolder):
    shutil.move(xryfolder, Sep)
        

#                                              ____________________________________________________________________________________________________________________
#                                             |    __________________________________________________________________________________________________________     |
#                                             |    |*********************************To decide about all unnecessary files***********************************|    |
#                                             |    |_________________________________________________________________________________________________________|    |
#                                             |___________________________________________________________________________________________________________________|


if batch =='delete' and xryScripts =='delete' and reports =='delete':
    shutil.rmtree(ff)
else:
    if batch =='delete':
        for i in range(1,n_carbon):
            os.remove(ff+'//c'+str(i)+'//d2k-batch.run')
    elif batch =='keep':
        print 'Batch files will be kept safe'
    if xryScripts =='delete':
        for i in range(1,n_carbon):
            os.remove(ff+'//c'+str(i)+'//xray.inp')
            os.remove(ff+'//c'+str(i)+'//xray2k.run')
            os.remove(ff+'//c'+str(i)+'//xray2k.x')
    elif batch =='keep':
        print 'X-ray script files will be kept safe'
    

            