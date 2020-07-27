from ase.io import read,write
import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import ast
import itertools

def PullandGraph(Temp):


    new_lad=np.zeros(shape=(95,4))
    count=0
    for i in range(1,6):
        #print(i)
        pathid="%sK"%Temp+"_Parallel"
        default_path="/panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh/ParallelStorage/"+"Run%s"%i #Default path where I want to store data
        home_path="/panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh/"
        dirname="{}".format(pathid)
        dirpath=os.path.join(default_path, dirname)


        os.chdir(dirpath)
        cwd=os.getcwd()
        for filename in os.listdir(cwd):
            try:
                with open(os.path.join(cwd, filename), 'r') as file:
       
                    contents=file.read()
                    dictionary=ast.literal_eval(contents)#This lets the file be read as a dictionary
                    U_Au=-1602.0827142442133/500
                    U_Cu=-1863.1513708142281/500

                    x_value=dictionary.get('Au_conc')
                    U=dictionary.get('energy')##Get's the mean energy from the dictionary
                    y_value=(U/1000-x_value*U_Au-((1-x_value)*U_Cu))
                    variance=dictionary.get('energy_var')
                    new_lad[count]=[Temp,x_value,y_value,variance]
                    count=count+1
            except IsADirectoryError:
                pass


    broken_bigboy=np.array_split(new_lad,5)


    sorted_bigboy=np.sort(broken_bigboy,axis=1)
    #print(sorted_bigboy)

    averaged_boy=np.mean(sorted_bigboy,axis=0)

    #print(averaged_boy)
    os.chdir(home_path)
    ColorValue=temp_color["%s"%Temp]
   #print()
    #print()
    #print(averaged_boy[:,1],averaged_boy[:,2])
    plt.plot(averaged_boy[:,1],averaged_boy[:,2], '%s'%ColorValue, markersize=6, linestyle="-.", label="%sK Energies"%Temp)
    plt.legend(loc='best')
    plt.xlabel('Concentration')
    plt.ylabel('Formation Energy')
   
    
    
def PullandGraphAverages(Temp):
    new_lad=np.zeros(shape=(95,4))
    count=0
    for i in range(1,6):
        pathid="%sK"%Temp+"_Parallel"
        default_path="/panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh/ParallelStorage/"+"Run%s"%i #Default 
        home_path="/panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh/"
        dirname="{}".format(pathid)
        dirpath=os.path.join(default_path, dirname)


        os.chdir(dirpath)
        cwd=os.getcwd()
        for filename in os.listdir(cwd):
            try:
                with open(os.path.join(cwd, filename), 'r') as file:
       
                    contents=file.read()
                    dictionary=ast.literal_eval(contents)#This lets the file be read as a dictionary
                    U_Au=-3204.1654284884266/1000
                    U_Cu=-3726.3027416284563/1000

                    x_value=dictionary.get('Au_conc')
                    U=dictionary.get('energy')##Get's the mean energy from the dictionary
                    y_value=(U/1000-x_value*U_Au-((1-x_value)*U_Cu))
                    variance=dictionary.get('energy_var')
                    new_lad[count]=[Temp,x_value,y_value,variance]
                    count=count+1
            except IsADirectoryError:
                pass
    broken_bigboy=np.array_split(new_lad,5)

    other_average=np.zeros(shape=(21,4))
    count=1
    tempstore=[]
    for j in range(0,19):
        tempstore_ave=[]
        tempstore_var=[]
        conc=(5+j*5)/100
        for array in broken_bigboy:
            array=array[array[:,1].argsort()]
            tempstore_ave.append(array[:,2][j])
            tempstore_var.append(array[:,3][j])
        energy_average=sum(tempstore_ave)/5
        variance_average=sum(tempstore_var)/5
        conc=(5+j*5)/100
        other_average[count]=[Temp,conc,energy_average,variance_average]
        count=count+1
    #print("your average is",other_average)
    other_average[0]=[Temp,0,0,0]
    other_average[-1]=[Temp,1,0,0]
    os.chdir(home_path)
    ColorValue=temp_color["%s"%Temp]

    
    ColorValue=temp_color["%s"%Temp]
#      plt.errorbar(other_average[:,1], other_average[:,2], yerr=other_average[:,3], fmt='-',
#                  ecolor='k', elinewidth=1, capsize=4, markersize=6)

    plt.plot(other_average[:,1],other_average[:,2], '%s'%ColorValue, markersize=6, linestyle="-",linewidth='4', label="%sK Averages"%Temp)
    plt.legend(loc='best')
    plt.xlabel('Concentration of Au')
    plt.ylabel('Average Formation Energy')
    plt.title('Average energies over 5 mc simulations')