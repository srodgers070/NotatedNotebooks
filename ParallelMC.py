from clease.montecarlo.observers import MCObserver
import numpy as np
from joblib import Parallel, delayed
from clease import Concentration
from clease import CEBulk
from clease import NewStructures
from ase.db import connect
import json
import os
import os.path
from ase.io import read,write
from clease.calculator import attach_calculator

class EnergyEvolution(MCObserver):
    """Trace the evolution of energy."""
    def __init__(self, mc, ignore_reset=False):
        self.mc = mc
        self.energies = []
        self.mean_energies = []
        self.sq_energies=[]
        MCObserver.__init__(self)
        self.name = "EnergyEvolution"
        self.ignore_reset = ignore_reset
        self.variance=[]
    def __call__(self, system_changes):
        """Append the current energy to the MC object.

        Parameters:

        system_changes: list
            System changes. See doc-string of
            `clease.montecarlo.observers.MCObserver`
        """
        quantities=self.mc.get_thermodynamic_quantities()
        self.energies.append(self.mc.current_energy + self.mc.energy_bias)
        self.mean_energies.append(self.mc.mean_energy.mean +
                                  self.mc.energy_bias)
        self.sq_energies.append((self.mc.current_energy +
                                  self.mc.energy_bias)**2)      
        variance=quantities["energy_var"]
        self.variance.append(variance)
    def reset(self):
        """Reset the history."""
        if self.ignore_reset:
            return
        self.energies = []

    def save(self, fname: str = "energy_evolution") -> None:
        """Save the energy evolution in .csv file.

        :param fname: File name without the extension (.csv)
        """
        full_fname = fname + '.csv'
        np.savetxt(full_fname, self.energies, delimiter=",")
        print(f"Energy evolution data saved to {full_fname}")
    
    def save_Mean(self,fname: str="mean_energy_evo") -> None:
        """Save the mean energy evolution in .csv file.

        :param fname: File name without the extension (.csv)
        """
        full_fname=fname+ '.csv'
        np.savetxt(full_fname, self.mean_energies, delimiter=",")
        print(f"Energy evolution data saved to {full_fname}")
        
    def save_square(self,fname: str="sqr_energy_evo") -> None:
        """Save the mean energy evolution in .csv file.

        :param fname: File name without the extension (.csv)
        """
        full_fname=fname+ '.csv'
        np.savetxt(full_fname, self.sq_energies, delimiter=",")
        print("Energy evolution data saved to {full_fname}")
        
    def save_variance(self, fname: str="Variance") -> None:
        """Compute thermodynamic quantities."""

        full_fname=fname+"Variance"+ '.csv'
        np.savetxt(full_fname, self.variance, delimiter=",")
        print(f"variance saved to {full_fname}")

        


big_boy=np.zeros(shape=(152,2))
count=0
for j in range(1,9):
    temp=100*j

    for i in range(0,19):
        conc=5+i*5
        big_boy[count]=[temp,conc]
        count=count+1

def McSim(row):
    
    Temp=int(row[:-1])
    conc_name=int(row[-1:])
    print("Temperature is %s"%Temp, ",Concentration is %s"%conc_name, "\nChecking for files and directories")


    db_name="FixedCell.db"


    pathid="%sK"%Temp+"_Parallel"
    default_path="/panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh/ParallelStorage" #Default path
    home_path="/panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh"
    dirname="{}".format(pathid)
    dirpath=os.path.join(default_path, dirname)
    os.chdir(home_path)


    conc=Concentration(basis_elements=[['Au', 'Cu']])
    settings=CEBulk(crystalstructure='fcc', a=3.8, size=[3,3,3], 
                    concentration=conc, db_name=db_name,
                    max_cluster_size=4, 
                    max_cluster_dia=[6.0,5.0,5.0])
        os.mkdir(default_path+"%s"%pathid)
        os.chdir(dirpath)
        os.mkdir(dirpath+"/MeanEnergy")
        os.mkdir(dirpath+"/EnergyEvol")
        os.mkdir(dirpath+"/Variance")
    except FileExistsError:
        #print("Directories already exist, proceeding to check if the files for this concentration exist")
        print()
    else:
        print('Directory {} created'.format(dirpath))
    
    ##Now We check if this concentration is completed. 
    os.chdir(dirpath)
    
    energyfile="%s.txt"%conc_name
    energyfilename=os.path.join(dirpath,energyfile)
    
    fileexistcheck=os.path.isfile(energyfilename)

    if (fileexistcheck):
        #print("%s"%conc_name+"%Concentration has been completed at this temperature, returning to home_path")
        os.chdir(home_path)
    else:
        #print("%s"%conc_name+"% Concentration has not been completed at this temperature, loading ECI...")
        os.chdir(home_path)
        
        
        file_name="eci_l1_fresh1.json"
        try:
            with open(file_name) as json_file:
                eci=json.load(json_file)
        except FileNotFoundError:
            #print("ECI file not found, returning to home directory and running again.")
            os.chdir(home_path)
            with open(file_name) as json_file:
                eci=json.load(json_file)
        else:
            print("ECI loaded correctly, proceeding to load structure")
       

        ##Now we read in our atoms and replace them to reach the appropriate concentration. 


        ReadIn=read('au1000.vasp')
        ReadIn=attach_calculator(settings, atoms=ReadIn, eci=eci)
        atoms_default=ReadIn
        atoms_to_replace=10*conc_name
        for j in range(0, atoms_to_replace):
                atoms_default[j].symbol='Cu'

        from clease.montecarlo import BinnedBiasPotential
        from clease.montecarlo import Montecarlo
        from clease.montecarlo import MetaDynamicsSampler

        ##We setup the monte carlo simulation for this temperature using the appropriate concentration of atoms
        mc=Montecarlo(atoms_default,Temp)
        obs=EnergyEvolution(mc)
        mc.attach(obs,interval=1000)
        energies=obs.energies
        steps=1000000


        #print("Running mc simulations...", atoms_default)
        #print()
        mc.run(steps)


        ##save all relavent quantities
        thermo=mc.get_thermodynamic_quantities()
        thermo["steps"]=steps
        os.chdir(dirpath)
        if conc_name==5:
            file=open("0%s"%conc_name+".txt", "w+")
            file.write("%s"%thermo)
            file.close()
        else:
            file=open("%s"%conc_name+".txt", "w+")
            file.write("%s"%thermo)
            file.close()
        os.chdir(dirpath+"/EnergyEvol")
        obs.save("EnergyEvol"+"_%s"%conc_name)
        os.chdir(dirpath+"/MeanEnergy")
        obs.save_Mean("MeanEnergyEvol"+"_%s"%conc_name)
        os.chdir(dirpath+"/Variance")
        obs.save_variance("Variance"+"_%s"%conc_name)
        print("%s"%conc_name+"%"+" Concentration complete at %sK"%Temp,
              "Moving to next concentration and returning to home_path")
        os.chdir(home_path)
        print()
        print()
print("setup completed, will now run parallel process")
Parallel(n_jobs=40, prefer="processes")(delayed(McSim)(row) for row in big_boy)