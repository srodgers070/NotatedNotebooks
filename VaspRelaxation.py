from ase.io import read, write
from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.constraints import UnitCellFilter
from ase.io import Trajectory
from ase.optimize import BFGS
from ase import Atoms



vasptorun=read('structure.vasp')
calc=Vasp(xc='pbe')
calc.set(kspacing=0.35)
calc.set(kpar=10)
calc.set(encut=400)
calc.set(ediff=10**-6)
vasptorun.set_calculator(calc)
ucf=UnitCellFilter(vasptorun)
traj=Trajectory ('structure.traj', 'w', vasptorun)
dyn= BFGS(ucf, trajectory=traj)
dyn.run(fmax=0.01)
traj=Trajectory('structure.traj', 'r')

file=open("name.txt", "r+")
Name=file.read()
file.close()


from ase.db import connect
db=connect('/home/s376r951/Latest_Project/VASPDFT/StartingFresh/FixedCell.db')
for row in db.select(name=Name, struct_type="initial"):
    initial_id=row.id
    final_id=row.final_struct_id
    db.update(id=initial_id, atoms=vasptorun, converged=True, ASE_Calculator="VASP")
    db.update(id=final_id, atoms=vasptorun, ASE_Calculator="VASP")