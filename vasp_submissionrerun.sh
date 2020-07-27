#!/bin/bash
#SBATCH --job-name=job_auto    # Job name
#SBATCH --partition=sixhour           # Partition Name (Required)
##SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, A$
##SBATCH --mail-user=email@ku.edu      # Where to send mail
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1                    # Run on a single CPU
##SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=06:00:00             # Time limit days-hrs:min:sec

module load anaconda
conda activate Research

~/.conda/envs/Research/bin/python  /panfs/pfs.local/home/s376r951/Latest_Project/VASPDFT/StartingFresh/VaspRelaxationRerun.py
