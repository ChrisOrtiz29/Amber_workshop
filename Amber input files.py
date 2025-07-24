#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os, sys


# In[ ]:


#Energy Minimization Stage 1
with open("min0.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Energy Minimization, Stage1, RESTRAINING ALL HEAVY ATOMS EXCEPT WATER AND IONS,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=1, !Energy Minimization
  maxcyc=5000, !Total Minimization Cycles to be run. Steepest Decent First, then Conjugate Gradient Method if ncyc < maxcyc
  ncyc=1000, !Number of Steepest Decent Minimization Steps to run before switching to Conjugate Gradient
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=1, !Complete Interactions are Calculated
  ntc=1, !SHAKE is NOT performed, DEFAULT
  ntpr=10, !Every 10 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=100, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=100, !Every 100 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="!(@H=|:WAT|@Na+|@Cl-)", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=25.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
  """)


# In[ ]:


#Energy Minimization Stage 2
with open("min1.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Energy Minimization, Stage2
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=1, !Energy Minimization
  maxcyc=5000, !Total Minimization Cycles to be run. Steepest Decent First, then Conjugate Gradient Method if ncyc < maxcyc
  ncyc=1000, !Number of Steepest Decent Minimization Steps to run before switching to Conjugate Gradient
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=1, !Complete Interactions are Calculated
  ntc=1, !SHAKE is NOT performed, DEFAULT
  ntpr=10, !Every 10 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=10, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=100, !Every 100 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@N,CA,C,O", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=10.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")


# In[ ]:


#Energy Minimization Stage 3
with open("min2.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Energy Minimization, Stage3
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=1, !Energy Minimization
  maxcyc=5000, !Total Minimization Cycles to be run. Steepest Decent First, then Conjugate Gradient Method if ncyc < maxcyc
  ncyc=1000, !Number of Steepest Decent Minimization Steps to run before switching to Conjugate Gradient
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=1, !Complete Interactions are Calculated
  ntc=1, !SHAKE is NOT performed, DEFAULT
  ntpr=10, !Every 10 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=10, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=100, !Every 100 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=5.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")

#Heating Stage

with open("heat.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Heating,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Do NOT Restart the Simulation; instead, run as a NEW Simulation
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=1, !Coordinates, but no Velocities, will be read. Formatted (ASCII) coordinate file is expected
  nstlim=5000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=10000, !Every 100 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=10000, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 1000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  get_ipython().system('temp0=310.00000, !Reference temperature at which the system is to be kept')
  tempi=0.00000, !Initial Temperature
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=5.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")


init_step = 0
total_step = 500000
step_inc  = 8928
temp_inc  = 5
init_temp = 50.0
norm_temp = 310.0
max_temp  = 320.0

heat_in = "heat.in"
topfile  = "temp_inc.in"
    
TempFile = open(topfile, 'w')
while (init_temp < 320.0):
    TempFile.writelines("&wt type='TEMP0', istep1=" + str(init_step) + ", istep2=" + str(init_step+step_inc) + ", value1=" + str(init_temp) + ", value2=" + str(init_temp+5) +", /" + "\n")
    init_step = init_step + step_inc
    init_temp = init_temp + 5

while (max_temp > 310.0):
    if max_temp == 310.0:
        TempFile.writelines("&wt type='TEMP0', istep1=" + str(init_step) + ", istep2=" + str(total_step) + ", value1=" + str(norm_temp) + ", value2=" + str(norm_temp) +", /" + "\n")
    else:
        TempFile.writelines("&wt type='TEMP0', istep1=" + str(init_step) + ", istep2=" + str(init_step+step_inc) + ", value1=" + str(max_temp) + ", value2=" + str(max_temp-5) +", /" + "\n")
        init_step = init_step + step_inc
        max_temp = max_temp - 5
        
TempFile.writelines("&wt type='END' /")
TempFile.close()

#merge the heating input file and the step by step temp increment
merge_file = "cat " + heat_in + " " + topfile + " > heating.in"
os.system(merge_file)


# NVT Equilibration

with open("equi_NVT.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: NVT Equilibration,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Do NOT Restart the Simulation; instead, run as a NEW Simulation
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=1, !Coordinates, but no Velocities, will be read. Formatted (ASCII) coordinate file is expected
  nstlim=5000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=0, !No pressure scaling (Default)
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=10000, !Every 100 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=50000, !Every 100 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 1000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  temp0=310.00000, !Reference temperature at which the system is to be kept
  tempi=310.00000, !Initial Temperature
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=2.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
 """)

# NPT Equilibration

with open("equi_NPT.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: NPT Equilibration,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Do NOT Restart the Simulation; instead, run as a NEW Simulation
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=5, !Coordinates and Velocities will be read; a formatted (ASCII) coordinate file is expected.
  nstlim=10000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=1, !MD with isotropic position scaling
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=1000, !Every 1000 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=50000, !Every 10000 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 10000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  temp0=310.00000, !Reference temperature at which the system is to be kept
  tempi=310.00000, !Initial Temperature
  pres0=1.01300, !Reference Pressure (in units of bars, where 1 bar = 0.987 atm) at which the system is maintained
  ntr=1, !Turn ON (Cartesian) Restraints
  restraintmask="@CA", !Atoms to be Restrained are specified by a restraintmask
  restraint_wt=2.00, !Force Constant for Restraint, kcal/(mol * A^2)
/
""")

# Production Run
with open("md_NPT.in", "w") as file:
    file.write(
"""Type of Simulation Being Done: Production Run,
 &cntrl
  ntxo=2, IOUTFM=1, !NetCDF Binary Format.
  imin=0, !MD Simulation
  irest=0, !Start new Simulation;
  ig=-1, !Pseudo-random number seed is changed with every run.
  ntx=5, !Coordinates and Velocities will be read; a formatted (ASCII) coordinate file is expected.
  nstlim=10000000, !Number of MD-steps to be performed. Default 1.
  dt=0.00200, !The time step (psec). Recommended MAXIMUM is .002 if SHAKE is used, or .001 if SHAKE is NOT used
  cut=10.00000, !Cut Off Distance for Non-Bounded Interactions
  igb=0, !No Generalized Born
  ntp=1, !MD with isotropic position scaling
  ntf=2, !Bond Interactions involving H-atoms omitted
  ntc=2, !Bonds involving Hydrogen are Constrained
  ntpr=1000, !Every 1000 steps, energy information will be printed in human-readable form to files "mdout" and "mdinfo"
  ntwx=50000, !Every 10000 steps, the coordinates will be written to the mdcrd file
  ntwr=10000, !Every 10000 steps during dynamics, the restart file will be written, ensuring that recovery from a crash will not be so painful. #If ntwr < 0, a unique copy of the file, "restrt_<nstep>", is written every abs(ntwr) steps
  ntt=3, !Use Langevin Dynamics with the Collision Frequency GAMA given by gamma_ln,
  gamma_ln=2.00000, !Collision Frequency, ps ^ (-1)
  temp0=310.00000, !Reference temperature at which the system is to be kept
  tempi=310.00000, !Initial Temperature
  pres0=1.01300, !Reference Pressure (in units of bars, where 1 bar = 0.987 atm) at which the system is maintained
/
 """)


# In[46]:


#Submission File
with open("job_submit.sh", "w") as file:
    file.write(
"""#!/bin/bash
#SBATCH --job-name=
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=user@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --gres=gpu:1g.5gb:1
#SBATCH --output=out.log
#SBATCH --partition=COMPUTE1Q
#SBATCH --account=yanglab

topfile = "protein_solv_ions.parm7"
coordfile = "protein_solv_ions.crd"

#Energy Minimization
singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda -O -i min0.in -p $topfile -c $coordfile -o min0.out -r min0.rst -ref $coordfile -x min0.nc -inf min0.info

singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda -O -i min1.in -p $topfile -c min0.rst -o min1.out -r min1.rst -ref min0.rst -x min1.nc -inf min1.info

singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda -O -i min2.in -p $topfile -c min1.rst -o min2.out -r min2.rst -ref min1.rst -x min2.nc -inf min2.info

#NVT Heating
singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda_SPFP -O -i heating.in -p $topfile -c min2.rst -o heat.out -r heat.rst -ref min2.rst -x heat.nc -inf heat.info

#NVT Equilibration
singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda_SPFP -O -i equi_NVT.in -p $topfile -c heat.rst -o equi_NVT.out -r equi_NVT.rst -ref heat.rst -x equi_NVT.nc -inf equi_NVT.info

#NPT Equilibration
singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda_SPFP -O -i equi_NPT.in -p $topfile -c equi_NVT.rst -o equi_NPT.out -r equi_NPT.rst -ref equi_NVT.rst -x equi_NPT.nc -inf equi_NPT.info

#Production Run
singularity exec --nv --bind /raid:/raid /raid/images/amber20.sif pmemd.cuda_SPFP -O -i md_NPT.in -p $topfile -c equi_NPT.rst -o md_NPT.out -r md_NPT.rst -ref equi_NPT.rst -x md_NPT.nc -inf md_NPT.info
""")






