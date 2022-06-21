import numpy as np
from chemdiff import *
from chemdiff.parallel import *
import os
import subprocess
from mpi4py import MPI

###### OPTIONAL #########
# if you don't have this module just comment out the `start` and `end` lines at beginning and end
from timeit import default_timer

start = default_timer()

################ INITIALIZE MPI ################
comm,nproc,rank = mpi_initialize()

############ SETUP AND CHEM PARAMS #############
f_pebout = 'pebble_composition.out'
if not os.path.isfile(f_pebout):
        with open(f_pebout,'w') as f:
                f.write('time    species    pebble_col_abundance\n')

grow_pebbles = True

# read in chemical file here
with open('./r00/z00/input.ini','r') as f:
        f.readline()
        f.readline()
        chm = f.readline().split()[-1]


############ INITIALIZE THE COLUMN #############
ti = 8.e5
tf = 1.e6 # yrs

# get midplane temperature
with open('./r00/z00/source.mdl','r') as f:
        f.readline()
        params = list(map(float,f.readline().split()))
tmid = params[3]
r = params[-2]*au

# model paramters
# should be same as priming for consistency!
alpha = 1e-3
nzs = 50
dt = 100 # yrs
chemtime = 500 # yrs
opacity = 1e3

nchems = int((tf-ti)/chemtime)
ndiffs = int(chemtime/dt)

# times to output results
touts = []
touts += [(i+1.)*pow(10,3) for i in range(9)]
touts += [(i+1.)*pow(10,4) for i in range(9)]
touts += [(i+2)/2*pow(10,5) for i in range(18)]
touts += [1.e6]
touts = np.array(touts)
touts = touts[touts<=tf]
nts = len(touts)

# create the column
col = Column(r,tmid,alpha,nzs)
col.set_diff_params(dt)
h = col.h

############# INITIALIZE THE CELLS #############

# read in the physical params
# these are constant for each cell (with the exception of grain_gas_mass_ratio,
# which will be read in later)
phys_params = {}
solver_params = {}
with open('./r00/z00/input.ini','r') as f:
        switch = None
        for line in f:
                line.rstrip('\n')
                if line[0] == '[' and line[-2] == ']':
                        switch = line[1:-2]
                elif switch == 'phys':
                        label,_,value = line.split()
                        phys_params[label] = float(value)
                elif switch == 'solver':
                        label,_,value = line.split()
                        solver_params[label] = float(value)

if rank == 0:
        if not os.path.exists('r00'):
                os.system('mkdir r00')

        # read info for each cell
        NH = 0
        for j in reversed(range(nzs)):
                dirr=f'r00/z{j:0>2}'
                if not os.path.exists(dirr):
                        os.system('mkdir '+dirr)
                z = col.dz*(j+0.5) # cm

                # read in cell info from source.mdl
                with open(f'./r00/z{j:0>2}/source_t{int(ti):0>9}.mdl','r') as f:
                        f.readline()
                        cell_params = f.readline().split()
                av,nh,tgas,tdust,NCO,NH2,NHD,xray = list(map(float,cell_params[1:-2]))

                # read in abundances from last astrochem output
                cell_abuns = get_final_abuns(f'./r00/z{j:0>2}/astrochem_output_t{int(ti):0>9}.h5')
                

                ### CREATE THE CELL
                chi = phys_params['chi']
                cosmic = phys_params['cosmic']
                grain_size = phys_params['grain_size']
                dg = grain_abun2dg(cell_abuns['grain'])
                print('dg = ',dg)
                NH += nh*col.dz
                col.cells[j] = Cell(r,z,chi=chi,cosmic=cosmic,grain_size=grain_size,
                                    dust_gas_ratio=dg,av=av,rho=mbar*nh/2,Tgas=tgas,
                                    Tdust=tdust,xray=xray,NCO=NCO,NH2=NH2,NHD=NHD,NH=NH,
                                    abundances = dict(cell_abuns))
if rank == 0:
        for i in range(1,nproc):
                comm.send(col,dest=i,tag=22*i)
else:
        comm.recv(source=0,tag=22*rank)
print('initialized column and cells -- rank :',rank)


################ MAIN LOOP #####################
peb_comp = {}
cwd = os.getcwd()
for t in range(nchems):
        time = ti+(t+1)*chemtime

        # do chemistry
        do_parallel_chemistry(col,comm,nproc,rank,
                              time,touts,chemtime=chemtime,network=chm)

        # rank 0 will do the pebble growth, diffusion, and update the cells
        goon = False
        if rank == 0:
                for diff_loop in range(ndiffs):
                        # change grain abundances
                        peb_comp = grow_grains(col,peb_comp,time,grow_pebbles)
                        do_diffusion(col)
                update_cells(col,opacity)

                if float(time) in touts:
                        print(time)
                        col_abunds = col.get_abundance_array()
                        with open(f_pebout,'a') as f:
                                for spec in peb_comp:
                                        peb_comp_norm = peb_comp[spec]/col.cells[0].NH
                                        f.write(f'{time}    {spec}    {peb_comp_norm:.10e}\n')

                goon = True
                # send the 'go on' signal and the updated column to other ranks
                for i in range(1,nproc):
                        comm.send(True,dest=i,tag=15*i)
                        comm.send(col,dest=i,tag=155*i)
        else:
                goon = comm.recv(source=0,tag=15*rank)
                col = comm.recv(source=0,tag=155*rank)

        while not goon:
                sleep(1)
                print('WAITING :: RANK ',rank)


##################### TIMER ####################
if rank == 0:
        end = default_timer()
        time_to_run = end-start
        hrs = time_to_run//3600
        mins = (time_to_run-(hrs*3600))//60
        secs = time_to_run-(hrs*3600+mins*60)
        with open('time_to_run.out','a') as f:
                f.write(f'time to run : {time_to_run:.1f} sec\n')
                f.write(f'            = {int(hrs)}h {int(mins)}m {secs:.1f}s\n')

print('done!')
