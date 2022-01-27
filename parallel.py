import numpy as np
import subprocess
import os
from .constants import *
from .wrapper import *
from mpi4py import MPI
from time import sleep
from timeit import default_timer

def mpi_initialize():
        '''
        Initialized MPI comm for parallel computing
        '''
        comm = MPI.COMM_WORLD
        size = comm.size
        rank = comm.rank
        return comm,size,rank

def do_parallel_chemistry(col,comm,size,rank,time,touts,chemtime=500,
        network='network.chm',abs_err=1.e-20,rel_err=1e-10,verbose=False):
        '''
        Runs Astrochem in parallel

        PARAMETERS
        ----------
        col : Column object to run astrochem on
        comm : MPI.COMM_WORLD for mpi4py
        size : int, number of processors to use
        rank : int, the rank of the cpu working on this
        time : float, the current time, check if output should be saved
        touts : list, times for which to save outputs
        astrochem_params : args to passed to run astrochem (see wrapper.write_chem_inputs)
                chemtime, network, abs_err, rel_err
        '''

        ### divide the cells among the processors
        ### THIS IS VERY INNEFICIENT IF size IS NOT A FACTOR OF col.ncells !!!
        ### THIS SHOULD BE FIXED
        my_cells = None
        nzs = col.ncells
        js_per_rank = nzs//size
        remain = nzs - (size*js_per_rank)
        if rank == 0:
                for i in range(1,size):
                        comm.send(col.cells[i*js_per_rank+remain:(i+1)*js_per_rank+remain],dest=i,tag=i*11)
                my_cells = col.cells[0:js_per_rank+remain]
        else:
                my_cells = comm.recv(source=0,tag=rank*11)

        ### Call astrochem to do the chemistry
        cwd = os.getcwd()
        wait = True
        if rank == 0:
                all_done = np.zeros(size,dtype='bool')
        for my_j,cell in enumerate(my_cells):
                # get the absolute j for each cells
                if rank == 0:
                        j = my_j
                else:
                        j = rank*js_per_rank+remain+my_j
                dirr = f'{cwd}/r00/z{j:0>2}'
                if verbose:
                        cell_start = default_timer()
                        print('beginning cell ',j)
                cell.write_chem_inputs(chemtime,abs_err=abs_err,rel_err=rel_err,
                        f_net=network,f_input=dirr+'/input.ini',f_source=dirr+'/source.mdl')
                subprocess.run(['astrochem','-q','input.ini'],cwd=dirr)
                if float(time) in touts:
                        subprocess.run(['cp','astrochem_output.h5',f'astrochem_output_t{int(time):0>9}.h5'],cwd=dirr)
                        subprocess.run(['cp','source.mdl',f'source_t{int(time):0>9}.mdl'],cwd=dirr)

                # update cell abundances
                out_times,d = get_abundict(f'{dirr}/astrochem_output.h5','all')
                for spec in d:
                        cell.abundances[spec] = d[spec][-1]
                if verbose:
                        cell_end = default_timer()
                        print(f'  done with cell {j} in {cell_end-cell_start:.1f} seconds')
                        print()
        ### update all the cells on rank 0
        if rank != 0:
                comm.send(my_cells,dest=0,tag=rank*99)
        else:
                for i in range(1,size):
                        col.cells[i*js_per_rank+remain:(i+1)*js_per_rank+remain] = comm.recv(source=i,tag=i*99)

        '''
        once all the cells have finished chemistry send a signal to rank 0
        once rank 0 has recieved the signals from everyone, send the all clear to go on
        I'm sure there is a better way to do this/this may not actually be necessary
        But it doesn't take very long to go through this and it's a good safeguard
        '''
        # first, all ranks except 0 send that they are done
        if rank != 0:
                comm.send(True, dest=0,tag=rank*12)
        else:
                all_done[0] = True
                for i in range(1,size):
                        all_done[i] = comm.recv(source=i,tag=i*12)
        # now rank zero says "everyone is done!" and sends to go ahead
        if rank == 0:
                if all(all_done):
                        wait = False
                        for i in range(1,size):
                                comm.send(False,dest=i,tag=i*13)
        else:
                wait = comm.recv(source=0,tag=13*rank)
        # now everyone should be together having both sent and recieved the ok
        while wait:
                print('WAITING :: RANK',rank)
                sleep(1)


def grow_grains(col,peb_comp,time,grow_pebbles=True,timescale_factor=1.,grow_height=1.,prime_time=0.):
        '''
        Grows pebbles near the midplane with timescale from Birnsteil 2012

        PARAMETERS
        ----------
        col : Column object
        peb_comp : dict, current pebble composition dictionary
        time : float, current time (yrs)
        grow_pebbles : bool, should pebbles be grown or not?
        timescale_factor : float, multiplicative factor for growth timescale 
                                                tau = timescale_factor/tau_0
        grow_height : float, height in scale heights below which pebbles should grow
        OUTPUT
        ------
        dict, update peb_comp dictionary
        '''
        if not prime_time:
                prime_time = 0.
        nzs = col.ncells
        for j in range(nzs):
                cell = col.cells[j]
                t_grow = timescale_factor/(cell.dust_gas_ratio*col.omega)
                if cell.z/col.h <= grow_height and time >=prime_time and grow_pebbles:
                        deps = -col.dt*sec/t_grow
                else:
                        deps = 0
                for spec in cell.abundances:
                        if spec[:5] == 'grain':
                                d_ice = deps*cell.abundances[spec] # X*dt/tau
                                cell.abundances[spec] += d_ice
                                if spec not in peb_comp:
                                        peb_comp[spec] = 0
                                peb_comp[spec] -= d_ice*cell.nh*col.dz # col averaged abundance
        return peb_comp

def do_diffusion(col):
        '''
        Do the diffusion calculation

        PARAMETERS
        ----------
        col : Column object
        '''
        nzs = col.ncells
        col_abunds = col.get_abundance_array()
        newarray = {}
        for spec in col_abunds:
                newarray[spec] = np.zeros(nzs)
                for j in range(nzs):
                        rhos_j = col.cells[j].rho
                        sp = col_abunds[spec]
                        fp = 0
                        if j < nzs-1:
                                rhos_j1 = col.cells[j+1].rho
                                fp = 0.5*(rhos_j+rhos_j1)*col.beta*(sp[j+1]-sp[j])/col.dz
                        elif j == nzs-1:
                                fp = 0.
                        if j==0:
                                fm = -fp
                        newarray[spec][j] = sp[j]+(fp-fm)/col.dz/rhos_j
                        fm = fp

        # save the new values
        for spec in newarray:
                for j in range(nzs):
                        col.cells[j].dust_gas_ratio = grain_abun2dg(newarray['grain'][j])
                        col.cells[j].abundances[spec] = newarray[spec][j]

def update_cells(col,opacity=1e5):
        '''
        Update the cells after diffusion

        PARAMETERS
        ----------
        col : Column Object
        opacity : total dust opacity
        '''
        nzs = col.ncells

        NCO = 0.
        NH2 = 0.
        NH = 0.
        tau = 0.
        for j in reversed(range(nzs)):
                cell = col.cells[j]
                nh = cell.nh
                nco = 0
                if 'CO' in list(cell.abundances.keys()):
                        nco = cell.abundances['CO']*nh
                nh2 = 0
                if 'H2' in list(cell.abundances.keys()):
                        nh2 = cell.abundances['H2']*nh
                NCO += nco*col.dz
                NH2 += nh2*col.dz
                NH += nh*col.dz
                tau += cell.rho*opacity*col.dz*cell.dust_gas_ratio
                cell.NCO = NCO
                cell.NH2 = NH2
                cell.NH = NH
                cell.av = tau/3.02
