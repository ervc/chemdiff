from .constants import *
import numpy as np
import h5py as h5
import sys

def grain_abun2dg(grain_abun,grain_density=1.81,grain_size=1e-5):
    dg = grain_abun*grain_density*(4./3.)*np.pi*pow(grain_size,3.)/mh
    return dg

def dg2grain_abun(dg,grain_density=1.81,grain_size=1e-5):
	grain_abun = dg*mh/(grain_density*(4./3.)*np.pi*pow(grain_size,3.))
	return grain_abun

def readfile(filename='astrochem_output.h5'):
	f = h5.File(filename,'r')
	spec = f['Species']
	abun = f['Abundances']
	time = f['TimeSteps']
	return f,spec,abun,time

def get_abundict(filename='astrochem_output.h5',specs = 'all'): 
    f,spec,abun,time = readfile(filename)

    d = {}
    time = np.array(time)
    for i in range(spec.size):
        specie = spec[i].decode('utf-8')
        if specie in specs or specs == 'all':
            abund = abun[0,:,i]
            d[specie] = abund
    f.close()
    return time,d

def get_final_abuns(f_name,specs='all'):
	final_abuns = {}
	times,d = get_abundict(f_name,specs=specs)
	for spec in d:
		final_abuns[spec]=d[spec][-1]
	return final_abuns

class Column:

	def __init__(self,r,tmid,alpha=1e-3,ncells=50):
		self.r = r # cm
		self.ncells = ncells
		self.tmid = tmid # K
		self.omega = np.sqrt(g*msun/r)/r #1/s
		self.cs = np.sqrt(bk*tmid/mbar) #cm/s
		self.h = self.cs/self.omega # cm

		self.alpha = alpha
		self.dz = 5*self.h/ncells # cm

		self.cells = np.empty(ncells,dtype='object')

	def set_diff_params(self,dt):
		self.dt = dt # yr

		self.diff = self.alpha*self.cs*self.h
		self.beta = 0.5*self.diff*self.dt*sec

	def get_abundance_array(self):
		all_abunds = {}
		for i in range(self.ncells):
			cell = self.cells[i]
			for spec in cell.abundances:
				if spec not in all_abunds:
					all_abunds[spec] = np.zeros(self.ncells)
				all_abunds[spec][i] = cell.abundances[spec]
		self.all_abunds = all_abunds
		return all_abunds

class Cell(object):
	'''
	A cell has:

	a physical location
		r,z

	physical parameters for the solver
		chi, cosmic, grain size, dust/gas ratio

	source model parameters
		given
			Av, nh, Tgas, Tdust, xrays
		calculated
			NCO, NH2, (density)

	abundances
		dict of abundances

	diffusion parameters
		given 
			omega, 
			alpha, dz, dt (solver params?)
		calculated
			cs, h, D, beta
	'''


	def __init__(self,r,z,
			chi=1,cosmic=1.3e-17,grain_size=0.1,dust_gas_ratio=0.01,
			av=1,rho=1e10,Tgas=50,Tdust=50,xray=0,NCO=1.,NH2=1.,NH=1.,
			abundances = {}):
		self.r = r # cm
		self.z = z # cm

		# phys
		self.chi = chi
		self.cosmic = cosmic # s-1
		self.grain_size = grain_size # micron
		self.dust_gas_ratio = dust_gas_ratio

		self.av = av # mag
		self.rho = rho # g.cm-3
		self.nh = 2*rho/mbar # cm-3
		self.Tgas = Tgas # K
		self.Tdust = Tdust # K
		self.xray = xray # s-1
		self.NCO = NCO # cm-2
		self.NH2 = NH2 # cm-2
		self.NH = NH # cm-2

		self.abundances = dict(abundances)


	def write_chem_inputs(self,tf,abs_err,rel_err,abun_out='all',f_net='network.chm',f_input='input.ini',f_source='source.mdl'):
		with open(f_input,'w') as f:
			f.write('[files]\n')
			f.write(f'source = {f_source}\n'+
				f'chem = {f_net}\n')

			f.write('[phys]\n')
			f.write(f'chi = {self.chi:.3e}\n'+
				f'cosmic = {self.cosmic:.2e}\n'+
				f'grain_size = {self.grain_size:.2e}\n'+
				f'grain_gas_mass_ratio = {self.dust_gas_ratio:.2e}\n')

			f.write('[solver]\n')
			f.write('ti = 1.00e-06\n'+
				f'tf = {tf:.2e}\n'+
				f'abs_err = {abs_err:.1e}\n'+
				f'rel_err = {rel_err:.1e}\n')

			f.write('[abundances]\n')
			for spec in self.abundances:
				try:
					if self.abundances[spec] != 0:
						f.write(f'{spec} = {self.abundances[spec]:.15e}\n')
				except:
					print('Problem with Abundances:\n',self.abundances)
					sys.exit()


			f.write('[output]\n')
			f.write(f'abundances = {abun_out}\n'+
				'time_steps = 128\n'+
				'trace_routes = 0')

		# print(f'Input written to {f_input}')

		with open(f_source,'w') as f:
			f.write('# self Av[mag] n(H)[cm-3] Tgas[K] Tdust[K] NCO[cm-2] NH2[cm-2] xray-ion[s-1] R[au] Z[au]\n')
			f.write(f'0    {self.av:.3e}    {self.nh:.3e}    {self.Tgas:.3e}    {self.Tdust:.3e}    {self.NCO:.3e}    {self.NH2:.3e}    {self.xray:.3e}    {self.r/au:.2f}    {self.z/au:.2f}')

		# print(f'Source written to {f_source}')

	def update_abundances(self,new_dict):
		self.abundances=dict(new_dict)







		