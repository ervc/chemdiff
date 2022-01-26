import chemdiff.constants as const

# default model paramters
model_inputs = {
	'chmfile' : 'umist12.chm',
	'pebfile' : 'pebble_composition.out',
	'ti' : 0.,
	'tf' : 1e6,
	'touts' : [1.0e3,2.0e3,3.0e3,4.0e3,5.0e3,6.0e3,7.0e3,8.0e3,9.0e3,
			   1.0e4,2.0e4,3.0e4,4.0e4,5.0e4,6.0e4,7.0e4,8.0e4,9.0e4,
			   1.0e5,1.5e5,2.0e5,2.5e5,3.0e5,3.5e5,4.0e5,4.5e5,5.0e5,
			   5.5e5,6.0e5,6.5e5,7.0e5,7.5e5,8.0e5,8.5e5,9.0e5,9.5e5,
			   1e6],
	'diff_dt' : 100,
	'chem_dt' : 500,
}

# default physical paramters
phys_inputs = {
	'r' : 30*const.au,
	'alpha' : 1.e-3,
	'chi' : 50,
	'cosmic' : 1.3e-17,
	'grain_size' : 0.1,
	'dg0' : 0.01,
	'opacity' : 1000,
}

# default initial chemical abundances
abundances = {
	'H2' : 0.5,
	'He' : 9.75e-2,
	'NH3': 1.45e-6,
	'H2O': 1.18e-4,
	'CO' : 6.00e-5,
	'N2' : 2.00e-5,
	'CH4': 2.00e-6,
	'CH3OH' : 1.00e-6,
	'H2S' : 1.91e-8,
	'CO2' : 5.00e-5,
	'HCN' : 3.50e-7,
	'grain' : 2.2e-12
}