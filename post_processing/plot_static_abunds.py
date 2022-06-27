import numpy as np
import matplotlib.pyplot as plt
import chemdiff as cd
import argparse

NZS = 50
DZ = 5/50
ZS = [j*DZ + 0.5*DZ for j in range(NZS)]

R = 30

def get_deuterated_form(spec):
	newspec = ''
	if 'H' not in spec:
		return spec
	hindex = spec.index('H')
	if hindex+1 >= len(spec):
		newspec = spec[:hindex]+'D'
		return newspec
	if spec[hindex+1].isnumeric():
		if spec[hindex+1] == '2':
			newspec = spec[:hindex]+'HD'+spec[hindex+2:]
		else:
			newspec = spec[:hindex]+f'H{int(spec[hindex+1])-1}D'+spec[hindex+2:]
	else:
		newspec = spec[:hindex]+'D'+spec[hindex+1]
	return newspec

def count_hydrogen(spec):
	if 'H' not in spec:
		return 0
	hindex = spec.index('H')
	if hindex+1 >= len(spec):
		return 1
	if spec[hindex+1].isnumeric():
		return int(spec[hindex+1])
	else:
		return 1

def count_deuterium(spec):
	if 'D' not in spec:
		return 0
	dindex = spec.index('D')
	if dindex+1 >= len(spec):
		return 1
	if spec[dindex+1].isnumeric():
		return int(spec[dindex+1])
	else:
		return 1

def format_directory(dir_name):
    ''' removes trailing / from directory name for consistency '''
    if dir_name[-1] == '/':
        dir_name = dir_name[:-1]
    return dir_name

def main(args):
	dirr = format_directory(args.directory)
	specs = args.species

	# x and y dictionaries for each species
	if args.xaxis=='abund':
		x = {}
		y = {}
	elif args.xaxis=='HD_to_H2':
		x = []
		y = []
	elif args.xaxis=='D_to_H':
		specs = ['HD','D','tot']
		x={}
		y={}

	ct = f't{int(args.time):0>9}'
	for j in range(NZS):
		if not args.quiet:
			print(j)
		times,abundict = cd.get_abundict(dirr+f'/r00/z{j:0>2}/astrochem_output_{ct}.h5')

		Y = ZS[j]

		if args.xaxis=='abund':
			for spec in args.species:
				X = abundict[spec][-1]
				if spec in x:
					x[spec].append(X)
					y[spec].append(Y)
				else:
					x[spec] = [X]
					y[spec] = [Y]
				dspec = get_deuterated_form(spec)
				X = abundict[dspec][-1]
				if dspec in x:
					x[dspec].append(X)
					y[dspec].append(Y)
				else:
					x[dspec] = [X]
					y[dspec] = [Y]


		elif args.xaxis=='HD_to_H2':
			X = abundict['HD'][-1]/abundict['H2'][-1]
			x.append(X)
			y.append(Y)


		elif args.xaxis=='D_to_H':
			hd = abundict['HD'][-1]
			hh = abundict['H2'][-1]
			d = abundict['D'][-1]
			h = abundict['H'][-1]

			totd = d+hd # = 2e-5
			toth = hd+2*hh+h # = 1

			if j==0:
				# x['HD'] = [hd/(hd+2*hh)]
				# x['D'] = [d/h]
				x['tot'] = [totd/toth]
				y['tot'] = [ZS[0]]
				for spec in args.species:
					hspec = spec
					dspec = get_deuterated_form(hspec)
					alld = count_deuterium(dspec)*abundict[dspec][-1] + count_deuterium(hspec)*abundict[hspec][-1]
					allh = count_hydrogen(dspec)*abundict[dspec][-1] + count_hydrogen(hspec)*abundict[hspec][-1]
					X = alld/allh
					x[spec] = [X]
					y[spec] = [ZS[0]]
			else:
				# x['HD'].append(hd/(hd+2*hh))
				# x['D'].append(d/h)
				x['tot'].append(totd/toth)
				y['tot'].append(ZS[j])
				for spec in args.species:
					hspec = spec
					dspec = get_deuterated_form(hspec)
					alld = count_deuterium(dspec)*abundict[dspec][-1] + count_deuterium(hspec)*abundict[hspec][-1]
					allh = count_hydrogen(dspec)*abundict[dspec][-1] + count_hydrogen(hspec)*abundict[hspec][-1]
					X = alld/allh
					x[spec].append(X)
					y[spec].append(ZS[j])





	fig,ax = plt.subplots()
	if args.xaxis=='abund':
		print(x)
		print(y)
		for spec in list(x.keys()):
			ax.plot(x[spec],y[spec],label=spec)
	elif args.xaxis=='HD_to_H2':
		ax.plot(x,y)
	elif args.xaxis=='D_to_H':
		for spec in args.species:
			label=f'[D/H]_{spec}'
			ax.plot(x[spec],y[spec],label=label)
		ax.plot(x['tot'],y['tot'],label='Total')
	title = f'R = {R} au\nt = {args.time:.1e} yr'
	ax.set(xscale='log',xlabel=args.xaxis,ylabel=args.yaxis,title=title,ylim=(0,5))
	ax.legend()
	plt.savefig(f'{args.directory}/shieldtest_{args.xaxis}.png',bbox_inches='tight')
	plt.show()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="plot abundances at a given time for the column")
	parser.add_argument('directory',help='path to directory',type=str)
	parser.add_argument('-t','--time',help='time to plot',type=float,default=1e6)
	parser.add_argument('-specs','--species',help='species to plot',type=str,default=['H2','CO','H2O'],nargs='+')
	parser.add_argument('-x','--xaxis',help='what to plot on xaxis',default='abund',type=str)
	parser.add_argument('-y','--yaxis',help='what to plot on yaxis',default='scale_height',type=str)
	parser.add_argument('-q','--quiet',help='suppress outputs',action='store_true')

	main(parser.parse_args())