import numpy as np
import chemdiff as cd
import argparse
import os
import subprocess

def directory_exists(dir_name):
    ''' Checks if a directory exists '''
    return os.path.isdir(dir_name)

def make_directory(dir_name):
    ''' Makes a directory '''
    subprocess.run(['mkdir',dir_name])
    return

def format_directory(dir_name):
    ''' removes trailing / from directory name for consistency '''
    if dir_name[-1] == '/':
        dir_name = dir_name[:-1]
    return dir_name

def check_model_directory(dir_name):
    ''' Checks if a directory exists. If it does, check if directory r00/ exists inside directory '''
    if not directory_exists(model_dir):
        raise FileNotFoundError(f'Directory does not exist: {model_dir}')
    if not directory_exists(model_dir+'/r00'):
        raise FileNotFoundError(f'No model results found in directory: {model_dir}')
    return 0

def get_specs_array(dir_name):
    ''' Get an array of specs output by model '''
    cti = f't{int(1e3):0>9}'
    if not os.path.isfile(dir_name+f'/r00/z00/astrochem_output_{cti}.h5'):
        cti = f't{int(1e3):0>7}'
    times,abundict = cd.get_abundict(dir_name+f'/r00/z00/astrochem_output_{cti}.h5')
    specs = np.empty(len(list(abundict.keys())),dtype='U25')
    for k,spec in enumerate(list(abundict.keys())):
        specs[k] = spec
    return specs

def get_time_outs(tf):
    ''' Create array of times with outputs '''
    touts = []
    exponent = 3
    final_exp = np.log10(tf)//1
    while exponent < final_exp:
        if exponent < 5:
            touts += [(i+1.)*10**exponent for i in range(9)]
        else:
            touts += [(i+2.)/2*10**exponent for i in range(18)]
        exponent += 1
    tout = 1*10**final_exp
    touts.append(tout)
    i = 1
    while tout < tf:
        tout = (i+2.)/2*10**final_exp
        touts.append(tout)
        i+=1
    touts = np.array(touts)
    return touts

def get_nh_and_av(source_file):
    ''' Read the gas density (nh) and visual extinction (av) from source.mdl file '''
    with open(source_file,'r') as f:
        f.readline()  # skip header
        params = f.readline().split()
    nh = float(params[2])
    av = float(params[1])
    return nh,av
        

def main(model_dir,output_name,tf,quiet):
    ''' Gather abundances, times, species, abundances, gas densities, and visual extinctions
    into an output file for faster reading later on.
    '''
    check_model_directory(model_dir)

    # get the time and height independent arrays
    nzs = 50
    specs = get_specs_array(model_dir)
    nspec = len(specs)
    touts = get_time_outs(tf)
    nts = len(touts)
    heights = np.linspace(0.05,4.95,nzs)

    # loop through for height and time dependet arrays
    gas_dens = np.zeros(nzs)
    extinctions = np.zeros((nts,nzs))
    abundances = np.zeros((nts,nzs,nspec))
    for j in range(nzs):
        dirr = model_dir+f'/r00/z{j:0>2}'
        if not quiet:
            print(dirr)
        for t,time in enumerate(touts):
            if not quiet:
                print(f'  {time}')
            ct = f't{int(time):0>9}'
            if not os.path.isfile(dirr+f'/source_{ct}.mdl'):
                ct = f't{int(time):0>7}'
            nh,av = get_nh_and_av(dirr+f'/source_{ct}.mdl')
            extinctions[t,j] = av
            if t==0:
                gas_dens[j] = nh
            final_dict = cd.get_final_abuns(dirr+f'/astrochem_output_{ct}.h5')
            for k,spec in enumerate(list(final_dict.keys())):
                abundances[t,j,k] = final_dict[spec]

    # save
    np.savez(model_dir+'/'+output_name,times=touts,species=specs,abundances=abundances,
             heights=heights,gas_dens=gas_dens,extinctions=extinctions)
    if not quiet:
        print(f'saved output to: {model_dir}/{output_name}')

    return
    

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='gather abundances into npz file. Abundances have shape (nts,nzs,nspec). Also saves specs(nspec), touts(nts), heights(nzs), gas_dens(nzs), avs(nts,nzs)')
    parser.add_argument('directory',help='Path to model directory to get abundances from.',type=str)
    parser.add_argument('-o','--output',help='Name of npz file to create (saved to directory). Default = {directory}/output_abuns.npz',
                        required=False,default='output_abuns.npz',type=str)
    parser.add_argument('-tf',metavar='tf',help='Final time to get abundances for. Default = 1e6 yr',
                        required=False,default=1e6,type=float)
    parser.add_argument('-q','--quiet',help='Silence output of loop through times and heights.',
                        default=False)

    args = parser.parse_args()
    model_dir = format_directory(args.directory)
    output_name = args.output
    quiet = args.quiet
    tf = args.tf
    main(model_dir,output_name,tf,quiet)
