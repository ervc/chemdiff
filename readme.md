# Chemdiff

## Overview

Chemdiff is the python based user interface for the CANDY package. Chemdiff repeatedly calls Astrochem, while allowing for vertical diffusion and dynamic pebble growth and ice sequestration in a protoplanetary disk. For a details on the physical structure of the disk (i.e. temperature/density structure) see Van Clepper et al. 2022.

At its core, chemdiff has three functions it calls repeatedly:
1) Update Chemistry
2) Diffuse Material
3) Grow Pebbles

These rely on two objects, located in the `wrapper.py` module: the Cell and a Column. Each Cell represents a small volume of constant density, temperature, composition, etc. located at some height above the midplane. A Column has a collection of Cells (default is 50 cells, ranging from midplane to 5 scale heights). Within a Column, material can be diffused between Cells. The `parallel.py` module contains functions to do these three steps, while the `run_parallel_growth.py` file located in the `python_examples/` directory contains a loop for these three functions.

## Running chemdiff

To run, make sure both chemdiff and astrochem are in your path, then place the `run_parallel_growth.py` and `cdinput.in` files into your desired directory and run:	

	python run_parallel_growth.py

or

	mpiexec python run_parallel_growth.py

An example sbatch submission script is also included.

## Inputs

The `run_parallel_growth.py` file is set to readin from `cdinput.in` file. The `cdinput.in` file should contain three sections labeled with a `#` at the start of the line: Model parameters (`# model`), Physical parameters (`# phys`), and intiail Abundances (`# abundances`).
Default parameters are given in the `chemdiff_io` module, and are included here.

### Model
```chmfile```

This should include the path to the chm file for astrochem. An absolute path should be used (i.e. starting from the root directory) so that even when astrochem is called in parallel from different directories it knows where to look.
*default:* `chmfile = /astrochem/networks/umist12.chm`

```pebfile```

This is the name of the file to output pebble abundances to. As Pebbles from they sequester some amount of ice on them, this file will keep track of the total abundance on pebbles as a function of time.
*default:* `pebfile = pebble_composition.out`

`ti`

Time to start the model in years.
*default:* `ti = 0`

`tf`

Time to end the model in years.
*default:* `tf = 1.e6`

`touts`

Times to save outputs in years. These should be multiples of your diffusion and/or chemical timestep.
*default:* `touts = [1e3, 2e3, 3e3, ..., 9e3, 1e4, 2e4, ..., 9e4, 1e5, 1.5e5, 2e5, 2.5e5, ..., 9.5e5, 1e6]`

`diff_dt`

Diffusion timestep in years. For stability, this should be less than $\frac{1}{100}(\alpha\Omega_K)^{-1}$, see Van Clepper et al. 2022 for details.
*default:* `diff_dt = 100`

`chem_dt`

Chemistry timestep in years. This should be greater than or equal to the diffusion timestep, and an integer multiple of the diffusion timestep.
*default:* `chem_dt = 500`


<!-- # Chemdiff

## Overview

Chemdiff is a python-based 1-D astrochemical code used to calculate the abundances of chemical species in protoplanetary disks. The code can account for vertical diffusion of species within the disk and pebble growth in addition to the astrochemical methods built on top of the Astrochem Code [Maret & Bergin (2015)][1]. See the [Astrochem documentation][2] for details on the ODE solver used to solve the chemical networks. In addition to the Astrochem solver, Chemdiff includes photodissociation self-shielding of CO, H2, and isotopologues, hydrogenation reactions of grains, xray ionization reactions, and Reactions with excited states of H2. Details of each of these reactions is given in the [Reaction rate calculations](#reaction-rate-calculations) section.

Chemdiff calculates diffusion and grain growth by defining a column at a given radial distance, $R$, from the central star with a given midplane temperature, $T_{mid}$, diffusion parameter, $\alpha$, a certain number of cells, $n_z$. Given this distance and midplane temperature, the column keplarian frequency, $\Omega$, sound speed, $c_s$, and scale height, $h$, are calculated (assuming a one solar mass star) using the formulas below. The column has a height of 5 scale heights, and is made of $n_z$ cells, each with a height of $\Delta z = 5/n_z$.

$$ \Omega = \sqrt{GM/R^3} $$

$$ c_s = \sqrt{\sigma T_{mid}/\bar{m}} $$

$$ h = c_s/\Omega $$

Within the column, cells are created with a given physical input parameters necessary to run astrochem.

## Physical parameter calculations

The cells are created with a given radial distance $r=R$, and height above the midplane $z_j = (j+0.5)\Delta z$, $0\leq j<n_z$. Each cell is also created with a given unattenuated UV flux, cosmic ionization rate, small grain size, dust-to-gas mass ratio, gas density, visual extinction, gas temperature, dust temperature, xray ionization rate, $CO$ column density, $H_2$ column density, $H$ column density, and abundance of each species in the cell. Many of these physical parameters have the same meaning as the parameters in the Astrochem documentation and are passed directly to the input or source files.

### UV flux

The unattenuated UV flux at each cell is given by the $\chi$ parameter and is kept constant for each cell in a column. This is the intensity of external UV field in Draine units. The default value is 1.

### Cosmic ionization rate

The ionization rate of molecular hydrogen due to cosmic rays in $s^{-1}$. The default value is $1.3\times10^{-17}s^{-1}$.

### Small grain size

This is the grain radius in microns for small grains that are used as surface area for grain surface reactions. As grains grow, the dust-to-gas mass ratio is adjusted, however this grain size parameter stays constant. The default value is $0.1 \mu m$.

### Dust-to-gas mass ratio

The small dust-to-gass mass ratio, $\epsilon$, value. This value is reduced with a characteristic timescale of $\tau_{grow} = 1/\Omega\epsilon$ for cells with $z \leq 1$. For more information on the dust growth process see the [Pebble Growth](#pebble-growth) section. The default initial value is 0.01.

### Gas density

The density of gas in $g$ $cm^{-3}$. This is calculated using hydrostatic equilibrium such that

$$ \rho_j = \rho_0 \exp(-z_j^2/2h^2) $$

where $\rho_0 = \Sigma / \sqrt{2\pi h^2}$ and $\Sigma$ is calculated according to [Aikawa & Herbst (2001)][3]

$$ \Sigma(R) = 7.2\times10^{23} m_h \left(\frac{R}{100AU}\right)^{-3/2} $$

The number density of hydrogen is also calculated using the relation

$$ n_{H,j} = 2\rho_j/\bar{m} $$

### Visual extinction

The visual extinction in magnitudes. This is calculated for each cell by numerically integrating down the column and using the relationship between visual and UV extinction $A_v = \tau_{UV}/3.02$ where

$$ \tau_{UV,j} = \sum_{i=j}^{n_z} \rho_i \kappa 100 \epsilon \Delta z$$

where $\kappa = 10^3 cm^2$ $g^{-1}$ is the UV dust opacity, and $\rho_i$ is the mass density of the $i$-$th$ cell. The $100 \epsilon$ term is used to scale the UV extinction with the dust-to-gas ratio. A value of $\kappa = 10^3 cm^2$ $g^{-1}$ is used as this gives good agreement with the relation $A_v = N_H/[1.8\times 10^{21} cm^{-2}mag^{-1}]$ [(Aikawa & Herbst, 2001)][3] when $\epsilon = 0.01$.

### Gas and dust temperature

The temperature of the gas and dust in each cell of the column. It is assumed the gas and dust are thermally coupled such that for the $j$-$th$ cell $T_{dust} = T_{gas} = T_j$. The temperature profile within the column is adopted from [Krijt (2018)][4], with $T_{atm} = 2T_{mid}$ and $z_q = 3h$

$$ T_j = T_{mid} + (T_{atm}-T_{mid})\left[\sin\left(\frac{\pi z_j}{2 z_q h}\right)\right]^4, \qquad z_j < z_q$$

$$ T_j = T_{atm}, \qquad \qquad \qquad \qquad \qquad \qquad \quad \; z_j \geq z_q $$

### X-ray ionization

The ionization rate of H due to x-rays. This is assumed to be zero throughout the column.

### Column Densities

The column densities of $CO$, $H_2$, and total $H$ are calculated using a numerical integration and the abundances, $X$, of each species

$$ N_{CO,j} = \sum_{i=j}^{n_z} X(CO)_i \Delta z $$

$$ N_{H_2,j} = \sum_{i=j}^{n_z} X(H_2)_i \Delta z $$

$$ N_{H,j} = \sum_{i=j}^{n_z} n_{H,i} \Delta z $$

Note that $N_{H,j}$ refers to the total column density of hydrogen atoms, and not the column density of molecular hydrogen.

### Abundances

Each cell has the time-dependent abundance of each species in the network stored as a dictionary of values.

## Pebble Growth

Pebble growth is simulating by removing small dust grains for cell with $z_j \leq h$ according to the growth timescale from [Birnstiel et al. (2012)][5]

$$ \tau_{grow} \approx \frac{1}{\Omega\epsilon} $$

With $\epsilon = \Sigma_d/\Sigma_g$ is the vertically integrated sut-to-gas mass ratio.

After each timestep from time $t$ to time $t+1$ of length $\Delta t$, the dust-to-gas mass ratio of the $j$-th cell is adjusted

$$ \epsilon_{j,t+1} - \epsilon_{j,t} = \Delta\epsilon_{j,t} = \frac{\epsilon_{j,t}}{\tau_{grow}} \Delta t = \epsilon_{j,t} f_{j,t}\Delta t$$

where $f_{j,t} = 1/\tau_{grow} = \Omega\epsilon_{j,t}$.





## Reaction rate calculations



[1]: <https://ui.adsabs.harvard.edu/abs/2015ascl.soft07010M/abstract> (Maret & Bergin, 2015)
[2]: <https://astrochem.readthedocs.io/en/latest/> (Astrochem documentation)
[3]: <https://ui.adsabs.harvard.edu/abs/2001A%26A...371.1107A/abstract> (Aikawa & Herbst, 2001)
[4]: <https://ui.adsabs.harvard.edu/abs/2018ApJ...864...78K/abstract> (Krijt et al., 2018)
[5]: <https://ui.adsabs.harvard.edu/abs/2012A%26A...539A.148B/abstract> (Birnstiel et al., 2012) -->