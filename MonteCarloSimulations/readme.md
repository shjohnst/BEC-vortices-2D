# Monte Carlo simulation code for "Evolution of large-scale flow from turbulence in a two-dimensional superfluid"
[(View on GitHub)](https://github.com/shjohnst/BEC-vortices-2D)

Summary of scripts needed to perform point-vortex Monte Carlo simulations.

**Author:** Andrew J. Groszek
**Date:**   18th October 2018


## Scripts
1. `MCMC_pointvortices.m`
	Runs the Markov Chain Monte Carlo simulations in the disk-shaped trap for a chosen vortex number and set of temperatures.
	The primary output is an array of multiple vortex configurations at each temperature.

2. `analyse_data.m`
	Takes the output of `MCMC_pointvortices.m` and classifies the vortex configurations.


## Functions (needed for scripts to run)
3. `calc_pointvortex_hamiltonian_disk`
        Calculates the point-vortex energy for a chosen vortex configuration.

4. `classify_vortices`
        Applies our vortex classification algorithm to assign vortices in a chosen configuration as clustered vortices, dipole vortices or free vortices. For details on how the algorithm works, see Section 6 of Valani et al., 'Einstein-Bose condensation of Onsager vortices', [New Journal of Physics 20, 053038 (2018)](https://doi.org/10.1088/1367-2630/aac0bb).