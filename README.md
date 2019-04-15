# VoxPhoton
Silicon photonics simulations with the volume integral equation method

To use, first enter the piecewise_constant subdirectory and run the build_const.m file. 
This will generate all the MEX files required.

The only example script currently available is example_Bragg_1500nm_to_1600nm.m 
This generates the transmission spectrum through a Bragg grating.

## To do:
1. Scripts to generate preconditioning results from paper [3], for straight waveguide, Bragg grating, disk resonator 
(frequency-domain simulations)

2. Scripts to generate adiabatic absorber results from [1]

2. Broadband disk and ring resonator simulations

3. Directional coupler

4. Y-branch splitter

* Write documentation for everything!

## Publications
There are currently three publications associated with this work:

[1] Tambova, A., Groth, S. P., White, J. K., & Polimeridis, A. G. (2018). 
    Adiabatic absorbers in photonics simulations with the volume integral equation method. 
    Journal of Lightwave Technology, 36(17), 3765-3777.
    
[2] Groth, S. P., Polimeridis, A. G., Tambova, A., & White, J. K. (2019). 
    Circulant preconditioning in the volume integral equation method for silicon photonics. 
    arXiv preprint arXiv:1903.07884.
    
[3] Groth, S., White, J., & Polimeridis, A. (2018, March). 
    Circulant preconditioning in the volume integral equation method for nanophotonics. 
    In 2018 International Applied Computational Electromagnetics Society Symposium (ACES) (pp. 1-2). IEEE.
    
The first details the truncation of the photonics structures via adiabatic absorbers. 
The second two discuss the acceleration of the individual frequency solves via circulant preconditioning.
A paper detailing the model-order reduction required for the broadband simulations is currently in preparation.
