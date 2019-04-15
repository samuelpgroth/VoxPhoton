# VoxPhoton
Silicon photonics simulations with the volume integral equation method

To use, first enter the piecewise_constant subdirectory and run the build_const.m file. 
This will generate all the MEX files required.

The only example script currently available is example_Bragg_1500nm_to_1600nm.m 
This generates the transmission spectrum through a Bragg grating.

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
