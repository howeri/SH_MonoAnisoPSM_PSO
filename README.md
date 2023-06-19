# Efficient Monostatic Anisotropic Point Scatterer Model for High Performance Computing
Copyright (c) 2023 Eric Huang  
School of Electrical and Computer Engineering  
Georgia Institute of Technology  

The Matlab code associates with the paper "Efficient Monostatic Anisotropic Point Scatterer Model for High Performance Computing".

High performance computing (HPC) electromagnetic (EM) emulators are used to simulate real-time EM wave interactions between numerous radar targets. The Radar Cross Section (RCS) data stores the reflection profiles of radar targets as a table; however, the needed storage size increases with the frequency sampling density, the aspect angle sampling density, and the number of target types used in the emulator. The large quantity of data needed often exceeds storage capability and limits the feasibility of manipulation and representation of the objects. The spherical harmonic based monostatic anisotropic point scatterer model is proposed specifically for HPC EM interaction simulations where the scattering response can be computed with a finite impulse response (FIR) filter. An efficient algorithm constructing this model with large scale RCS data is discussed. The scatterer position and the reflection profile of each scatterer are solved using least squares methods and particle swarm optimization (PSO). In addition, the function evaluations in PSO are accelerated by taking advantage of the matrix structure, making the algorithm 22 times faster comparing to the naive approach. The results show that the point scatterer model can effectively represent the RCS data of a radar target.

For questions, queries and bug reports, please feel free to contact: huangeric@gatech.edu

## Examples:
To compute the monostatic point scatterer model of the aircraft model, run "pso_mono.m".

## Note:
The code is tested using Matlab R2019b. The spherical harmonic functions are generated by Javier Montalt Tordera (2021). Spherical Harmonics (https://github.com/jmontalt/harmonicY/releases/tag/v2.0.1), GitHub. Retrieved March 19, 2021. The aircraft geometry is obtained from R. Okada, “B787-8 dreamliner.” Online. https://grabcad.com/library/ b787-8-dreamliner-1 Accessed October 27, 2020.
