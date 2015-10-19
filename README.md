# mipsp
Program for Monte Carlo simulations and Free Energy calculations of molecular imprinting as described in T. Curk, J. Dobnikar and D. Frenkel, Rational design of molecularly imprinted polymers, Soft Matter (2015), DOI: 10.1039/C5SM02144H

The simulation code is written in Fortran 90 and can be used freely for Free energy calculations of analytes in imprinted cavities, or for simulating imprinting and rebinding process. The guide to using the program is found in respective folders. Any enquiry, questions or comments are welcome, please send me an email to tcurk5@gmail.com and I'll be happy to help.

The program can be used for:

1. Free energy calculations. The Free energy calculation program is found in FreeE_calc folder. In the input file we need to provide the geometry of the cavity (positions of imprinted ligands - functional monomers) and the analyte geometry (positions of receptors - binding sites). We also need to provide the ligand-receptor interaction matrix, i.e. which receptor binds to which ligand with binding energy. The program then uses Monte Carlo sampling to compute the Free energy of binding the analyte, or the "Avidity" equilibrium association constant.

2. Simulation of imprinting/rebinding. The program in SIMS folder can be used to simulate both the imprinting process and subsequent analyte rebinding of a larger system.
