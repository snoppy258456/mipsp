# mipsp
Program for Monte Carlo simulations and Free Energy calculations of molecular imprinting as described in T. Curk, J. Dobnikar and D. Frenkel, Rational design of molecular imprinting, Soft Matter (2015), DOI: 10.1039/C5SM02144H

The simulation code is written in Fortran 90 and can be used freely for Free energy callculations of analytes in imprinted cavities, or for simulating imprinting and rebinding process. 

1. The free energy calculation program is found in FreeE_calc folder. In the input file we need to provide the geometry of the cavity (positions of imprinted ligands - functional monomers) and the analyte geometry (positions of receptors - binding sites). We also need to provide the ligand-receptor interaction matrix, i.e. which receptor binds to which ligand with binding energy. The program then uses Monte Carlo sampling to compute the Free energy of bidning the analyte, or the "Avidity" equlibrium association constant.

2. The program in SIMS folder can be used to simulate both the imprinting process and subsequent analyte rebinding of a large system.