********************************************************************
******************* MIPS SIMULATION README *************************
********************************************************************
by Tine Curk, tc387@cam

Molecularly imprinted polymers code, you need:
mcmips_SIMS.f90    -program
input-par.dat  	   -input file
 
The code can be compiled using a a standard Fortran compiler, e.g. 
gfortran mcmips_SIMS.f90 -O3

This code is essentially the same as the one provided for free energy calculations, only the screen output is different and provides a calculated binding affinity of an analaytes/MIP system rather than the binding free energy of an analyte into a single cavity.  

We can use the code either to simulate imprinting process (folder IMP), in this case the number of templates and ligands is fixed in the system and ligands are free to move around the system (effectively a canonical Monte Carlo simulation). Or we can use the code to simulate re-binding of analytes into a system with fixed ligand anchor positions (folder MIP) (Grand canonical Monte Carlo). Ligand anchors can be either randomly distributed (NIP) or can be initialised with an .xyz configuration file obtained from a  previous imprinting simulation. 

The parameters of the input file are already described in the FreeE_calc folder. Here I will only highlight the parameters relevant to the imprinting/rebinding simulation, other parameters should be the same in both cases.In the following "IMP" will designate imprinting simulations and "MIP" will designate analyte rebinding to mips simulation parameters set up.

======================================================
************ INPUT PARAMETERS DESCRIPTION ************
-----------------------------------------------------------
100	     	max # of colloids (per species, one number for all) - for IMP this number should be the total number of templates in the system. For MIP the number should be very large (~1000) , larger than would be ever reached, as it designates the maximum number of analytes in the system. 
100 0.0		kspring,dspring ; spring stiffness, length - ligand-receptor spring constant, for IMP it should be large (~1000) and it designates the fidelity of imprinted site (cavity) as compared to the template. For MIP it should be whatever matrix stiffness we wish to simulate.
-5		mu, chemical potential for each colloid species log(vol_fraction) - for IMP it should be a huge number (~100) in order to push the number of templates in the system to max # of colloids, so that we effectively perform a canonical simulations. For IMP it should be whatever the chemical potential of analytes we wish to simulate.
------------------------------------------------------------
.false.		mobile anchors - ligands are free to move in the system, should be .true. for IMP, and .false. for MIP calculations
.false. 	imprint_anchors, otherwise random positions in the box - always .false. 
.false.	.false.	read initial configuration, read only anchor positions - should be both .true. for MIP and both .false. for IMP. The program will read only ligand positions from the .xyz configuration file.
'imp_mu2050_ncol1020_out1012.xyz'	initial configuration file name - self explanatory, the name of initial configuration file if we perform MIP simulations. 
-----------------------------------------------------------
IMPRINTED ANCHOR POSITIONS vectors  - this entry is not needed if "imprint anchors" above is set to .false.

********** END INPUT PARAMETERS DESCRIPTION ***********
======================================================



======================================================
**************** PROGRAM OUTPUT IMP *****************
The program output should look something like this for IMP simulation. First parameters are written to screen, than simulation starts and every number of icycles system status is written to screen: 
icycle - MC cycle number, exe time: execution time in seconds for this cycle
average number of colloids(analytes) in the system, average number of bound analytes
energy of the system (old and new energy should be the same, otherwise there is a bug in the code)
acceptance ratio of different MC moves (move, inert/delete, bond form/break)
Binding affinity calculated from the number of bound analytes, system volume and analyte chemical potential
In the case of imprinting (IMP) simulation binding affinity is meaningless. 

Example here for an imprinting simulation:

================ MC MIPS SIM ================
 =============================================
 ================= INPUT PAR =================
 boxsize     12.6000000000000        12.6000000000000        12.6000000000000     
 init # of colloids             0
 max # of colloids            20
 # of anchors(ligands)            40           0
 fraction of insert/delete, bond create/destroy moves    0.200000000000000       0.400000000000000     
 activity    5.184705528587072E+021
 tot n cycles              100000000
 nout                10000000
 ---------------------------------------------
 rcol    0.500000000000000     
 max hop/rot col    0.100000000000000       0.200000000000000     
 # of rx sites per col            2
 k spring    1000.00000000000     
 single bond energy rxobondene(2,1)  -14.5048173205732     
 make/break bond distance cut off    1.00000000000000     
 max # of particles in each cell             8
 size of cell -colloids:   1.05000000000000        1.05000000000000        1.05000000000000     
 size of cell -anchors:   1.05000000000000        1.05000000000000        1.05000000000000     
 ---------------------------------------------
 random seed   T
 mobile anchors  T
 read initial conf F F
 init conf filename  imprinted_ncol1020_nanc1040.xyz                   
 outfilename  imprinted                                         
 ============== END INPUT PAR ================
 =============================================
 ================= START SIM =================

icycle  10000000, exe time (s)    15.5898
average avncol:  19.999771400  avnbcol:  15.918987200
tot ene old:  -386.639   tot ene:  -386.639
moveacc:0.00644   excacc:0.00000   bondacc:0.00000
Average number of bound analytes   15.9190
Binding Affinity    0.0000

icycle  20000000, exe time (s)    15.4175
average avncol:  20.000000000  avnbcol:  19.670586200
tot ene old:  -382.290   tot ene:  -382.290
moveacc:0.00273   excacc:0.00000   bondacc:0.00000
Average number of bound analytes   19.6706
Binding Affinity    0.0000

icycle  30000000, exe time (s)    15.4474
average avncol:  20.000000000  avnbcol:  19.541982900
tot ene old:  -432.919   tot ene:  -432.919
moveacc:0.00306   excacc:0.00000   bondacc:0.00000
Average number of bound analytes   19.5420
Binding Affinity    0.0000


************** END PROGRAM OUTPUT IMP ****************
======================================================



======================================================
**************** PROGRAM OUTPUT MIP ****************
The program output should look something like this for MIP simulation..

================ MC MIPS SIM ================
 =============================================
 ================= INPUT PAR =================
 boxsize     12.6000000000000        12.6000000000000        12.6000000000000     
 init # of colloids             0
 max # of colloids          2000
 # of anchors(ligands)            40           0
 fraction of insert/delete, bond create/destroy moves    0.200000000000000       0.400000000000000     
 activity    1.013009359863071E-005
 tot n cycles             1000000000
 nout               100000000
 ---------------------------------------------
 rcol    0.500000000000000     
 max hop/rot col    0.100000000000000       0.200000000000000     
 # of rx sites per col            2
 k spring    100.000000000000     
 single bond energy rxobondene(2,1)  -11.0509396810821     
 make/break bond distance cut off    1.00000000000000     
 max # of particles in each cell             8
 size of cell -colloids:   1.05000000000000        1.05000000000000        1.05000000000000     
 size of cell -anchors:   1.05000000000000        1.05000000000000        1.05000000000000     
 ---------------------------------------------
 random seed   T
 mobile anchors  F
 read initial conf T T
 init conf filename  imprinted_ncol1020_nanc1040.xyz                   
 outfilename  mip                                               
 ============== END INPUT PAR ================
 =============================================
 ================= START SIM =================

icycle 100000000, exe time (s)   161.5007
average avncol:   3.675517140  avnbcol:   3.655524680
tot ene old:   -70.284   tot ene:   -70.284
moveacc:0.05655   excacc:0.00238   bondacc:0.00001
Average number of bound analytes    3.6555
Binding Affinity  180.3951

icycle 200000000, exe time (s)   162.6105
average avncol:   4.190091550  avnbcol:   4.170100000
tot ene old:   -64.513   tot ene:   -64.513
moveacc:0.05055   excacc:0.00211   bondacc:0.00001
Average number of bound analytes    4.1701
Binding Affinity  205.7886





