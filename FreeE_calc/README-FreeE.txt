********************************************************************
************** FREE ENERGY CALCULATIONS README *********************
********************************************************************

Molecularly imprinted polymers code, you need:
mcmips_FE.f90  	   -program
input-par.dat  	   -input file
 
The code can be compiled using a a standard Fortran compiler, e.g. 
gfortran mcmips_FE.f90 -O3


In the input file we need to provide the geometry of the cavity (positions of imprinted ligands - functional monomers) and the analyte geometry (positions of receptors - binding sites). We also need to provide the ligand-receptor interaction matrix, i.e. which receptor binds to which ligand with binding energy. The program then uses Monte Carlo sampling to compute the Free energy of binding the analyte, or the "Avidity" equilibrium association constant. When the program is run it will write the results to screen. The current parameters in input-par.dat correspond to calculation of binding free energy of enantiomers shown on Figure S5 in the SI of the paper.

In the program we can assign "reactive" types called RX. These can be ligands or receptors the program doesn't distinguish between them implicitly. We need to choose the number of distinct RX types, their positions in the cavity (ligands) and on the analyte (receptors) and we provide an interaction energy matrix.

In the input-par.dat file provided we wish to calculate binding free energies of enantiomers. We have 8 RX types in total, we assign 4 to be ligands(anchors) and the next 4 to be receptors on the analyte. We therefore have 8x8 interaction matrix. An interaction of +100 is deemed no interaction. 

In the following I will describe the parameters in the input file: 
******************************************************
************ INPUT PARAMETERS DESCRIPTION ************
-----------------------------------------------------------
10.0,10.0,10.0  box size x,y,z  - system size should be large enough such that colloid cannot form bonds through periodic boundary, 10x10x10 is large enough unless very soft matrices are simulated (kspring<1) in which case it should be increased.
1 8	        # of colloid species, # of rx species - colloid species should be 1, # of RX species can be changed, it specifies the total number of distinct (receptors + ligands) "reactive" types in the system.
0	     	init # of colloids (for each specie) - should be always 0
1	     	max # of colloids (per species, one number for all) - should be always 1 for Free energy calculations	  
1 1 1 1 0 0 0 0		tot # of anchors  per rx specie - here we specify the number of anchor per each RX type, with 8 types we need to provide 8 integers
4		# of rx sites per colloid - total number of all RX sites on colloid, all types together
100 0.0		kspring,dspring ; spring stiffness, length - ligand-receptor spring constant, dspring is the zero energy spring length and should be 0	
-15		mu, chemical potential for each colloid species log(vol_fraction) - chemical potential of analytes, for Free energy calculations should be in the range -20 - -5, depending on other parameters to get good statistics.	
200000000	tot # of MC cycles - total simulation time in Monte Carlo cycles
10000000	output configuration timer - frequency of writing results to screen and writing a configuration file. Should be large enough >~1e6 so that Free Energy calculation is averaged over sufficient number of cycles.   
------------------------------------------------------------
0.5		colloid radius (one for all species) - should be 0.5, so that radius is 1 which is defined as the unit length scale in simulations
0.1		max hop col - maximum MC displacement move for colloids(analytes), 0.1 is a good choice but could be optimised further. 
0.2		max_rot_col (actually sin(.)) - maximum MC rotation move for colloids(analytes), 0.2 is a good choice but could be optimised further. 
0.2 0.4		fraction of insert/delete, create/destroy bond moves ; else hop moves - These numbers designate the relative probabilities of performing different MC moves: 1. colloid exchange(insert/delete), 2. bond formation/breaking, 3. move colloid. Only the first two should be specified and the third one is calculated as the three have to sum to 1.0 
1		bond create/destroy cutoff - how far away ligand/receptor bond can be formed destroyed, 1.0 is a good choice, but can be increased particularly if very soft matrices (kspring < 5) are simulated. This number does not affect equilibrium properties (but obviously must be >0) as analytes can still diffuse further away when the bond is formed.
8		mnpic+1 -- max # of colloids in each cell + 1  - maximum number of particles in the cell list, this number has to be provided for book-keeping cell list array initialization, 8 is a good number as the cell size is given by the colloid diameter. If this is too low the program will crash or give nonsensical results.
-----------------------------------------------------------
.true.		random seed - initialize Monte Carlo with a random seed, should be .true. unless debugging
.false.		mobile anchors - ligands are free to move in the system, should be .false. for Free energy calculations
.true. 		imprint_anchors, otherwise random positions in the box - ligand anchors have specified positions (below) if .true., otherwise ligand anchors are initiated with random positions in the simulation box
.false.	.false.	read initial configuration, read only anchor positions - should be both .true. if the system anchor position are to be initialized from an initial configuration file. For free energy calculations both should be .false. 
'imp_mu2050_ncol1020_out1012.xyz'	initial configuration file name - self explanatory
'test'	output file name - the forename of output configuration files
-----------------------------------------------------------
RX BOND ENE INTERACTION MATRIX	 - interaction energy matrix between rx types, here 1-4 types are ligands and 5-8 types are anchors, Here we wished that 1 binds to 5, 2 to 6, 3 to 7 and 4 to 8. Hence, these entries are filled (matrix must be symmetric!) other values are set to a large positive value which practically means no binding.   	
100 100 100 100 -3 100 100 100
100 100 100 100 100 -3 100 100
100 100 100 100 100 100 -3 100
100 100 100 100 100 100 100 -3
-3 100 100 100 100 100 100 100
100 -3 100 100 100 100 100 100
100 100 -3 100 100 100 100 100
100 100 100 -3 100 100 100 100
-----------------------------------------------------------
RX SPECIE ON COLLOID GEOMETRY (specie,vector3) normalised to particle radius - receptors on analytes normalised vector positions, the first number is the RX type, the following three specify the vector. Each receptor must be in its own line, the RX types must follow in a monotonic sequence, from lowest to highest (for example, 5,6,7,8). 
5  -1.0  0.0   0.0      
6  0.0   1.0   0.0
7  1.0   0.0   0.0
8  0.0  -1.0   0.0
-----------------------------------------------------------
IMPRINTED ANCHOR POSITIONS vectors  - the vector positions of imprinted ligand anchors in the system, 3 number for each vector. Vectors are written on the same line, and RX types corresponding to these vectors follow from lowest to highest depending on the specified number of anchors per RX type above 
2.0 2.5 2.5  2.5 3.0 2.5  3.0 2.5 2.5  2.5 2.0 2.5  
-----------------------------------------------------------

********* END INPUT PARAMETERS DESCRIPTION ***********
******************************************************

******************************************************
******************* PROGRAM OUTPUT *******************
The program output on screen should look something like this:
First simulation parameters are written to screen. Than simulation starts and every cycle steps the system status is written to screen, for example in the example below we calculated Free energy of binding to be about 7.9kT

 ================ MC MIPS SIM ================
 =============================================
 ================= INPUT PAR =================
 boxsize     10.0000000000000        10.0000000000000        10.0000000000000     
 init # of colloids             0
 max # of colloids             1
 # of anchors(ligands)             1           1           1           1           0           0           0           0
 fraction of insert/delete, bond create/destroy moves    0.200000000000000       0.400000000000000     
 activity    4.539992976248485E-005
 tot n cycles              100000000
 nout                10000000
 ---------------------------------------------
 rcol    0.500000000000000     
 max hop/rot col    0.100000000000000       0.200000000000000     
 # of rx sites per col            4
 k spring    10.0000000000000     
 single bond energy rxobondene(2,1)   99.3029379584089     
 make/break bond distance cut off    1.00000000000000     
 max # of particles in each cell             8
 size of cell -colloids:   1.00000000000000        1.00000000000000        1.00000000000000     
 size of cell -anchors:   1.00000000000000        1.00000000000000        1.00000000000000     
 ---------------------------------------------
 random seed   T
 mobile anchors  F
 read initial conf F F
 init conf filename  imp_mu2050_ncol1020_out1012.xyz                   
 outfilename  test                                              
 ============== END INPUT PAR ================
 =============================================
 ================= START SIM =================

icycle  10000000, exe time (s)    13.2985
average avncol:   0.151509900  avnbcol:   0.112815000
tot ene old:     0.000   tot ene:     0.000
moveacc:0.09518   excacc:0.02959   bondacc:0.00533
Analyte binding associaton constant K_A =      2928.06
Analyte binding Free Energy F_cav =  -7.982095

icycle  20000000, exe time (s)    13.1928
average avncol:   0.142210500  avnbcol:   0.103095800
tot ene old:     0.000   tot ene:     0.000
moveacc:0.09130   excacc:0.03027   bondacc:0.00499
Analyte binding associaton constant K_A =      2646.81
Analyte binding Free Energy F_cav =  -7.881109

