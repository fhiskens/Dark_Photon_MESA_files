This directory contains the MESA files and data required to simulate stars experiencing energy-loss to dark photons. They were written/calculated as part of the paper `Constraining dark photons with self-consistent simulations of globular cluster stars' (https://arxiv.org/abs/2306.13335) by Matthew J. Dolan, FJH and Raymond R. Volkas.

The Data directory contains tabulated values of the function G (Eq. A2). The structure of this directory and its subdirectories is described in Sec. A.1.

The simulations performed in this paper were broken up into three stages:

1. pre-MS to helium ignition,
2. horizontal branch, and
3. asymptotic giant branch.
 
The files required to perform each one of these stages have been stored in different MESA work directories, which we call the 'base' directories. Specifically, stage 1 is performed by 'base_RGB', stage 2 by 'base_SO' and stage 3 by 'base_AGB'. Because our simulations were performed on a cluster, input physics was specified in two different files. The fixed (base) input physics, which did not change between different simulations of each phase, is contained in the 'inlist_lowinter' file. Input specifications which are varied between different simulations are instead housed in a separate inlist, called 'inlist_cluster'. Examples of these can be found in each of the 'base directories'.

Energy-loss to dark photons has been included by modifying the module responsible for calculating energy-loss to thermal neturinos (neu). These modifications are housed in 'src/run_star_extras.f' file within each base directory. The 'run_star_extras.f' files in each of the base directories are identical.

If one wishes to use these files, one must specify a number of file paths. The corresponding files and lines are written below.

1. base_*/inlist: lines 10, 22, 33
2. base_*/inlist_lowinter: lines 12 and 36
3. base_*/src/run_star_extras.f: line 303.

Several features in the run_star_extras file and inlists originate from the MIST input physics (J. Choi, A. Dotter, C. Conroy, M. Cantiello, B. Paxton and B.D. Johnson, Mesa Isochrones and Stellar Tracks (MIST). I. Solar-scaled Models, Astrophys. J. 823 (2016) 102 [1604.08592]).
