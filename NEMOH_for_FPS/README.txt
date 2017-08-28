Readme to excecute NEMOH boundary element method code.
1. To compute frequency domain coefficients run the 'axiSymmesh.m' script
	1A. Use 36 for the number of point for angular discretisation
	1B. Enter 'storage' for the directory name
	1C. Enter 0.05 for the vertical position of gravity center
	1d. Enter 1000 for the target number of panels
2. Run the script 'runningNemoh.m'
3. Run 'GenerateSSreal.m'
4. Move data file 'ABCD_r_SSreal.dat' to the working directory.

Please not that all examples already have the NEMOH data that was used to generate the
figures in the paper.
