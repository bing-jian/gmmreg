##=====================================================================
## $RCSfile: readme.txt,v $
## $Author: bing.jian $
## $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
## $Revision: 109 $
##=====================================================================

This directory contains the MATLAB code for the robust point-set
registration algorithm discribed in the ICCV'05 paper:

"Bing Jian and Baba C. Vemuri, 
A Robust Algorithm for Point Set Registration Using Mixture of Gaussians."

It is a part of software package which can be freely downloaded from
http://www.cise.ufl.edu/research/cvgmi/Software.php#gmmreg

Yet another Matlab implementation of this registration algorithm 
can be found at
http://www2.cs.man.ac.uk/~hous1/#downloads


Files in this "MATLAB" directory are organized as follows:

MATLAB/  

	gmmreg_demo.m
	    A Matlab test script of the C++ implementation.
	    See http://www.cise.ufl.edu/research/cvgmi/Software.php#gmmreg

	gmmreg_L2.m
	    The main entry for the MATLAB implementation.

	initialize_config.m
	    Generate the configuration struct used in gmmreg_L2.m
	
	GaussTransform/
	    MEX-files for implementing the GaussTransform

	registration/
	    Functions used in the Matlab implementation of the GMMReg algorithm, 
	    requiring 'GaussTransform' and the optimization toolbox.

	auxiliary/
	    Some supporting functions, mostly for displaying the results.
		

Use 'addpath(genpath(pwd))' to add the "MATLAB" directory 
and its subdirectories to the MATLAB search path.

Please addresse any questions to bing.jian@gmail.com

	


     
    	