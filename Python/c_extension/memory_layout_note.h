/*
 Note:  the assumed memory layout is point-wise contiguous, 
 i.e. the j-th coordinate of the i-th point is at location  
 (i*d + j) where d is the dimensionality. Thus, the MATLAB 
 mex-function expects an input matrix where each point is a 
 d-dimensional column vector as MATLAB is column-major while 
 the input array to the NumPy should place each point as a 
 d-dimensional row vector since NumPy is row-major.
*/
