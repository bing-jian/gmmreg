% [f, grad] = GaussTransform(A,B,scale)
% The inner product between two spherical Gaussian mixtures computed using the Gauss Transform.
% The centers of the two mixtures are given in terms of two point sets A and B (of same dimension d)
% represented by an mxd matrix and an nxd matrix, respectively.
% It is assumed that all the components have the same covariance matrix represented by 
% a scale parameter (scale).  Also, in each mixture, all the components are equally weighted.

%%  $Author: bjian $
%%  $Date: 2008/04/06 03:59:15 $
%%  $Revision: 1.1 $  
function [f,g] = GaussTransform(A,B,scale)	

if exist('mex_GaussTransform','file')
    [f,g] = mex_GaussTransform(A',B',scale);
    g = g';
else
    message = ['Precompiled GaussTransform module not found.\n' ...
        'If the corresponding MEX-functions exist, run the following command:\n' ...
        'mex mex_GaussTransform.c GaussTransform.c -output mex_GaussTransform'];
    message_id = 'MATLAB:MEXNotFound';
    error (message_id, message);
    
end
