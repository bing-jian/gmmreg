%% Perform thin-plate spline warping
%% Input:
%%       landmarks:   source pts stored in nxd matrix.  
%%       parameters:  parameters in nxd matrix where first (d+1) rows are
%%       affine parameters corresponding to <1,x,y>
%% Output:
%%       warped_pts:  target pts in nxd matrix
%%       energy:      bending energy

function [warped_pts, bending_energy] = transform_by_tps(param, landmarks, ctrl_pts)
%%=====================================================================
%% $RCSfile: transform_by_tps.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-30 18:09:59 -0500 (Sun, 30 Nov 2008) $
%% $Revision: 116 $
%%=====================================================================
if (nargin==2)
    [n,d] = size(landmarks);
    [B,lambda] = compute_basis(landmarks);
    warped_pts = B*param;
    tps_param = param(d+2:n,:);
    bending_energy = trace(tps_param'*diag(lambda)*tps_param);
else
    [m,d] = size(landmarks);
    [n,d] = size(ctrl_pts);
    [K,U] = compute_kernel(ctrl_pts,landmarks);
    Pm = [ones(m,1) landmarks];
    Pn = [ones(n,1) ctrl_pts];
    PP = null(Pn'); B = [Pm U*PP]; 
    warped_pts = B*param;
    tps_param = param(d+2:n,:);
    bending_energy = trace(tps_param'*PP'*K*PP*tps_param);
end

