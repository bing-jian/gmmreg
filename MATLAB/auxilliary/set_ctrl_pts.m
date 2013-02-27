function [ctrl_pts] = set_ctrl_pts(M,S,interval)
%%=====================================================================
%% Project:   Point Set Registration using Gaussian Mixture Model
%% Module:    $RCSfile: set_ctrl_pts.m,v $
%% Language:  MATLAB
%% Author:    $Author: bing.jian $
%% Date:      $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% Version:   $Revision: 109 $
%%=====================================================================

[axis_limits] = determine_border(M, S);
x_min = axis_limits(1,1);
x_max = axis_limits(1,2);
y_min = axis_limits(2,1);
y_max = axis_limits(2,2);

[m,d] = size(M);
if (nargin<3)
    interval = 5;
end
if (d==2)
    [x,y] = ndgrid(linspace(x_min,x_max,interval), linspace(y_min,y_max,interval));
    n_pts = interval*interval;
    ctrl_pts = [reshape(x,n_pts,1) reshape(y,n_pts,1)];
end
if (d==3)
    z_min = axis_limits(3,1);
    z_max = axis_limits(3,2);
    [x,y,z] = ndgrid(linspace(x_min,x_max,interval), linspace(y_min,y_max,interval),linspace(z_min,z_max,interval));
    n_pts = interval*interval*interval;
    ctrl_pts = [reshape(x,n_pts,1) reshape(y,n_pts,1) reshape(z,n_pts,1)];    
end