% Perform a spatial tranformation on a given pointset
% motion:  the motion model represented by string, can be
%         'rigid2d',    'rigid3d', 'affine2d',  'affine3d', 'tps'
% parameter: a row vector
function [transformed_pointset] = transform_pointset(pointset, motion, parameter, varargin)
%%=====================================================================
%% $RCSfile: transform_pointset.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================

switch lower(motion)
    case 'rigid2d'
        transformed_pointset = transform_by_rigid2d(pointset, parameter);
    case 'rigid3d'
        transformed_pointset = transform_by_rigid3d(pointset, parameter);
    case 'affine2d'
        transformed_pointset = transform_by_affine2d(pointset, parameter);
    case 'affine3d'
        transformed_pointset = transform_by_affine3d(pointset, parameter);
    case 'tps'
        ctrl_pts = varargin{1};
        init_affine = varargin{2};
        [n,d] = size(ctrl_pts);
        p = reshape([init_affine parameter],d,n); p = p'; 
        transformed_pointset = transform_by_tps(p, pointset, ctrl_pts);
    otherwise
        error('Unknown motion type');
end;


