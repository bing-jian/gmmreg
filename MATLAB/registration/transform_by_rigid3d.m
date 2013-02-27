% [result] = transform_by_rigid3d(pointset, param)
% perform a 3D rigid transform on a pointset and
% return the transformed pointset
% Note that here 3D rigid transform is parametrized by 7 numbers
% [quaternion(1 by 4) translation(1 by 3)]
%
% See also: quaternion2rotation, transform_by_rigid2d,
% transform_by_affine3d

function [result] = transform_by_rigid3d(pointset, param)
%%=====================================================================
%% $RCSfile: transform_by_rigid3d.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-30 18:09:59 -0500 (Sun, 30 Nov 2008) $
%% $Revision: 116 $
%%=====================================================================
n = size(pointset,1);
r = quaternion2rotation(param(1:4));
t = ones(n,1)*param(5:7);
result = pointset*r' + t;

