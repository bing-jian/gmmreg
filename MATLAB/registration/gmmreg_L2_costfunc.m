function [f,g] = gmmreg_L2_costfunc(param, config)
%%=====================================================================
%% $RCSfile: gmmreg_L2_costfunc.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================
model = config.model;
scene = config.scene;
motion = config.motion;
scale = config.scale;
[transformed_model] = transform_pointset(model, motion, param);
switch lower(config.motion)
    case 'rigid2d'
        [f, grad] = rigid_costfunc(transformed_model, scene, scale);
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        grad = grad*model;
        theta = param(3);
        r = [-sin(theta) -cos(theta);
             cos(theta)  -sin(theta)];
        g(3) = sum(sum(grad.*r));
    case 'rigid3d'
       [f,grad] = rigid_costfunc(transformed_model, scene, scale);
        [r,gq] = quaternion2rotation(param(1:4));
        grad = grad';
        gm = grad*model; 
        g(1) = sum(sum(gm.*gq{1}));
        g(2) = sum(sum(gm.*gq{2}));
        g(3) = sum(sum(gm.*gq{3}));
        g(4) = sum(sum(gm.*gq{4}));        
        g(5) = sum(grad(1,:));
        g(6) = sum(grad(2,:));
        g(7) = sum(grad(3,:));
    case 'affine2d'
        [f,grad] = general_costfunc(transformed_model, scene, scale);
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        g(3:6) = reshape(grad*model,1,4);
    case 'affine3d'
        [f,grad] = general_costfunc(transformed_model, scene, scale);
        grad = grad';
        g(1) = sum(grad(1,:));
        g(2) = sum(grad(2,:));
        g(3) = sum(grad(3,:));
        g(4:12) = reshape(grad*model,1,9);
    otherwise
        error('Unknown motion type');
end;


function [f, g] = rigid_costfunc(A, B, scale)
[f, g] =  GaussTransform(A,B,scale);
f = -f; g = -g;


function [f, g] = general_costfunc(A, B, scale)
[f1, g1] = GaussTransform(A,A,scale);
[f2, g2] = GaussTransform(A,B,scale);
f =  f1 - 2*f2;
g = 2*g1 - 2*g2;