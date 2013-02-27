%function [param, tt] = GMMReg(model, scene, scale, motion, display, init_param);
%   'model' and 'scene'  are two point sets
%   'scale' is a free scalar parameter
%   'motion':  the transformation model, can be
%         ['rigid2d', 'rigid3d', 'affine2d', 'affine3d']
%         The default motion model is 'rigid2d' or 'rigid3d' depending on
%         the input dimension
%   'display': display the intermediate steps or not.
%   'init_param':  initial parameter

function [param, transformed_model, history, config] = gmmreg_L2(config)
%%=====================================================================
%% $Author: bing.jian $
%% $Date: 2009-02-10 02:13:49 -0500 (Tue, 10 Feb 2009) $
%% $Revision: 121 $
%%=====================================================================

% todo: use the statgetargs() in statistics toolbox to process parameter name/value pairs
% Set up shared variables with OUTFUN
history.x = [ ];
history.fval = [ ];
if nargin<1
    error('Usage: gmmreg_L2(config)');
end
[n,d] = size(config.model); % number of points in model set
if (d~=2)&&(d~=3)
    error('The current program only deals with 2D or 3D point sets.');
end

options = optimset( 'display','iter', 'LargeScale','off','GradObj','on', 'TolFun',1e-010, 'TolX',1e-010, 'TolCon', 1e-10);
options = optimset(options, 'outputfcn',@outfun);
options = optimset(options, 'MaxFunEvals', config.max_iter);
options = optimset(options, 'GradObj', 'on');

tic
switch lower(config.motion)
    case 'tps'
        scene = config.scene;
        scale = config.scale;
        alpha = config.alpha;
        beta = config.beta;

        [n,d] = size(config.ctrl_pts);
        [m,d] = size(config.model);
        [K,U] = compute_kernel(config.ctrl_pts, config.model);
        Pm = [ones(m,1) config.model];
        Pn = [ones(n,1) config.ctrl_pts];
        PP = null(Pn');  % or use qr(Pn)
        basis = [Pm U*PP];
        kernel = PP'*K*PP;

        init_tps = config.init_tps;  % it should always be of size d*(n-d-1)
        if isempty(config.init_affine)
            % for your convenience, [] implies default affine
            config.init_affine = repmat([zeros(1,d) 1],1,d);
        end
        if config.opt_affine % optimize both affine and tps
            init_affine = [ ];
            x0 = [config.init_affine init_tps(end+1-d*(n-d-1):end)];
        else % optimize tps only
            init_affine = config.init_affine;
            x0 = init_tps(end+1-d*(n-d-1):end);
        end
        param = fminunc(@(x)gmmreg_L2_tps_costfunc(x, init_affine, basis, kernel, scene, scale, alpha, beta, n, d), x0,  options);
        transformed_model = transform_pointset(config.model, config.motion, param, config.ctrl_pts, init_affine);
        if config.opt_affine
            config.init_tps = param(end+1-d*(n-d-1):end);
            config.init_affine = param(1:d*(d+1));
        else
            config.init_tps = param;
        end
    otherwise
        x0 = config.init_param;
        param = fmincon(@gmmreg_L2_costfunc, x0, [ ],[ ],[ ],[ ], config.Lb, config.Ub, [ ], options, config);
        transformed_model = transform_pointset(config.model, config.motion, param);
        config.init_param = param;
end
toc

    function stop = outfun(x,optimValues,state,varargin)
     stop = false;
     switch state
         case 'init'
             if config.display>0
               set(gca,'FontSize',16);
             end
         case 'iter'
               history.fval = [history.fval; optimValues.fval];
               history.x = [history.x; reshape(x,1,length(x))];
               if config.display>0
                   hold off
                   switch lower(config.motion)
                       case 'tps'
                           transformed = transform_pointset(config.model, config.motion, x, config.ctrl_pts,init_affine);
                       otherwise
                           transformed = transform_pointset(config.model, config.motion, x);
                   end
                   dist = L2_distance(transformed,config.scene,config.scale);
                   DisplayPoints(transformed,config.scene,d);
                   title(sprintf('L2distance: %f',dist));
                   drawnow;
               end
         case 'done'
              %hold off
         otherwise
     end
    end


end



function [dist] = L2_distance(model, scene, scale)
    dist = GaussTransform(model,model,scale) + GaussTransform(scene,scene,scale) - 2*GaussTransform(model,scene,scale);
end


