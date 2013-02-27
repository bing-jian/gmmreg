function [config] = initialize_config(model, scene, motion)

config.model = model;
config.scene = scene;
config.motion = motion;
% estimate the scale from the covariance matrix
[n,d] = size(model);
config.scale = power(det(model'*model/n), 1/(2^d));
config.display = 0;
config.init_param = [ ];
config.max_iter = 100;
config.normalize = 0;
switch lower(motion)
    case 'tps'
        interval = 5;
        config.ctrl_pts =  set_ctrl_pts(model, scene, interval);
        config.alpha = 1;
        config.beta = 0;
        config.opt_affine = 1;
        [n,d] = size(config.ctrl_pts); % number of points in model set
        config.init_tps = zeros(n-d-1,d);
        init_affine = repmat([zeros(1,d) 1],1,d);
        config.init_param = [init_affine zeros(1, d*n-d*(d+1))];
        config.init_affine = [ ];
    otherwise
        [x0,Lb,Ub] = set_bounds(motion);
        config.init_param = x0;
        config.Lb = Lb;
        config.Ub = Ub;
end

