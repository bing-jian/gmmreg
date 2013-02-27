%% compute the kernel and basis matrix for thin-plate splines
%% reference:  Landmark-based Image Analysis, Karl Rohr, p195
function [K] = tps_compute_kernel(x,z)
%%=====================================================================
%% $RCSfile: tps_compute_kernel.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-30 18:09:59 -0500 (Sun, 30 Nov 2008) $
%% $Revision: 116 $
%%=====================================================================
[n, d] = size (x);
[m, d] = size (z);

% calc. the K matrix.
% 2D: K = r^2 * log r
% 3D: K = -r
K = zeros (n,m);

for i=1:d
  tmp = x(:,i) * ones(1,m) - ones(n,1) * z(:,i)';
  tmp = tmp .* tmp;
  K = K + tmp;
end;

if d == 2
  mask = K < 1e-10; % to avoid singularity.
  K = 0.5 * K .* log(K + mask) .* (K>1e-10);
else
  K = - sqrt(K);
end;


