% convert quaternion to rotation matrix
% Note:  [0,0,0,1] --> eye(3);
%
% See also: transform_by_rigid3d
function [R, g] = quaternion2rotation(q)
%%=====================================================================
%% $RCSfile: quaternion2rotation.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================


x = q(1); y = q(2); z=q(3); r = q(4);
x2 = q(1) * q(1);
y2 = q(2) * q(2);
z2 = q(3) * q(3);
r2 = q(4) * q(4);

R(1,1) = r2 + x2 - y2 - z2;		% fill diagonal terms
R(2,2) = r2 - x2 + y2 - z2;
R(3,3) = r2 - x2 - y2 + z2;

xy = q(1) * q(2);
yz = q(2) * q(3);
zx = q(3) * q(1);
rx = q(4) * q(1);
ry = q(4) * q(2);
rz = q(4) * q(3);

R(1,2) = 2 * (xy + rz);			% fill off diagonal terms
R(1,3) = 2 * (zx - ry);
R(2,3) = 2 * (yz + rx);
R(2,1) = 2 * (xy - rz);
R(3,1) = 2 * (zx + ry);
R(3,2) = 2 * (yz - rx);

ss = (x2+y2+z2+r2);
R = R/ss;

if (nargout>1)
    ssss = ss*ss;
    % derivative of R(1,1) = r2 + x2 - y2 - z2;
    g1(1,1) = 4*x*(y2+z2)/ssss; g2(1,1) = -4*y*(x2+r2)/ssss;
    g3(1,1) = -4*z*(x2+r2)/ssss; g4(1,1) = 4*r*(y2+z2)/ssss;
    % derivative of R(2,2) = r2 - x2 + y2 - z2; 
    g1(2,2) = -4*x*(y2+r2)/ssss; g2(2,2) = 4*y*(x2+z2)/ssss;
    g3(2,2) = -4*z*(y2+r2)/ssss; g4(2,2) = 4*r*(x2+z2)/ssss;
    % derivative of R(3,3) = r2 - x2 - y2 + z2;
    g1(3,3) = -4*x*(z2+r2)/ssss; g2(3,3) = -4*y*(r2+z2)/ssss;
    g3(3,3) = 4*z*(x2+y2)/ssss; g4(3,3) = 4*r*(x2+y2)/ssss;

    % fill off diagonal terms
    % derivative of R(1,2) = 2 * (xy + rz);			
    g1(1,2) = 2*y/ss - 2*x*R(1,2)/ssss;
    g2(1,2) = 2*x/ss - 2*y*R(1,2)/ssss;
    g3(1,2) = 2*r/ss - 2*z*R(1,2)/ssss;
    g4(1,2) = 2*z/ss - 2*r*R(1,2)/ssss;
    % derivative of R(1,3) = 2 * (zx - ry);
    g1(1,3) = 2*z/ss - 2*x*R(1,3)/ssss;
    g2(1,3) = -2*r/ss - 2*y*R(1,3)/ssss;
    g3(1,3) = 2*x/ss - 2*z*R(1,3)/ssss;
    g4(1,3) = -2*y/ss - 2*r*R(1,3)/ssss;
    % derivative of R(2,3) = 2 * (yz + rx);
    g1(2,3) = 2*r/ss - 2*x*R(2,3)/ssss;
    g2(2,3) = 2*z/ss - 2*y*R(2,3)/ssss;
    g3(2,3) = 2*y/ss - 2*z*R(2,3)/ssss;
    g4(2,3) = 2*x/ss - 2*r*R(2,3)/ssss;
    % derivative of R(2,1) = 2 * (xy - rz);
    g1(2,1) = 2*y/ss - 2*x*R(2,1)/ssss;
    g2(2,1) = 2*x/ss - 2*y*R(2,1)/ssss;
    g3(2,1) = -2*r/ss - 2*z*R(2,1)/ssss;
    g4(2,1) = -2*z/ss - 2*r*R(2,1)/ssss;
    % derivative of R(3,1) = 2 * (zx + ry);
    g1(3,1) = 2*z/ss - 2*x*R(3,1)/ssss;
    g2(3,1) = 2*r/ss - 2*y*R(3,1)/ssss;
    g3(3,1) = 2*x/ss - 2*z*R(3,1)/ssss;
    g4(3,1) = 2*y/ss - 2*r*R(3,1)/ssss;
    % derivative of R(3,2) = 2 * (yz - rx);
    g1(3,2) = -2*r/ss - 2*x*R(3,2)/ssss;
    g2(3,2) = 2*z/ss - 2*y*R(3,2)/ssss;
    g3(3,2) = 2*y/ss - 2*z*R(3,2)/ssss;
    g4(3,2) = -2*x/ss - 2*r*R(3,2)/ssss;
    
    g{1} = g1; g{2} = g2; g{3} = g3; g{4} = g4;
end