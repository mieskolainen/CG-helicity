% Active rotate in SO(3) by (theta, phi)
%
% input:     x = 3x1 vector
%        theta = radians, [0,pi]
%          phi = radians, [-pi,pi] or [0, 2pi]
%
% output:    y = 3x1 vector
%
% mikael.mieskolainen@cern.ch, 2017

function y = rot_thetaphi(x, theta, phi)

x = x(:);

cthe = cos(theta);
sthe = sin(theta);
cphi = cos(phi);
sphi = sin(phi);

R = [cthe*cphi, -sphi, sthe*cphi;
     cthe*sphi,  cphi, sthe*sphi;
    -sthe, 0, cthe];

y = R*x; 

end