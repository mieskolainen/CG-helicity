% Spherical Harmonics  Y_lm (theta,phi)
%
% Normalization by the usual Quantum Mechanics conventions
%
% mikael.mieskolainen@cern.ch, 2017

function Y = sphericalharmonics(theta,phi,l,m)

% Get the associated legendre polynomials
P = legendre(l, cos(theta));

Y = (-1)^m*sqrt( (2*l+1)/(4*pi)*(factor(l-m)/factor(l+m)) ) ...
     * P(m) * exp(1i*m*phi);

end