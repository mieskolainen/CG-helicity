% Wigner large D
%
% Wigner small d_{mm'}^J(theta) ~ rotation around y_z-axis
% Wigner large D(theta,phi) = exp(iJ_y theta) exp(iJ_z phi)
%
%
% mikael.mieskolainen@cern.ch, 2017

function value = wignerD(theta,phi,m,mp,J)

value = wignersmalld(theta,m,mp,J).*exp(1i*mp*phi);

end

