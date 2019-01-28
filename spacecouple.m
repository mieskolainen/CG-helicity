% Couple vector spaces
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

% Upper proton-propagator-proton
% --------
%    \
%     \
J_a  = 0.5;
s1_a = 0.5;
s2_a = 2;

fa   = sym('fa', [(2*s1_a+1)*(2*J_a+1) (2*s2_a+1)]); % \lambda \lambda \M
rhoa = sym('rhoa', 2*J_a+1);

% Lower proton-propagator-proton
%     /
%    /
% --------
J_b  = 0.5;
s1_b = 0.5;
s2_b = 2;

fb   = sym('fb', [(2*s1_b+1)*(2*J_b+1) (2*s2_b+1)]); % \lambda \lambda \M
rhob = sym('rhob', 2*J_b+1);

% Propagator-Propagator-Central
% \
%  -------
% /
J_c  = 2;
s1_c = s2_a;
s2_c = s2_b;

fc = sym('fc', [(2*J_c+1) (2*s1_c+1)*(2*s2_c+1)]); % 

%
% Central-Daughter1-Daughter2
%       /
% ------
%       \
J_d  = J_c;
s1_d = 0;
s2_d = 0;

fd = sym('fd', [(2*s1_d+1)*(2*s2_d+1) 2*J_d+1]);
%}

% Initial state density matrix
rho_i = kron(rhoa, rhob);

ftot = kron(fa, fb) * fc'





