% Positivity constraint for spin-1 density matrix with parity conservation
%
% rho = [0.5*(1-a), b+1i*c, d;
%        b-1i*c,   a, -b+1i*c;
%        d, -b-1i*c, 0.5*(1-a)];
%
% mikael.mieskolainen@cern.ch, 2017

function [c,ceq] = pconst1(x)

a = x(1);
b = x(2);
c = x(3);
d = x(4);

% Positivity conditions (from symbolic code)
pc(1) = - 4*b^2 - 4*c^2 - 2*d^2 + a - (3*a^2)/2 + 1/2;
pc(2) = -(3*(2*d - a + 1)*(2*a*d - a + a^2 + 4*b^2 + 4*c^2))/2;

if (pc(1) < 0 || pc(2) < 0)
    c = 1e6; % Give penalty
else
    c = 0;   % No penalty
end

ceq = [];
end