% Positivity constraint for spin-2 density matrix with parity conservation
%
% rho  = [a,       f+1i*c,   g+1i*d,      h+1i*e,      m;
%         f-1i*c,   b,        j+1i*k,      l,          -h+1i*e;
%         g-1i*d,   j-1i*k,   1-2*(a+b),  -j+1i*k,      g-1i*d;   
%         h-1i*e,   l,       -j-1i*k,      b,          -f+1i*c;
%         m,       -h-1i*e,   g+1i*d,     -f-1i*c,      a];     
%        
%
% mikael.mieskolainen@cern.ch, 2017

function [c,ceq] = pconst2(x)

a = x(1);
b = x(2);
c = x(3);
d = x(4);
e = x(5);
f = x(6);
g = x(7);
h = x(8);
j = x(9);
k = x(10);
l = x(11);
m = x(12);

% Positivity conditions (from symbolic code)
pc(1) = - 6*a^2 - 8*a*b + 4*a - 6*b^2 + 4*b - 4*c^2 - 4*d^2 - 4*e^2 - 4*f^2 - 4*g^2 - 4*h^2 - 4*j^2 - 4*k^2 - 2*l^2 - 2*m^2;
pc(2) = - 12*a^3 - 48*a^2*b + 6*a^2 - 48*a*b^2 + 24*a*b + 12*a*c^2 - 12*a*d^2 + 12*a*e^2 + 12*a*f^2 - 12*a*g^2 + 12*a*h^2 - 24*a*j^2 - 24*a*k^2 + 12*a*m^2 - 12*b^3 + 6*b^2 + 12*b*c^2 - 24*b*d^2 + 12*b*e^2 + 12*b*f^2 - 24*b*g^2 + 12*b*h^2 - 12*b*j^2 - 12*b*k^2 + 12*b*l^2 - 12*c^2 + 24*c*d*j + 24*c*e*l - 24*c*e*m - 24*c*g*k + 12*d^2*m - 24*d*e*j + 24*d*f*k - 24*d*h*k - 12*e^2 + 24*e*g*k - 12*f^2 + 24*f*g*j + 24*f*h*l - 24*f*h*m + 12*g^2*m - 24*g*h*j - 12*h^2 - 12*j^2*l - 12*k^2*l - 6*l^2 - 6*m^2;
pc(3) = 48*a^2*b - 168*a^2*b^2 - 96*a^3*b + 96*a^2*c^2 + 96*a^2*e^2 + 96*a^2*f^2 + 96*a^2*h^2 - 48*a^2*j^2 - 48*a^2*k^2 + 72*a^2*l^2 - 96*a*b^3 + 48*a*b^2 + 144*a*b*c^2 - 96*a*b*d^2 + 144*a*b*e^2 + 144*a*b*f^2 - 96*a*b*g^2 + 144*a*b*h^2 - 96*a*b*j^2 - 96*a*b*k^2 + 96*a*b*l^2 + 96*a*b*m^2 - 48*a*c^2 + 96*a*c*d*j - 96*a*c*e*l + 192*a*c*e*m - 96*a*c*g*k - 96*a*d*e*j + 96*a*d*f*k - 96*a*d*h*k - 48*a*e^2 + 96*a*e*g*k - 48*a*f^2 + 96*a*f*g*j - 96*a*f*h*l + 192*a*f*h*m - 96*a*g*h*j - 48*a*h^2 - 96*a*j^2*l - 96*a*k^2*l - 48*a*l^2 + 96*b^2*c^2 - 48*b^2*d^2 + 96*b^2*e^2 + 96*b^2*f^2 - 48*b^2*g^2 + 96*b^2*h^2 + 72*b^2*m^2 - 48*b*c^2 + 96*b*c*d*j - 192*b*c*e*l + 96*b*c*e*m - 96*b*c*g*k + 96*b*d^2*m - 96*b*d*e*j + 96*b*d*f*k - 96*b*d*h*k - 48*b*e^2 + 96*b*e*g*k - 48*b*f^2 + 96*b*f*g*j - 192*b*f*h*l + 96*b*f*h*m + 96*b*g^2*m - 96*b*g*h*j - 48*b*h^2 - 48*b*m^2 + 24*c^4 + 48*c^2*d^2 - 48*c^2*e^2 + 48*c^2*f^2 + 48*c^2*g^2 + 48*c^2*h^2 + 48*c^2*j^2 + 48*c^2*k^2 + 48*c^2*l*m + 96*c*d^2*e + 96*c*d*j*l - 96*c*d*j*m - 192*c*e*f*h + 96*c*e*g^2 + 96*c*e*j^2 + 96*c*e*k^2 + 96*c*e*l - 96*c*e*m - 96*c*g*k*l + 96*c*g*k*m + 48*d^2*e^2 + 48*d^2*f^2 + 96*d^2*f*h + 48*d^2*h^2 + 48*d^2*l^2 - 96*d*e*j*l + 96*d*e*j*m + 96*d*f*k*l - 96*d*f*k*m - 96*d*h*k*l + 96*d*h*k*m + 24*e^4 + 48*e^2*f^2 + 48*e^2*g^2 + 48*e^2*h^2 + 48*e^2*j^2 + 48*e^2*k^2 + 48*e^2*l*m + 96*e*g*k*l - 96*e*g*k*m + 24*f^4 + 48*f^2*g^2 - 48*f^2*h^2 + 48*f^2*j^2 + 48*f^2*k^2 + 48*f^2*l*m + 96*f*g^2*h + 96*f*g*j*l - 96*f*g*j*m + 96*f*h*j^2 + 96*f*h*k^2 + 96*f*h*l - 96*f*h*m + 48*g^2*h^2 + 48*g^2*l^2 - 96*g*h*j*l + 96*g*h*j*m + 24*h^4 + 48*h^2*j^2 + 48*h^2*k^2 + 48*h^2*l*m + 48*j^2*m^2 + 48*k^2*m^2 + 24*l^2*m^2;
pc(4) = 120*(c^2 + 2*c*e + e^2 + f^2 + 2*f*h + h^2 - a*b - a*l + b*m + l*m)*(a*l - 2*c*e - a*b - 2*f*h - b*m + l*m + 2*a*b^2 + 2*a^2*b - 2*a*c^2 - 2*b*c^2 - 2*a*e^2 + 2*b*d^2 - 2*a*f^2 - 2*b*e^2 - 2*b*f^2 - 2*a*h^2 + 2*b*g^2 - 2*b*h^2 + 2*a*j^2 + 2*a*k^2 - 2*a^2*l + 2*b^2*m - 2*d^2*l - 2*g^2*l + 2*j^2*m + 2*k^2*m + c^2 + e^2 + f^2 + h^2 + 4*a*c*e + 4*b*c*e - 2*a*b*l + 4*a*f*h + 2*a*b*m + 4*b*f*h - 4*c*d*j + 4*d*e*j + 4*c*g*k - 4*d*f*k + 4*d*h*k - 4*e*g*k - 4*f*g*j + 4*g*h*j - 2*a*l*m - 2*b*l*m);
        
if (pc(1) < 0 || pc(2) < 0 || pc(3) < 0 || pc(4) < 0)
    c = 1e6; % Give penalty
else
    c = 0;   % No penalty
end

ceq = [];
end