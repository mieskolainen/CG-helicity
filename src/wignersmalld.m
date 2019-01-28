% Wigner small-d matrix elements,
%
% Real valued representation (z-y-z rotation sequence)
% other choises give half of the elements imaginary.
%
% d_{mm'}^J(theta) for
% - spin projections m, m'
% - total angular momentum J
%
% https://en.wikipedia.org/wiki/Wigner_D-matrix
%
% mikael.mieskolainen@cern.ch, 2017

function element = wignersmalld(theta, m, mp, J)

factor = sqrt(factorial(J+m)*factorial(J-m)*factorial(J+mp)*factorial(J-mp));

s = 0;
for k = 0:2*J
    
    % Ignore terms with negative factorials
    if (J-mp-k >= 0 && J+m-k >= 0 && k+mp-m >= 0 && k >= 0)
    s = s + ((-1)^k/(factorial(J-mp-k)*factorial(J+m-k)*factorial(k+mp-m)*factorial(k)) )...
           *(cos(theta/2))^(2*J+m-mp-2*k)*(-sin(theta/2))^(mp-m+2*k);
    end
end
element = factor*s;

end