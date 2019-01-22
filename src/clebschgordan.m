% SU(2) Glebsch-Gordan coefficients for the angular momentum couplings
%
% <j1,j2,m1,m2 | (j1,j2,) j,m>
%
% j1,j2,m1,m2 in the uncoupled basis
% j,m in the coupled basis
% 
% |j1,j2,j,m> = |j1,j2,m1,m2> < j1,j2,m1,m2 | j1,j2,j,m > 
%             = \sum_{m1,m2} C(j1,j2,m1,m2,j,m) |j1,j2,m1,m2>
%
% http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html
%
%
% Output:   cg  = Value
%          num  = Value^2 (Numerator)
%          den  = Value^2 (Denominator)
%
% mikael.mieskolainen@cern.ch, 2017

function [cg,num,den] = clebschgordan(j1,j2,m1,m2, j,m)

cg  = 0;
num = 0;
den = 0;

% Check selection rules
if (ggSR(j1,j2,m1,m2,j,m) == false)
    return;
end

% Use Wigner-3j symbol to evaluate
cg = (-1)^(m+j1-j2)*sqrt(2*j+1) * W3j(j1,j2,j,m1,m2,-m); % !! Note minus on m

% Rational representation of the square
[num,den] = rat(cg^2);

end

% Selection rules
function state = ggSR(j1,j2,m1,m2,j,m)

if ( ~chint([j1,j2,j]) )
    warning('glebschgordan:: (j1,j2,j) not integer or half-integer');
    state = false;
    return;
end
if ( ~cint(j1 + j2 + j) ) % "Integer perimeter rule"
    warning('glebschgordan:: j1 + j2 + j not integer');
    state = false;
    return;
end
if ( ~chint([m1,m2,m]) )
    warning('glebschgordan:: (m1,m2,m) not integer or half-integer');
    state = false;
    return;
end
if (m1 + m2 ~= m)
    warning('glebschgordan:: m1 + m2 ~= m');
    state = false;
    return;
end
if (j1 - m1 ~= floor(j1 - m1))
    warning('glebschgordan:: parity of 2j_1 and 2m_1 does not match');
    state = false;
    return;
end
if (j2 - m2 ~= floor(j2 - m2))
    warning('glebschgordan:: parity of 2j_2 and 2m_2 does not match');
    state = false;
    return;
end
if (j - m ~= floor(j - m))
    warning('glebschgordan:: parity of 2j and 2m does not match');
    state = false;
    return;
end
if (j1 + j2 < j || j < abs(j1 - j2))
    warning('glebschgordan:: j over valid domain');
    state = false;
    return;
end
if (abs(m1) > j1)
    warning('glebschgordan:: |m1| > j1');
    state = false;
    return;
end
if (abs(m2) > j2)
    warning('glebschgordan:: |m2| > j2');
    state = false;
    return;
end
if (abs(m) > j)
    warning('glebschgordan:: |m| > j');
    state = false;
    return;
end
state = true;
end