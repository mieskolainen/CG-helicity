% Wigner 3j symbol
%
% (j1 j2 j3)
% (m1 m1 m3)
%
% http://mathworld.wolfram.com/Wigner3j-Symbol.html
%
% mikael.mieskolainen@cern.ch, 2017

function wc = W3j(j1,j2,j3,m1,m2,m3)

% Check selection rules
if (W3jSR(j1,j2,j3,m1,m2,m3) == false)
    return;
end

% Create coefficient structure
c1 = -j2 + m1 + j3;
c2 = -j1 - m2 + j3;
c3 =  j1 + j2 - j3;
c4 =  j1 - m1;
c5 =  j2 + m2;

% Boundaries for which factorials are non-negative
lower = max([0,  -c1, -c2]);
upper = min([c3,  c4,  c5]);
fsum = 0;

% Evaluate sum term
for k = lower:upper
    fsum = fsum + (-1)^k / ...
        ( factorial(k) ...
        * factorial(c1+k) ...
        * factorial(c2+k) ...
        * factorial(c3-k) ...
        * factorial(c4-k) ...
        * factorial(c5-k) );
end

% Evaluate by Racah formula
wc = (-1)^(j1-j2-m3) ...
        * sqrt( triangle(j1,j2,j3) ...
        * factorial(j1+m1) * factorial(j1-m1) ...
        * factorial(j2+m2) * factorial(j2-m2) ...
        * factorial(j3+m3) * factorial(j3-m3) ) * fsum;

end

% Triangle coefficient
function d = triangle(j1,j2,j3)

d = factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) ...
    / factorial(j1+j2+j3+1);

end

% Selection rules
function state = W3jSR(j1,j2,j3,m1,m2,m3)

if ( ~cint(j1 + j2 + j3) )
    warning('W3j:: j1 + j2 + j3 not integer');
    state = false;
    return;
end
if (m1 + m2 + m3 ~= 0)
    warning('W3j:: m1 + m2 + m3 != 0');
    state = false;
    return;
end
if (j1 - m1 ~= floor(j1 - m1))
    warning('W3j:: parity of 2j_1 and 2m_1 does not match');
    state = false;
    return;
end
if (j2 - m2 ~= floor(j2 - m2))
    warning('W3j:: parity of 2j_2 and 2m_2 does not match');
    state = false;
    return;
end
if (j3 - m3 ~= floor(j3 - m3))
    warning('W3j:: parity of 2j_3 and 2m_3 does not match');
    state = false;
    return;
end
if (j1 + j2 < j3 || j3 < abs(j1 - j2))
    warning('W3j:: j over valid domain');
    state = false;
    return;
end
if (abs(m1) > j1)
    warning('W3j:: |m1| > j1');
    state = false;
    return;
end
if (abs(m2) > j2)
    warning('W3j:: |m2| > j2');
    state = false;
    return;
end
if (abs(m3) > j3)
    warning('W3j:: |m3| > j3');
    state = false;
    return;
end

state = true;
end