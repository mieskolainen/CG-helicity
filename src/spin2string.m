% Input:    J = 0,1/2,1,3/2,2,
% Output: str = corresponding string
%
% mikael.mieskolainen@cern.ch, 2018

function str = spin2string(J)

[n,d] = rat(J);

if (d ~= 1)
    str = sprintf('%d/%d', n,d);
else
    str = sprintf('%d', n);
end

end