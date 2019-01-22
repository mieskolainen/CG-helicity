% Return +,- string for parity
%
% mikael.mieskolainen@cern.ch, 2017

function out = pm2string(value)

if (value > 0)
    out = '+';
else
    out = '-';
end

end