% Check if integer or half-integer
%
% mikael.mieskolainen@cern.ch, 2017

function y = chint(x)

x = x(:);
for i = 1:length(x)
    if ( 2*x(i) ~= floor(2*x(i)) )
        y = false;
        return
    end
end
y = true;
end