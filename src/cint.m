% Check if integer
%
% mikael.mieskolainen@cern.ch, 2017

function y = cint(x)

x = x(:);
for i = 1:length(x)
    if ( x(i) ~= floor(x(i)) )
        y = false;
        return
    end
end
y = true;
end