% Clebsch-Gordan value^2 to a string
%
% Input:  cg        = Clebsch-Gordan value
%         num, den  = Numerator, Denominator of Clebsch-Gordan-value^2
%
% mikael.mieskolainen@cern.ch, 2018

function str = cg2string(cg, num, den)

if (sign(cg) > 0)
    sgn = '';
else
    sgn = '-';
end
if (num == 1 && den == 1)
    str = sprintf('1');
    return;
end
if (den == 1)
    str = sprintf('%s\\sqrt{%d}', sgn, num);
else
    str = sprintf('%s\\sqrt{%d/%d}', sgn, num, den);
end

end