% Parameter vector to parity conserving helicity density matrix
%
% mikael.mieskolainen@cern.ch, 2017

function rho = xtorho(x, J)

if (J == 1)
    a = x(1);
    b = x(2);
    c = x(3);
    d = x(4);

    % Construct initial state density matrix
    rho = [0.5*(1-a), b+1i*c, d;
             b-1i*c,   a, -b+1i*c;
             d, -b-1i*c, 0.5*(1-a)];
end

if (J == 2)
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

    rho  = [a,       f+1i*c,   g+1i*d,      h+1i*e,      m;
            f-1i*c,   b,        j+1i*k,      l,          -h+1i*e;
            g-1i*d,   j-1i*k,   1-2*(a+b),  -j+1i*k,      g-1i*d;   
            h-1i*e,   l,       -j-1i*k,      b,          -f+1i*c;
            m,       -h-1i*e,   g+1i*d,     -f-1i*c,      a];
end

end