% Generate random Parity Conserving spin density matrix for
% spin J = 1 or J = 2
%
% Input:          J = 1 or 2 (spin) 
%           measure = Random measure, 0 for Haar, 1 for Bures
% Output:       rho = Spin density matrix of size (2J+1)
%
% For analytic positivity conditions see:
% P. Minnaert, Physical Review, 1966
%
% Otherwise, one can check numerically all eigenvalues > 0
%
% mikael.mieskolainen@cern.ch, 2017

function rho = randpcrho(J, measure)

trials = 0;

while (true)

    % Get random density matrix
    rho = randomrho(2*J+1, measure, 0);
    
    % Spin 1
    if (J == 1)
    
        % 1.1 Purely real elements by parity and hermiticity
        rho(1,3) = real( rho(1,3) );
        
        % 1.2 Parity constrained
        rho(2,3) = -real( rho(1,2) ) + 1i*imag(rho(1,2) );

        % 2. Lower triangle = hermiticity constrained
        for i = 1:3
            for j = 1:i
                rho(i,j) = conj( rho(j,i) );
            end
        end

        % 3. Diagonal parity constrained
        rho(3,3) = real( rho(1,1) );
        
        %{
        rho = [0.5*(1-a), b+1i*c, d;
                 b-1i*c,   a, -b+1i*c;
                 d, -b-1i*c, 0.5*(1-a)];
        %}
        
        % Important, renormalize at this point (=trace(rho) = 1 condition)
        rho = rho / trace(rho);
        
        a = real(rho(2,2));
		b = real(rho(1,2));
		c = imag(rho(1,2));
		d = real(rho(1,3));
        
        % FINAL CHECK: Positivity conditions (from symbolic code)
        pc(1) = - 4*b^2 - 4*c^2 - 2*d^2 + a - (3*a^2)/2 + 1/2;
        pc(2) = -(3*(2*d - a + 1)*(2*a*d - a + a^2 + 4*b^2 + 4*c^2))/2;
        
        if (pc(1) > 0 && pc(2) > 0)
            return;
        end
        trials = trials + 1;
    end
    
    % Spin 2
    if (J == 2)

        % Make it parity symmetric

        % 1. First right semi-triangle

        % 1.1 Purely real
        rho(1,5) = real( rho(1,5));
        rho(2,4) = real( rho(2,4));

        % 1.2 Parity constrained
        rho(2,5) = -real( rho(1,4)) + 1i*imag(rho(1,4));
        rho(3,4) = -real( rho(2,3)) + 1i*imag(rho(2,3));
        rho(3,5) = conj(rho(1,3));
        rho(4,5) = -real( rho(1,2)) + 1i*imag(rho(1,2));

        % 2. Lower triangle = hermiticity constrained
        for i = 1:5
            for j = 1:i
                rho(i,j) = conj( rho(j,i) );
            end
        end

        % 3. Diagonal parity constrained
        rho(4,4) = rho(2,2);
        rho(5,5) = rho(1,1);
        %{
        rho  = [a,       f+1i*c,   g+1i*d,      h+1i*e,      m;
                f-1i*c,   b,        j+1i*k,      l,          -h+1i*e;
                g-1i*d,   j-1i*k,   1-2*(a+b),  -j+1i*k,      g-1i*d;   
                h-1i*e,   l,       -j-1i*k,      b,          -f+1i*c;
                m,       -h-1i*e,   g+1i*d,     -f-1i*c,      a];     
        %}

        % Important, renormalize at this point (=trace(rho) = 1 condition)
        rho = rho / trace(rho);
        
        a = rho(1,1);
        b = rho(2,2);
        c = imag(rho(1,2));
        d = imag(rho(1,3));
        e = imag(rho(1,4));
        f = real(rho(1,2));
        g = real(rho(1,3));
        h = real(rho(1,4));
        j = real(rho(2,3));
        k = imag(rho(2,3));
        l = rho(2,4);
        m = rho(1,5);
        
        % FINAL CHECK: Positivity conditions (from symbolic code)
        pc(1) = - 6*a^2 - 8*a*b + 4*a - 6*b^2 + 4*b - 4*c^2 - 4*d^2 - 4*e^2 - 4*f^2 - 4*g^2 - 4*h^2 - 4*j^2 - 4*k^2 - 2*l^2 - 2*m^2;
        pc(2) = - 12*a^3 - 48*a^2*b + 6*a^2 - 48*a*b^2 + 24*a*b + 12*a*c^2 - 12*a*d^2 + 12*a*e^2 + 12*a*f^2 - 12*a*g^2 + 12*a*h^2 - 24*a*j^2 - 24*a*k^2 + 12*a*m^2 - 12*b^3 + 6*b^2 + 12*b*c^2 - 24*b*d^2 + 12*b*e^2 + 12*b*f^2 - 24*b*g^2 + 12*b*h^2 - 12*b*j^2 - 12*b*k^2 + 12*b*l^2 - 12*c^2 + 24*c*d*j + 24*c*e*l - 24*c*e*m - 24*c*g*k + 12*d^2*m - 24*d*e*j + 24*d*f*k - 24*d*h*k - 12*e^2 + 24*e*g*k - 12*f^2 + 24*f*g*j + 24*f*h*l - 24*f*h*m + 12*g^2*m - 24*g*h*j - 12*h^2 - 12*j^2*l - 12*k^2*l - 6*l^2 - 6*m^2;
        pc(3) = 48*a^2*b - 168*a^2*b^2 - 96*a^3*b + 96*a^2*c^2 + 96*a^2*e^2 + 96*a^2*f^2 + 96*a^2*h^2 - 48*a^2*j^2 - 48*a^2*k^2 + 72*a^2*l^2 - 96*a*b^3 + 48*a*b^2 + 144*a*b*c^2 - 96*a*b*d^2 + 144*a*b*e^2 + 144*a*b*f^2 - 96*a*b*g^2 + 144*a*b*h^2 - 96*a*b*j^2 - 96*a*b*k^2 + 96*a*b*l^2 + 96*a*b*m^2 - 48*a*c^2 + 96*a*c*d*j - 96*a*c*e*l + 192*a*c*e*m - 96*a*c*g*k - 96*a*d*e*j + 96*a*d*f*k - 96*a*d*h*k - 48*a*e^2 + 96*a*e*g*k - 48*a*f^2 + 96*a*f*g*j - 96*a*f*h*l + 192*a*f*h*m - 96*a*g*h*j - 48*a*h^2 - 96*a*j^2*l - 96*a*k^2*l - 48*a*l^2 + 96*b^2*c^2 - 48*b^2*d^2 + 96*b^2*e^2 + 96*b^2*f^2 - 48*b^2*g^2 + 96*b^2*h^2 + 72*b^2*m^2 - 48*b*c^2 + 96*b*c*d*j - 192*b*c*e*l + 96*b*c*e*m - 96*b*c*g*k + 96*b*d^2*m - 96*b*d*e*j + 96*b*d*f*k - 96*b*d*h*k - 48*b*e^2 + 96*b*e*g*k - 48*b*f^2 + 96*b*f*g*j - 192*b*f*h*l + 96*b*f*h*m + 96*b*g^2*m - 96*b*g*h*j - 48*b*h^2 - 48*b*m^2 + 24*c^4 + 48*c^2*d^2 - 48*c^2*e^2 + 48*c^2*f^2 + 48*c^2*g^2 + 48*c^2*h^2 + 48*c^2*j^2 + 48*c^2*k^2 + 48*c^2*l*m + 96*c*d^2*e + 96*c*d*j*l - 96*c*d*j*m - 192*c*e*f*h + 96*c*e*g^2 + 96*c*e*j^2 + 96*c*e*k^2 + 96*c*e*l - 96*c*e*m - 96*c*g*k*l + 96*c*g*k*m + 48*d^2*e^2 + 48*d^2*f^2 + 96*d^2*f*h + 48*d^2*h^2 + 48*d^2*l^2 - 96*d*e*j*l + 96*d*e*j*m + 96*d*f*k*l - 96*d*f*k*m - 96*d*h*k*l + 96*d*h*k*m + 24*e^4 + 48*e^2*f^2 + 48*e^2*g^2 + 48*e^2*h^2 + 48*e^2*j^2 + 48*e^2*k^2 + 48*e^2*l*m + 96*e*g*k*l - 96*e*g*k*m + 24*f^4 + 48*f^2*g^2 - 48*f^2*h^2 + 48*f^2*j^2 + 48*f^2*k^2 + 48*f^2*l*m + 96*f*g^2*h + 96*f*g*j*l - 96*f*g*j*m + 96*f*h*j^2 + 96*f*h*k^2 + 96*f*h*l - 96*f*h*m + 48*g^2*h^2 + 48*g^2*l^2 - 96*g*h*j*l + 96*g*h*j*m + 24*h^4 + 48*h^2*j^2 + 48*h^2*k^2 + 48*h^2*l*m + 48*j^2*m^2 + 48*k^2*m^2 + 24*l^2*m^2;
        pc(4) = 120*(c^2 + 2*c*e + e^2 + f^2 + 2*f*h + h^2 - a*b - a*l + b*m + l*m)*(a*l - 2*c*e - a*b - 2*f*h - b*m + l*m + 2*a*b^2 + 2*a^2*b - 2*a*c^2 - 2*b*c^2 - 2*a*e^2 + 2*b*d^2 - 2*a*f^2 - 2*b*e^2 - 2*b*f^2 - 2*a*h^2 + 2*b*g^2 - 2*b*h^2 + 2*a*j^2 + 2*a*k^2 - 2*a^2*l + 2*b^2*m - 2*d^2*l - 2*g^2*l + 2*j^2*m + 2*k^2*m + c^2 + e^2 + f^2 + h^2 + 4*a*c*e + 4*b*c*e - 2*a*b*l + 4*a*f*h + 2*a*b*m + 4*b*f*h - 4*c*d*j + 4*d*e*j + 4*c*g*k - 4*d*f*k + 4*d*h*k - 4*e*g*k - 4*f*g*j + 4*g*h*j - 2*a*l*m - 2*b*l*m);
        
        if ( pc(1) > 0 && pc(2) > 0 && pc(3) > 0 && pc(4) > 0)
            return;
        end
        
        % Alternative, check all eigenvalues are real >= 0
        %lambda = eigs(rho);
        %if (sum(lambda < 0) == 0)
        %    return;
        %end
        trials = trials + 1;
    end
end

end