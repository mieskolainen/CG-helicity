% Random von Neumann density matrix generator, with positive condition
% automatically taken into account
%
% input:       N = Dimension
%        measure = 0 for Haar/Hilbert-Schmidt, 1 for Bures measure density
%          field = 0 for complex and 1 for purely real entries
%
% output:    rho = Density matrix of size NxN
%
% mikael.mieskolainen@cern.ch, 2017

function rho = randomrho(N, measure, field)

% State matrix
r = 0;
if (field == 0) % complex
    r = randn(N) + 1i*randn(N);
elseif (field == 1) % real
    r = randn(N);
end

% Bures measure density by map r -> (M + I)r
if (measure == 1)
    r = (randomumat(N,field) + eye(N))*r;
end

% Outer product (guarantees positivity of rho) and normalize Tr = 1
rho = r*r';
rho = rho / trace(rho);

%{
% Check positivity
lambda = eigs(rho);
if (sum(lambda < 0) > 0)
   warning('randomrho:: Some eigenvalues are not >= 0!\n');
   disp(lambda(:)');
end
%}

end