% von Neumann entropy of a density matrix
%
% S = -tr(rho ln(rho) ), where ln is a matrix logarithm,
%       also in terms of eigenvalues, S = -sum_j lambda_j ln(lambda_j)
%
% Input:   rho  =  Density matrix
% Output     S  =  Entropy
%
% mikael.mieskolainen@cern.ch, 2017

function S = vnentropy(rho)

% Calculate eigenvalues
lambda = eigs(rho);

% Check all positive
if (sum(lambda+1e-12 < 0) > 0)
   warning('vnentropy:: Eigenvalues not all positive real!');
   disp(lambda(:)');
end

% Equivalent to below
%S = -trace(rho * logm(rho));

% Check for pure state
% (pure states have S = 0 -> Schmidt/SVD rank 1 <-> only one eigenvalue == 1)
if (abs(sum(lambda) - 1) < 1e-12 && sum(abs(lambda - 1) < 1e-12) == 1)
    S = 0;
else
    S = -sum(lambda .* log(lambda));
end

end