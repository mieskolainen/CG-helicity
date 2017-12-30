% Random unitary (orthogonal) rotation matrix generator by QR decomposition
%
% input:      N = dimension of the matrix
%         field = 0 for complex, 1 for real
%
% output:     M = Random unitary matrix size NxN
%
% mikael.mieskolainen@cern.ch, 2017

function M = randomumat(N, field)

r = 0;
if (field == 0)     % Complex -> Unitary
    r = randn(N) + 1i*randn(N);
elseif (field == 1) % Real -> Orthogonal
    r = randn(N);
end

% QR decompose
[Q,R] = qr(r);
d = diag(R);
d = d ./ abs(d);
M = Q*diag(d);

end