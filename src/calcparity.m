% Calculate parity matching elements of a spin density matrix
%
% input:   J = spin (1,2,3,...) boson
% output:  parity matching elements
%
% mikael.mieskolainen@cern.ch, 2017

function calcparity(J)

% (-1)^(l-l')
rho_ = @(l,lp) (-1)^(l-lp);

fprintf('\n');
for l = J:-1:-J
    for lp = J:-1:-J
        fprintf('rho_%d%d\'' \t= %d x rho_%d%d\'' \n', l,lp, rho_(l, lp), -l, -lp);
    end
end

end