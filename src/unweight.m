% Unweight a set of events
%
% input  X  =  A set of weighted events (N x K), K = the observables
%        W  =  Event weights
%
% output Y  =  Unweighted events << N (statistics reduced, naturally)
%
% mikael.mieskolainen@cern.ch, 2017

function Y = unweight(X, W)

wmax = max(W(:));

Y = [];

% von Neumann acceptance-rejection
for i = 1:size(X,1)
    u = rand(1);
    if (u < W(i)/wmax)
        Y(end+1,:) = X(i,:);
    end
end

fprintf('Unweighting efficiency: %0.5f (%d/%d) \n', size(Y,1) / size(X,1), size(Y,1), size(X,1));

end
