% Ordinary Legendre polynomials
% 
% Input:     x = input values, e.g., -1...1 (d x 1) vector
%            N = order, P_0(x), P_1(x), ... P_N(X)
% Output:    y = output values (N+1 x d)
% 
% Via recursive definition:
% https://en.wikipedia.org/wiki/Legendre_polynomials
%
% mikael.mieskolainen@cern.ch, 2017

function P = legendrepoly(x,N)

x = x(:)'; % Make it row

P  = zeros(N+1,length(x));
P(1,:) = 1;

if (N > 0) % Order n = 1
    P(2,:) = x;
end
if (N > 1)  % Order n = 2,3,4,...
    for n = 1:N
        k = n+1; % Matlab indexing for indexing P
        P(k+1,:) = ( (2*n+1)*x.*P(k,:)-n*P(k-1,:) ) / (n+1);
    end
end

% Cut to order requested
P = P(1:N+1,:);

end