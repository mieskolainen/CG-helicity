% Multipole parameters t_LM^* (note ^*) to density matrix rho_{lambda lambda_'}
% 
% See rho2multipole
%
% input:   T_LM  =  Matrix of multipole parameters (n x 3), [L M t_LM^*]
% 
%
% For reference see: "Spin in Particle Physics", E. Leader
%
% mikael.mieskolainen@cern.ch, 2017

function rho = multipole2rho(T_LM)

% Calculate spin dimension
J = real(T_LM(end,1) / 2); % 0 <= L <= 2s
rho = zeros(2*J+1);
lambda_values = -J:1:J;

s = J;
for i = 1:2*J+1
    for j = 1:2*J+1
        lambda       = lambda_values(i);
        lambda_prime = lambda_values(j);
        summ   = 0;
        for k = 1:size(T_LM,1)
           L = T_LM(k,1);
           M = T_LM(k,2);
           t = T_LM(k,3);
           summ = summ + (2*L+1)*clebschgordan(s,L,lambda_prime,M,   s,lambda) * t;
        end
        rho(i,j) = summ;
    end
end
% Normalize by the number of states
rho = rho / (2*s + 1);

end