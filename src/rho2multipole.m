% Helicity spin density matrix to helicity multipole t_LM^* (note ^*) parameters
% 
% Forward:
% t_LM^* = sum_{mm'} <sm | sm';LM> rho_{mm'}
% 
% Inverse:
% rho_{mm'} = 1/(2s+1) sum_{LM} (2L+1) <sm | sm'; LM> t_{LM}^*
% 
% The number L is
% 0 <= L <= 2s
% 
% The number M is bounded by
% -L <= M <= L
% 
% Output   T_LM  =  n x 3, format [L M t_LM^*]
%
% For reference see: "Spin in Particle Physics", E. Leader
%
% mikael.mieskolainen@cern.ch, 2017

function T_LM = rho2multipole(rho)

J = (size(rho,1)-1)/2;
n = 2*J+1;

% t_LM^* = \sum_{lambda, lambda'} <s lambda | s lambda'; LM> rho_{lambda,lambda'}
lambda_val       = -J:1:J;
lambda_prime_val = -J:1:J;
s = J;

k = 1;
for L = 0:2*J
    for M = -L:L
        summ = 0;
        for i = 1:size(rho,1)
            for j = 1:size(rho,2)
                lambda       = lambda_val(i);
                lambda_prime = lambda_prime_val(j);
                
                % <j1,j2,m1,m2 | (j1,j2,) j,m>
                summ = summ + clebschgordan(s,L,lambda_prime,M,   s,lambda) * rho(i,j);
            end
        end
        T_LM(k,:) = [L M summ]; k = k + 1;
    end
end
end
