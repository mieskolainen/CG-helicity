% Create table of spin-couplings using Jacob-Wick helicity formalism
%
% See for example:
%
% Jacob, M., and Gr C. Wick.
% "On the general theory of collisions for particles with spin."
% Annals of Physics 7.4 (1959): 404-428.
%
%
% mikael.mieskolainen@cern.ch, 2017

addpath ./src
clear; close all;

% Basics:
%
% Intrinsic parity of fermions is + (even), and - (odd) for anti-fermions,
% this is by convention. Bosons and anti-bosons have the same parity.
% This is due to antisymmetric (fermions) versus symmetric wavefunctions
% (bosons). Parity for a system of particles is given
% by P_1 x P_2 x ... x P_N x (-1)^L
%
% Spectroscopic notation (2S + 1)L_J
% - S is the total spin quantum number,
% - L is the orbital angular momentum quantum number,
% - J is the total angular momentum quantum number
% 
% Example: parapositronium:  1S_0 (anti-parallel spins of electron/positron)
%          orthopositronium: 3S_1 (parallel spins of electron/positron)
% 
% Orbital angular momentum notation: L: S = 0, P = 1, D = 2, F = 3
% 
% Charge conjugation changes the sign of all additive quantum numbers.
% 
% C-Parity for a system of free particles is the product of C-parities,
% and C-parity is defined only for neutral particles.
% 
% C-Parity of a bound state of bosons is
% (-1)^L
% 
% C-Parity of a bound state of fermions is
% (-1)^L(-1)^(S+1)(-1) = (-1)^(L+S)
%
% \vec{J} = \vec{L} + \vec{S}

spins = [0 1/2 1 3/2 2];

qn = [];
parity = [];

%% Generate different combinations

% This could/should be made more algebraic/compact..
for i = 1:length(spins)
    for j = 1:length(spins)
        for k = 1:length(spins)
            if (j>=k) % Remove redundancy of final states (upper triangle)
                if (spins(j) == spins(k)) % Same spin for final states, remove redundancy
                    % 000
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 -1 -1];
                    % 001
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 -1 1];
                    % 010
                    %qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    %parity(end+1,:) = [-1 1 -1];
                    % 011
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 1 1];
                    % 100
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 -1 -1];
                    % 101
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 -1 1];
                    % 110
                    %qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    %parity(end+1,:) = [1 1 -1];
                    % 111
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 1 1];
                else % No redundancy here
                    % 000
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 -1 -1];
                    % 001
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 -1 1];
                    % 010
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 1 -1];
                    % 011
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [-1 1 1];
                    % 100
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 -1 -1];
                    % 101
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 -1 1];
                    % 110
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 1 -1];
                    % 111
                    qn(end+1,:) = [spins(i) spins(j) spins(k)];
                    parity(end+1,:) = [1 1 1];
                end
            end
        end
    end
end


%% Now calculate all possible A->1+2 decays

T = zeros(size(qn,1),1);

fileID = fopen('QNtable.txt','w');

for i = 1:size(qn,1)

fprintf('Decay %d/%d \n', i, size(qn,1));

fprintf(fileID, '\nDecay %d/%d \t\t\t\t \n', i, size(qn,1));
fprintf(fileID, '[%0.1f^(%s) -> %0.1f^(%s) %0.1f^(%s)]\n', ...
         qn(i,1), pm2string(parity(i,1)), qn(i,2), pm2string(parity(i,2)), qn(i,3), pm2string(parity(i,3)));
% Initial and final state tot.ang.momentum
J  = qn(i,1);
s1 = qn(i,2);
s2 = qn(i,3);

values1 = [];
values2 = [];

% Vector relations
% \vec{s} = \vec{s1} + \vec{s2}   % Jacob & Wick
% \vec{l} = \vec{J}  - \vec{s};   % By definition
warning off;

statistics = 0.5;   % both Fermions & Bosons

for s = 0:statistics:s1+s2
    for l = 0:statistics:J+s
        for lambda1 = -s1:statistics:s1
            for lambda2 = -s2:statistics:s2
                
                lambda = lambda1 - lambda2;
                
                % 1. Construct Glebsch-Gordans
                % clebschgordan(j1,j2,m1,m2,J,M)
                cg1 = clebschgordan(l,s,0,lambda,J,lambda);
                cg2 = clebschgordan(s1,s2,lambda1,-lambda2,s,lambda);
                
                % 2. Check parity conservation
                % Parity of two meson system P_tot = P_1*P_2*(-1)^L
                P_tot = parity(i,2)*parity(i,3)*(-1)^l;
                
                if (parity(i,1) == P_tot)
                    parity_conservation = 'OK';
                else
                    parity_conservation = 'FALSE';
                end
                
                % 3. Angular Momentum Conservation of z-axis quantities
                
                if (abs(lambda1 - lambda2) <= J)
                    J_conservation = 'OK'; 
                else
                    J_conservation = 'FALSE';
                end
                
                if (cg1*cg2 ~= 0 && strncmp(J_conservation, 'OK', 1))
                  fprintf(fileID, 'l = %0.1f, s = %0.1f : T = <cg1> x <cg2> = %0.3f x %0.3f \t: lambda = %0.1f, lambda1 = %0.1f lambda2 = %0.1f,    \tP = %s, \tJ = %s \n', ...
                  l, s, cg1, cg2, lambda,lambda1,lambda2, parity_conservation, J_conservation);  
                end
                
                values1(end+1) = cg1;
                values2(end+1) = cg2;
            end
        end
    end
end

T(i) = sum(values1 .* values2);

end
fclose(fileID);
