% Create table of relative decay helicity amplitude couplings in the
% Jacob-Wick helicity formalism
%
% See:
% Create table of spin-couplings using Jacob-Wick helicity formalism
%
% See for example:
%
% Jacob, M., and Gr C. Wick.
% "On the general theory of collisions for particles with spin."
% Annals of Physics 7.4 (1959): 404-428.
%
% and
%
% Amsler, Bizot, Simulation of Angular Distributions, Computer Physics
% Communications, (1981)
%
% In addition to the CG-couplings, one needs (free) parameters \alpha_ls
% (theory/experiment)
% in the helicity matrix elements T_ij. These give the weighted expansion:
% 
% T_{\lambda_1, \lambda_2}
% = \sum_{ls} \alpha_{ls} <J\lambda|ls 0 \lambda > <s\lambda | s_1 s_2 \lambda_1, -\lambda_2>,
% := \sum_{ls} \alpha_{ls} <Clebsch-Gordan 1> <Clebsch-Gordan 2>
%
% mikael.mieskolainen@cern.ch, 2018

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

%spins = [0 1/2 1 3/2 2];
spins = [0 1 2];

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


%% Now calculate all possible A->1+2 decay helicity amplitudes

T    = zeros(size(qn,1),1);
T_PC = zeros(size(qn,1),1);

fileID = fopen('QNtable.txt','w');

fprintf(fileID, 'Mikael Mieskolainen, 2018 \n');
fprintf(fileID, '\n');
fprintf(fileID, 'Helicity decay matrix T (2s_1 + 1)x(2s_2 + 1) elements defined as: \n\n');
fprintf(fileID, 'T_{\\lambda_1, \\lambda_2} \n');
fprintf(fileID, '= \\sum_{ls} \\alpha_{ls} <J \\lambda | l s 0 \\lambda> <s \\lambda | s_1 s_2 \\lambda_1 -\\lambda_2>\n');
fprintf(fileID, ':= \\sum_{ls} \\alpha_{ls} <clebsch-gordan1> <clebsch-gordan2>,');
fprintf(fileID, '\n');
fprintf(fileID, 'where \\alpha_{ls} are free parameters (from theory/experiment).\n');
fprintf(fileID, '\n');
fprintf(fileID, 'Helicity is \\lambda = \\vec{J}.\\vec{p}/|p| = \\vec{l}.\\vec{p}/|p| + m_s = m_s, \n');
fprintf(fileID, 'since (\\vec{l}.\\vec{p} = 0 always), -s <= m_s <= s \n');
fprintf(fileID, '\n');
fprintf(fileID, '\\lambda := \\lambda_1 - \\lambda_2 (Jacob-Wick) \n');
fprintf(fileID, '\\vec{s} = \\vec{s_1} + \\vec{s_2}, 0 <= |s| <= |s_1| + |s_2| \n');
fprintf(fileID, '\\vec{l} = \\vec{J} - \\vec{s}, 0 <= |l| <= |J| + |s| (by triangle ineq.) \n');


for i = 1:size(qn,1)

fprintf('%d/%d \n', i, size(qn,1));

fprintf(fileID, '\nDecay %d/%d \t\t\t\t \n', i, size(qn,1));

fprintf(fileID, '[%s^%s > %s^%s %s^%s (J^P > s_1^P_1 s_2^P_2)]\n', ...
         spin2string(qn(i,1)), pm2string(parity(i,1)), ...
         spin2string(qn(i,2)), pm2string(parity(i,2)), ...
         spin2string(qn(i,3)), pm2string(parity(i,3)));

% Initial and final state tot.ang.momentum
J  = qn(i,1);
s1 = qn(i,2);
s2 = qn(i,3);

values1 = [];
values2 = [];

values1_PC = [];
values2_PC = [];

% Vector relations
% \vec{s} = \vec{s1} + \vec{s2}   % Jacob & Wick
% \vec{l} = \vec{J}  - \vec{s};   % By definition
warning off;

%statistics = 0.5;   % both Fermions & Bosons
statistics = 1.0;    % Bosons only

for s = -(s1+s2):statistics:s1+s2
    for l = -(J+s):statistics:J+s
        for lambda1 = -s1:statistics:s1
            for lambda2 = -s2:statistics:s2
                
                % --------------------------------------------------------
                % [0. Definition]
                lambda = lambda1 - lambda2;
                
                % --------------------------------------------------------
                % [1. Construct Clebsch-Gordans]
                
                % <J lambda | ls 0 lambda>
                [cg1,cg1num,cg1den] = clebschgordan(l,s,0,lambda,  J,lambda);
                
                % <s lambda | s1 s2 lambda1 -lambda2> (note minus on
                % lambda2!)
                [cg2,cg2num,cg2den] = clebschgordan(s1,s2,lambda1,-lambda2,  s,lambda);
         
                % --------------------------------------------------------
                % 2. Check parity conservation
                P_tot = parity(i,2)*parity(i,3)*(-1)^l;
                
                if (parity(i,1) == P_tot)
                    parity_conservation = 'OK';
                else
                    parity_conservation = 'FALSE';
                end
                
                % --------------------------------------------------------
                % [3. Angular Momentum Conservation of z-axis quantities]
                %if (abs(lambda1 - lambda2) <= J)
                    J_conservation = 'OK';
                %else
                %    J_conservation = 'FALSE';
                %end
                
                if (cg1*cg2 ~= 0 && strncmp(J_conservation, 'OK', 1))

                  fprintf(fileID, 'l = %s, s = %s, \\lambda_1 = %2s \\lambda_2 = %2s : \\lambda = %2s, P = %5s, <cg1> x <cg2> = %s x %s \n', ...
                  spin2string(l), spin2string(s), spin2string(lambda1), spin2string(lambda2), spin2string(lambda), ...
                  parity_conservation, cg2string(cg1,cg1num,cg1den), cg2string(cg2,cg2num,cg2den));
                end

                % Save values
                values1(end+1) = cg1;
                values2(end+1) = cg2;
                if (strncmp(parity_conservation, 'OK', 1))
                    values1_PC(end+1) = cg1;
                    values2_PC(end+1) = cg2;
                end
            end
        end
    end
end

T(i)    = sum(values1 .* values2);
T_PC(i) = sum(values1_PC .* values2_PC);

fprintf(fileID, 'info: \\sum_{ls\\lambda1\\lambda2} <cg1> x <cg2> = %0.2f \n', T(i));
fprintf(fileID, 'info: \\sum_{ls\\lambda1\\lambda2} <cg1> x <cg2> = %0.2f (P conserving) \n', T_PC(i));

end
fclose(fileID);
