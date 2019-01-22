% Generate spin-decay distributions for two-body decays J -> s1 + s2 
%
% input:       rho_i = Initial state (mother) spin density matrix (2J+1 x 2J+1)
%              s1    = Daughter 1 spin
%              s2    = Daughter 2 spin
%              T     = Helicity coupling matrix (2*s1+1)*(2*s2+1)
%              N     = Number of weighted events
%       theta_mother = Mother particle theta angle (set 0, if no rotation)
%       phi_mother   = Mother particle phi angle   (set 0, if no rotation)
%
% output:          X = [cos(theta),phi] matrix (N x 2) for each event
%            weights = weights for each event (N x 1)
%
% See e.g. 
% C. Amsler, J.C. Bizot, Simulation of angular distributions, 1981
%
% N.B. Function has been tested only with T = 1,
% which is the case for particular "magic" decays for which both GG's 
% <.><.> are 1, and there is only one overall \alpha_{ls} parameter
% (see below).
%
% The particular quantum number combinations, and their level of magic,
% can be found using the table created with createJW.m.
%
% ------------------------------------------------------------------------
% It is for the user to construct T matrix of size (2*s1+1)*(2*s2+1):
%
% T_{\lambda_1,\lambda_2}
%   = \sum_{ls} \alpha_{ls} x 
%          < J\lambda|ls0\lambda> <s\lambda|s1s2\lambda1,-\lambda2>  ,
%
% where two <.> terms are Glebsch-Gordans and \alpha_{ls}
% are free (data-driven/model-driven) coupling parameters.
%
% - See code below how helicities \lambda_1, \lambda_2 run in matrices.
%
% - See createJW.m for Glebsch-Gordan couplings.
%
% - We use convention where helicities run form (positive to negative)
%   such as 1,0,-1
%
% mikael.mieskolainen@cern.ch, 2017


function [X,weights] = spindecay(rho_i, s1, s2, T, N, theta_mother, phi_mother)

% Check proper normalization tr(rho) = 1
if ( abs( trace(rho_i) - 1) > 1e-6 )
   warning('Trace(rho) = 1 does not hold!');
   return;
end

%{
% Test eigenvalues
lambda0 = eigs(rho_i);

% Test that we have a valid density matrix
if (sum(lambda0 < 0) > 0)
    fprintf('Density matrix has negative eigenvalues: \n');
    disp(lambda0)
    fprintf('Fix it!\n');
    return;
end
%}

% ------------------------------------------
% von Neumann entropy
S = vnentropy(rho_i);

% Mixedness measure
mixedness = trace(rho_i^2);
fprintf('spindecay:: Tr[rho] = %0.2f, Mixedness Tr[rho^2] = %0.2f, von Neumann entropy S = %0.2f \n', ...
    trace(rho_i), mixedness, S);
% ------------------------------------------

weights = zeros(N,1);
X = zeros(N,2);

% Fully parallel event generation loop
parfor n = 1:N
    
    % Pick random decay daughter theta and phi
    costheta = -1 + 2*rand(1); % Flat in costheta
    theta = acos(costheta);    % 
    phi = 2*pi*rand(1);        % Flat in phi
    
    % Get event weight by helicity amplitude formalism
    weights(n) = gethelamp(theta, phi, rho_i, s1, s2, T, theta_mother, phi_mother);
    
    % Save event data
    X(n,:) = [costheta, phi];
end

end
