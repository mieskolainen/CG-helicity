% Get Jacob-Wick helicity amplitude weight
%
% See spindecay() for related information
%
% mikael.mieskolainen@cern.ch, 2017

function W = gethelamp(theta, phi, rho_i, s1, s2, T, theta_mother, phi_mother)

% Mother spin
J = (size(rho_i,1)-1)/2;

% Construct initial state angular momentum quantization values
% -J <= M <= J 
M_values = -J:J; % negative to positive

% Construct final state helicity configurations
lambda_values = zeros((2*s1+1)*(2*s2+1), 2);
lambda_index  = size(lambda_values);            % For indexing T matrix
i = 1;
mu = 1;
for lambda1 = -s1:s1 % negative to positive
    nu = 1;
    for lambda2 = -s2:s2 % negative to positive
        lambda_values(i,:) = [lambda1, lambda2];
        lambda_index(i,:)  = [mu,nu];
        i = i + 1;
        nu = nu + 1;
    end
    mu = mu + 1;
end

% Construct transition amplitude matrix
f = zeros((2*s1+1)*(2*s2+1), 2*J+1);

% Rows = final state spin projections
for i = 1:size(f,1)

    lambda1 = lambda_values(i,1);
    lambda2 = lambda_values(i,2);

    % Construct total lambda
    lambda = lambda1 - lambda2;

    % Columns = initial state spin projections
    for j = 1:size(f,2)

        % Pick mother helicity
        M_ = M_values(j);

        % Decay amplitudes for -J <= M <= J
        f(i,j) = wignerD(theta,phi,lambda,M_,J) * T(lambda_index(i,1),lambda_index(i,2));
    end
end

% First construct the D-matrix for spin rotation
D_ = zeros(length(M_values), length(M_values));
for i = 1:size(D_,1)
    for j = 1:size(D_,2)
        D_(i,j) = wignerD(theta_mother, phi_mother, M_values(i), M_values(j), J);
    end
end

% Rotate the initial state density matrix in the mother helicity frame
% to the pp-frame -> Mixing of pure states
rho_i_rot = D_*rho_i*D_';

% Weight of the event by the density matrix formalism
W = real( trace(f*rho_i_rot*f') );
% real() due to matlab floating point numerics

end
