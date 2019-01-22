% Spin simulation 3
%
% mikael.mieskolainen@cern.ch, 2017

clear;
close all;
addpath ../../#matlabcodes/
addpath ./src/


%% Final state quantum numbers

s1 = 0;   % Daughter 1 spin
s2 = 0;   % Daughter 2 spin
T  = 1;   % Helicity coupling matrix (2*s1+1) x (2*s1+1) (user sets this up)


%% Spin-1 or Spin-2 random parity conserving

% MC samples
N = 1e4;

% Plot different random J = 1, parity conserving density matrices
K = 50; % Number of random density matrices

% Choose spin (1 or 2)
for J = 1:2

    statistics = 1;

    % Legend texts
    legs = cell(2*J+1, 1);
    xedge = linspace(-1,1,75);   % Bin edges
    yedge = linspace(0,2*pi,75);
    
    for kk = 1:K
        
        fprintf('Matrix %d/%d \n', kk, K);
        
        % Get random parity conserving spin density matrix
        rho = randpcrho(J,0);
        
        % MC the distribution, histogram, plot
        [X,weights] = spindecay(rho, s1, s2, T, N, 0, 0);
        H = hist2w(X,weights,xedge,yedge);
        plotrho(xedge, yedge, H, rho);
        
        % Save
        filename = sprintf('./J%dfigs/%d.pdf', J, kk);
        cmd = sprintf('print -dpdf %s', filename);
        eval(cmd); pause(2);
        system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s',filename,filename));
        close all;

    end % Random matrix loop
end % J-loop

