% Spin simulation 2
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


%% Plot pure spin states

% Number of events
N = 1e4;

% Spin
for J = 1:2

    statistics = 1; % Bose

    fig1 = figure;
    fig2 = figure;

    % Legend texts
    legs = cell(2*J+1, 1);
    xedge   = linspace(-1,1,75);   % cos(theta)
    yedge   = linspace(0,2*pi,75); % phi
    xcenter = (xedge(1:end-1) + xedge(2:end))/2; % bin centers
    ycenter = (yedge(1:end-1) + yedge(2:end))/2; % bin centers
    
    % Simulate each pure states
    for kk = 1:2*J+1

        % Loop over diagonal [J <= M <= -J]
        rho_i = zeros(2*J+1); rho_i(kk,kk) = 1; 
        [X,weights] = spindecay(rho_i, s1, s2, T, N, 0, 0);
        
        M = J:-statistics:-J;
        H = hist2w(X,weights,xedge,yedge);
        
        % Plot different histograms
        
        % 2D histogram
        figure(fig1);
        subplot(floor(sqrt(2*J+1) + 0.5), ceil(sqrt(2*J+1)), kk);
        imagesc(xcenter, ycenter, H' / sum(H(:))); colormap(hot);% colorbar;
        xlabel('$\cos(\theta)$','interpreter','latex');
        ylabel('$\phi$ (rad)','interpreter','latex');
        title(sprintf('$|%d,%d \\rangle$', J, M(kk)),'interpreter','latex'); 
        axis square;
        
        % 1D histogram
        figure(fig2);
        stephistedge(xedge, sum(H',1) ); hold on;
        xlabel('$\cos(\theta)$','interpreter','latex');
        title(sprintf('Pure states $|J,\\lambda \\rangle$'),'interpreter','latex'); 
        legs{kk} = sprintf('$|%d,%d \\rangle$', J, M(kk));
    end

    figure(fig1);
    filename = sprintf('./figs/J%d_pure_2D.pdf', J);
    cmd = sprintf('print -dpdf %s', filename);
    eval(cmd); pause(2); close(fig1);
    system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename)); pause(2);
    
    figure(fig2);
    l = legend(legs);
    set(l,'interpreter','latex','location','southeast');
    axis tight; axis square;
    set(gca,'yscale','log');
    
    filename = sprintf('./figs/J%d_pure_1D.pdf', J);
    cmd = sprintf('print -dpdf %s', filename);
    eval(cmd); pause(2); close(fig2);
    system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename)); pause(2);
end
