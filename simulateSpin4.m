% Density matrix simulation as a function of Spin Dimension
%
% mikael.mieskolainen@cern.ch, 2017

clear;
close all;
addpath ../../#matlabcodes/
addpath ./src/

measure = 0;    % Random measure (0 or 1)
N_max   = 12;   % Maximum N

NMC = 1e3;      % Number of random matrices per N
S_values  = zeros(NMC,N_max);
MX_values = zeros(NMC,N_max); 
NN_values = zeros(NMC,N_max);

tic;
for N = 1:N_max
    fprintf('N = %d/%d \n', N, N_max);
    parfor i = 1:NMC
        
        % Get random density matrix
        rho = randomrho(N,measure,0);
        
        % Get von Neumann entropy
        S_values(i,N)  = vnentropy(rho);
        
        % Mixedness
        MX_values(i,N) = trace(rho^2);
    end
end
toc;

%%
% Scaling statistics
close all;

figure;
plot(1:N_max, mean(S_values),'s-'); hold on;
plot(1:N_max, mean(MX_values), 's-');
xlabel('$N$','interpreter','latex');
%set(gca,'yscale','log');
set(gca,'xtick',1:2:N_max);
axis square; axis tight;
axis([1 N_max 0 1.8]);
l = legend('$\langle S \rangle$','$\langle$Tr$(\rho^2)\rangle$');
set(l,'interpreter','latex','location','northeast');

figure;
plot(1:N_max, var(S_values),'s-'); hold on;
plot(1:N_max, var(MX_values),'s-');
xlabel('$N$','interpreter','latex');
set(gca,'yscale','log');
set(gca,'xtick',1:2:N_max); axis square; axis tight;
l = legend('$\langle S^2 \rangle - \langle S \rangle^2$', ...
           '$\langle$Tr$(\rho^2)^2\rangle$ - $\langle$Tr$(\rho^2)\rangle^2$');
set(l,'interpreter','latex','location','northeast');


%% 1D Statistics
close all;

N_stop = min(12,N_max);
bins = 50;
legs = {};
for N = 2:N_stop
    legs{N-1} = sprintf('$N = %d$', N); 
end

fig1 = figure;
for N = 2:N_stop
    xedge  = linspace(min(S_values(:,N)), max(S_values(:,N)), bins);
    % 1d histogram
    counts = hist1m(S_values(:,N), xedge);
    delta  = xedge(2)-xedge(1);
    stephistedge(xedge, counts / sum(counts) / delta); hold on;
end
xlabel('$S$','interpreter','latex');
ylabel('$f(S)$','interpreter','latex');

l = legend(legs);
set(l,'interpreter','latex','location','northwest');
legend('boxoff'); axis square;

fig2 = figure;
for N = 2:N_stop
    xedge  = linspace(min(MX_values(:,N)), max(MX_values(:,N)), bins);
    % 1d histogram
    counts = hist1m(MX_values(:,N), xedge);
    delta  = xedge(2)-xedge(1);
    stephistedge(xedge, counts / sum(counts) / delta); hold on;
end
xlabel('Tr($\rho^2$)','interpreter','latex');
ylabel('$f$(Tr($\rho^2$))','interpreter','latex');

l = legend(legs);
set(l,'interpreter','latex','location','northeast');
legend('boxoff'); axis square;


%% 2D statistics (mixedness vs von Neumann entropy)
close all;

N_stop = min(8,N_max);

for N = 1:N_stop
    figure;
    % x and y edges
    xedge = linspace(min(MX_values(:,N)), max(MX_values(:,N)), 100);
    yedge = linspace(min(S_values(:,N)), max(S_values(:,N)), 100);
    
    % 2d histogram
    counts = hist2m([MX_values(:,N) S_values(:,N)], xedge, yedge);
    imagesc((xedge(1:end-1)+xedge(2:end))/2, ...
            (yedge(1:end-1)+yedge(2:end))/2, counts / sum(counts(:)));
    
    xlabel('Tr($\rho^2$)','interpreter','latex');
    ylabel('$S$','interpreter','latex');
    
    title(sprintf('$N = %d$', N),'interpreter','latex');
    %axis([min(MX_values(:,J)) max(MX_values(:,J)) min(S_values(:,J)) max(S_values(:,J))]);
    axis square;
    colorbar; colormap('hot'); 
    set(gca,'ydir','normal');
end

