% Visualize the spin density matrix
%
% input:  xedge = x axis bin EDGES
%         yedge = y axis bin EDGES
%             H = 2D histogram 
%           rho = spin density matrix (2J+1)x(2J+1)
%
% output:   fig = figure handle 
%
% mikael.mieskolainen@cern.ch, 2017

function fig = plotrho(xedge, yedge, H, rho)

% BIN centers for imagesc
xcenter = (xedge(1:end-1) + xedge(2:end))/2;
ycenter = (yedge(1:end-1) + yedge(2:end))/2;

% von Neumann entropy
S = vnentropy(rho);

% Mixedness measure
mixedness = trace(rho^2);
fprintf('plotrho:: Tr[rho] = %0.2f, Mixedness Tr[rho^2] = %0.2f, von Neumann entropy S = %0.2f \n', ...
    trace(rho), mixedness, S);    

J = (size(rho)-1) / 2;

%% 2D histogram
fig = figure;
subplot(1,3,1);
imagesc(xcenter, ycenter, H' / sum(H(:))); colormap(hot);
xlabel('$\cos(\theta)$','interpreter','latex');
ylabel('$\phi$ (rad)','interpreter','latex');
axis square;
set(gca,'xtick', round(linspace(-1,1,5),1));
set(gca,'ytick', [0 2 4 6]);
axis([-1.01 1.01 0 2*pi]);

%% 1D histograms

delta = (xedge(2) - xedge(1));

% Cos(Theta)
subplot(1,3,2);
sm = sum(H',1); sm = sm/sum(sm)/delta;
stephistedge(xedge, sm); hold on;
set(gca,'xtick', round(linspace(-1,1,5),1));
set(gca,'ytick', [0 0.25 0.5 0.75 1]);
axis([-1 1 0 1]);
axis square;
xlabel('$\cos(\theta)$','interpreter','latex');

if (J == 1)
    title(sprintf('$S = %0.2f$, Tr$[\\rho^2] = %0.2f : \\rho_{00} = %0.2f,$ Re$[\\rho_{10}] = %0.2f,$ Im$[\\rho_{10}] = %0.2f, \\rho_{1-1} = %0.2f$', ... 
        S, mixedness, rho(2,2), real(rho(1,2)), imag(rho(1,2)), rho(1,3) ),'interpreter','latex'); 
end
if (J == 2)
    title(sprintf('$S = %0.2f$, Tr$[\\rho^2] = %0.2f : \\rho_{22} = %0.2f, \\rho_{11} = %0.2f, \\rho_{00} = %0.2f \\, ...$', ...
        S, mixedness, rho(1,1), rho(2,2), rho(3,3)),'interpreter','latex'); 
end

subplot(1,3,3);
% Phi
sm = sum(H,1); sm = sm / sum(sm) / delta;
stephistedge(yedge, sm); hold on;
xlabel('$\phi$ (rad)','interpreter','latex');
axis square;
set(gca, 'xtick', [0 2 4 6]);
set(gca,'ytick', [0 0.25 0.5 0.75 1]);
axis([0 2*pi+0.05 0 1]);

end