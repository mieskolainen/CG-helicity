% Spin simulation 1
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


%% Compare spin densities

close all;

for k = 1:2

if (k == 1) % True

rho = [0.476+1i*0.000, 	0.019+1i*0.058, -0.183+1i*0.000;		
       0.019-1i*0.058,  0.047+1i*0.000, -0.019+1i*0.058;	
      -0.183+1i*0.000, -0.019-1i*0.058,  0.476+1i*0.000];
end
if (k == 2) % Fit result

rho = [0.498+1i*0.000, 0.024+1i*0.027,  -0.150+1i*0.000;      
       0.024-1i*0.027, 0.004+1i*0.000,  -0.024+1i*0.027;      
      -0.150+1i*0.000, -0.024-1i*0.027,  0.498+1i*0.000];

end

% Normalize to tr(rho) = 1 (due to rounding of elements above)
rho = rho / trace(rho);

% MC the distribution, histogram, plot
xedge = linspace(-1,1, 70);   % costheta bin edges
yedge = linspace(0,2*pi, 70); % phi bin edges

N = 1e5;
tic;
theta_mother = 0; phi_mother = 0;
[X,weights] = spindecay(rho, s1, s2, T, N, theta_mother, phi_mother);
toc;

H = hist2w(X, weights, xedge, yedge); 
plotrho(xedge, yedge, H, rho);

filename = sprintf('./figs/spinsim_1_%d.pdf', k);
eval(sprintf('print -dpdf %s', filename)); 
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));

end

%% Create unweighted sample

Y = unweight(X,weights);
H = hist2m(Y, xedge, yedge);
plotrho(xedge, yedge, H, rho);

