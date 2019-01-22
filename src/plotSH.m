%% Plot Spherical Harmonics
%
% input:   lmax  =  Maximum order (e.g. 4)
% output:  fig   =  Figure handle
%
% mikael.mieskolainen@cern.ch, 2017

function fig1 = plotSH(lmax)

fig1 = figure;

% Visualization parameters       
amplitude = 1;
radius = 1e-3;

% Spherical coordinate discretization
step  = 40;
theta = 0:pi/step:pi;                 % Polar angle
phi   = 0:pi/(step/2):2*pi;           % Azimuth angle
[phi,theta] = meshgrid(phi,theta);    % Grid

% Spherical -> Cartesian coordinates
r = radius.*sin(theta);
x = r.*cos(phi);
y = r.*sin(phi);
z = radius.*cos(theta);

for l = 0:lmax
   for m = 0:l  % +- l
        
        % Get the Associated Legendre Polynomials P_l^m(cos(theta))
        Ymn = legendre(l,cos(theta(:,1)));
        Ymn = Ymn(m+1,:)';
        
        % Complex Spherical Harmonic Y_l^m(cos(theta), phi)
        c = repmat(Ymn, 1, size(theta,1));
        c = c.*exp(1i*m*phi);
        
        % Get spherical representation
        rho = radius + amplitude*real(c);
        
        % Plot it
        subplot(lmax+1,lmax+1,l*(lmax+1)+m+1)
        surf(x,y,z, rho);
        title(sprintf('$\\ell=%0.0f,m=\\pm%0.0f$', l, m),'interpreter','latex', 'fontsize', 8);
        shading interp; colormap hot; axis equal; axis off; axis square;
        view(0,30);
   end
end

end