% This code calculates the effective index of a dielectric ring and
% visualize the mode profiles
clear; 
close all;
% 
% Reference: L. Prkna, M. Hubalek, and J. Ctyroky, “Vectorial eigenmode solver for
% bent waveguides based on mode matching,” IEEE Photonics Technology
% Letters, vol. 16, no. 9, pp. 2057–2059, 2004.
% comsol = [ 1    1.912366819227168   1550;
%                   2   1.9041339628267908  1544;
%                   3   1.8936137407819036  1535;
%                   4   1.8791037648567879  1523];
% 
% our_results = [1.9148 1.9017 1.8913 1.8759].';
% p = [comsol(:,2)  our_results(:,1)];
% pp = [p (p(:,1)-p(:,2))./p(:,1)*100] % outputs
% 
%     1.9124    1.9148   -0.1272
%     1.9041    1.9017    0.1278
%     1.8936    1.8913    0.1222
%     1.8791    1.8759    0.1705

%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Prkna2004';
lambda=1550e-9;                   %% Wavelength (meter)
nmodes_search=10;                %% number of solutions asked

%%% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=151;                %% number of points in the direction x
Ny=101;                %% number of points in the direction y
Nyshift = -20;
dX = 50e-9;
dY = 50e-9;

R_central = 200e-6;
ring_width =  1.0e-6;
ring_height = 0.5e-6;
film_thickness = 1e-6;

n_ring = 1.99;
n_film = n_ring;
n_substrate = 1.45;
n_background = 1;

Nr_ring = round(ring_width/dX);
Nz_ring = round(ring_height/dX);
Nz_film = round(film_thickness/dY);

%%% Grid definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = (-(Nx-1)/2:(Nx+1)/2-1)*dX  + R_central;
y = (-(Ny-1)/2:(Ny+1)/2-1)*dY;

dx = dX;
dy = dY;

[X,Y] = meshgrid(x,y);
n = n_background*ones(size(X));

Nsubs_boundary = Nyshift+round(Ny-1)/2;
n(1:Nsubs_boundary,:) = n_substrate;        % substrate
n(Nsubs_boundary+1:Nsubs_boundary+Nz_film,:) = n_film;      % thin film
n(Nsubs_boundary+Nz_film+1: Nsubs_boundary+Nz_film+Nz_ring, ...
    round((Nx-1)/2)+1-Nr_ring+floor(Nr_ring/2):round((Nx-1)/2)+floor(Nr_ring/2)) = n_ring;  % ring
% permittivities
eps=n.^2;

% % Plot Geometry  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_geometry_grid;
print('-dpng','-r100',[filename '_mesh'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 2*pi/lambda;
plot_fields = 0;
neff_min=n_ring*0.8;           %% filter the solutions where the effective index is superior than
neff_max=n_ring*0.99;           %% filter the solutions where the effective index is inferior than


beta_guess = 2*pi/lambda*0.97*n_ring*(1+1e-6i);
[Er, Ephi, Ez, Hr, Hphi, Hz, betas, alphas, neffs, spurious, mvals] = Ring_Solver(R_central,dx,dy,x,y,lambda,eps,nmodes_search, neff_min, neff_max,plot_fields,filename,beta_guess);

disp(neffs(alphas<0));

save('Prkna2004')
!mv *.mat ./mat_files/
!mv *.png ./png_files/