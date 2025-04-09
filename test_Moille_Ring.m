% This code calculates the effective index of a dielectric ring and
% visualize the mode profiles
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)

% COMSOL  This work    diff (%)
% 1.7909    1.7901    0.0434
% 1.7524    1.7491    0.1910
% 1.6257    1.6229    0.1753
% 1.6092    1.6049    0.2655

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'MoilleRing';
lambda=1550e-9;       %% Wavelength (meter)

R_central = 23.0e-6;
ring_width =  1.5e-6;
ring_height = 0.7e-6;
wg_width = 0.5e-6;
wg_height = ring_height;
gap_width = 0.7e-6;
n_background = 1.444;
n_ring = 1.9781;

nmodes_search=8;               %% number of solutions asked

%%% Optical index Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=201;                %% number of points in the direction x (rho)
Ny=101;                %% number of points in the direction y (z)
Nx_ring = 75;
Ny_ring = 35;

dx = ring_width / Nx_ring;
dy = ring_height / Ny_ring;

%%%% Grid definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = (-(Nx-1)/2*dx:dx:(Nx-1)/2*dx)+R_central;
y = -(Ny-1)/2*dy:dy:(Ny-1)/2*dy;
[X,Y] = meshgrid(x,y);

%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neff_min =n_background;
neff_max = n_ring*0.99;

n = n_background*ones(size(X));     % background
l1 = round(Ny/2)-round(Ny_ring/2);
l2 = round(Nx/2)-round(Nx_ring/2);
n(l1:l1+Ny_ring-1, l2:l2+Nx_ring-1) = n_ring; % ring
% permittivity array
eps=n.^2;
% Plot Geometry  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_geometry_grid; axis equal;
print('-dpng','-r100',[filename '_mesh'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_fields = 1;
k0 = 2*pi/lambda;
m = round(k0*R_central*neff_max);

beta_guess = 2*pi/lambda*0.9*n_ring*(1+1e-6i);
[Er, Ephi, Ez, Hr, Hphi, Hz, betas, alphas, neffs, spurious, mvals] = Ring_Solver(R_central,dx,dy,x,y,lambda,eps,nmodes_search, neff_min, neff_max,plot_fields,filename,beta_guess);

disp(neffs(alphas<0));

save(['results_' filename]);
!mv *.mat ./mat_files/
!mv *.png ./png_files/