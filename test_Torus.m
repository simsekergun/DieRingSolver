% This code calculates the effective index of a dielectric torus and
% visualize the mode profiles

clear; 
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Torus';
lambda=1550e-9;       %% Wavelength (meter)

nmodes_search=30;                %% number of solutions asked

%%% Optical index Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=121;                %% number of points in the direction x
Ny=121;                %% number of points in the direction y

R_major = 20e-6;
r_minor = 0.6e-6;

n_ring = 1.9761;
n_background = 1.444;
neff_max=0.99*n_ring;    %% filter the solutions where the effective index is superior than
neff_min=n_background;                     %% filter the solutions where the effective index is inferior than

dx = 30e-9;
dy = 30e-9;

%%% Grid definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = (-(Nx-1)/2:(Nx+1)/2-1)*dx  + R_major;
y = (-(Ny-1)/2:(Ny+1)/2-1)*dy;

[X,Y] = meshgrid(x,y);
n = n_background*ones(size(X));                          % background
n((X-R_major).^2+Y.^2<r_minor.^2) = n_ring;     % ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps=n.^2;
% Plot Geometry  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_geometry_grid;
print('-dpng','-r100',[filename '_mesh'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_central = R_major;
k0 = 2*pi/lambda;
plot_fields = 1;

beta_guess = 2*pi/lambda*0.9*n_ring*(1+1e-6i);
[Er, Ephi, Ez, Hr, Hphi, Hz, betas, alphas, neffs, spurious, mvals] = Ring_Solver(R_central,dx,dy,x,y,lambda,eps,nmodes_search, neff_min, neff_max,plot_fields,filename,beta_guess);

disp(neffs(alphas<0));

save(['results_' filename]);
!mv *.mat ./mat_files/
!mv *.png ./png_files/