function [Erho, Ephi, Ez, Hrho, Hphi, Hz, betas,alphas, neffs, spurious, mvals] = Ring_Solver(R_central,drho,dz,rho,z,lambda,eps,nmodes, neff_min, neff_max,plot_fields, filename,beta_guess)
mu0 = 4e-7 * pi;
eps0 = 8.854e-12;
freq = 2.998e8/lambda;
omega = 2*pi*freq;
k0 = 2*pi/lambda;

[Nz, Nr] = size(eps);
Nrz = Nr * Nz;
[rhos, ~] = meshgrid(rho,z);

boundary_idx = false(Nz, Nr);
boundary_idx(1,:) = true;        % bottom boundary (z=const)
boundary_idx(end,:) = true;      % top boundary (z=const)
boundary_idx(:,1) = true;        % left boundary (rho=const)
boundary_idx(:,end) = true;      % right boundary (rho=const)
boundary_idx = boundary_idx(:);  % flatten to vector for use with matrices


AA=ones(1,Nr*Nz);
BB=ones(1,Nr*Nz-1);
Arz=ones(1,(Nr-1)*Nz);

% First derivative (fourth-order accurate)
DR1 = (-1/12) * spdiags(Arz, -2*Nz, Nrz, Nrz) + (8/12) * spdiags(Arz, -Nz, Nrz, Nrz) ...
    - (8/12) * spdiags(Arz, Nz, Nrz, Nrz) + (1/12) * spdiags(Arz, 2*Nz, Nrz, Nrz);

DZ1 = (-1/12) * spdiags(BB, -2, Nrz, Nrz) + (8/12) * spdiags(BB, -1, Nrz, Nrz) ...
    - (8/12) * spdiags(BB, 1, Nrz, Nrz) + (1/12) * spdiags(BB, 2, Nrz, Nrz);

% Second derivative (fourth-order accurate)
DR2 = (-1/12) * spdiags(Arz, -2*Nz, Nrz, Nrz) + (16/12) * spdiags(Arz, -Nz, Nrz, Nrz) ...
    + (-30/12) * spdiags(AA, 0, Nrz, Nrz) + (16/12) * spdiags(Arz, Nz, Nrz, Nrz) ...
    + (-1/12) * spdiags(Arz, 2*Nz, Nrz, Nrz);

DZ2 = (-1/12) * spdiags(BB, -2, Nrz, Nrz) + (16/12) * spdiags(BB, -1, Nrz, Nrz) ...
    + (-30/12) * spdiags(AA, 0, Nrz, Nrz) + (16/12) * spdiags(BB, 1, Nrz, Nrz) ...
    + (-1/12) * spdiags(BB, 2, Nrz, Nrz);

%% Building of the Hamiltonian
eps_diag = spdiags(eps(:), 0, Nrz, Nrz);
inv_eps_diag = spdiags(1./eps(:), 0, Nrz, Nrz);

DR1 = DR1/drho;
DR2 = DR2/drho^2;

DZ1 = DZ1/dz;
DZ2 = DZ2/dz^2;

% We need to handle internal points and apply fourth-order stencil
deps_drho_inveps= [zeros(Nz,2) -eps(:,1:end-4)+8*eps(:,2:end-3)-8*eps(:,4:end-1) + eps(:,5:end) zeros(Nz,2) ]./eps / (12 * drho);
deps_dz_inveps = [zeros(2,Nr); -eps(1:end-4, :)+8*eps(2:end-3, :)-8*eps(4:end-1, :) + eps(5:end, :); zeros(2,Nr) ]./eps / (12 * dz);
depsdr_inveps = spdiags(deps_drho_inveps(:), 0, Nrz, Nrz);
depsdz_inveps = spdiags(deps_dz_inveps(:), 0, Nrz, Nrz);

norR = spdiags(rhos(:), 0, Nrz, Nrz);

rh2Rc2 =norR.^2/R_central^2;                           % rho^2/R_c^2
rhRc2 =norR/R_central^2;                                  % rho/R_c^2
oneRc2 =spdiags(ones(size(rhos(:))), 0, Nrz, Nrz)/R_central^2;   % 1/R_c^2
oneRc =spdiags(ones(size(rhos(:))), 0, Nrz, Nrz)/R_central;   % 1/R_c^2

%% Building of the Hamiltonian 
Laplacian = rh2Rc2*(DR2 + DZ2 + eps_diag*k0^2)+ rhRc2*DR1 ;
M1 = Laplacian - oneRc2 +rh2Rc2*DR1*depsdr_inveps;
M2 = rh2Rc2*DR1*depsdz_inveps;
M5 = rh2Rc2*DZ1*depsdr_inveps;
M6 = Laplacian + rh2Rc2*DZ1*depsdz_inveps;
M7 = zeros(size(M1));
M11 = Laplacian-oneRc2-rh2Rc2*depsdz_inveps*DZ1;
M12= rh2Rc2*depsdz_inveps*DR1;
M15 = rh2Rc2*depsdr_inveps*DZ1;
M16 =  Laplacian-rh2Rc2*depsdr_inveps*DR1;

M3 = 2/omega/eps0*inv_eps_diag*oneRc*DZ1;
M4 = -2/omega/eps0*inv_eps_diag*oneRc*DR1;
M9 = -2/omega/mu0*oneRc*DZ1;
M10 = 2/omega/mu0*oneRc*DR1;

%% Diagonalization of the Hamiltonian
M = sparse([M1 M2 M7 M7; M5 M6 M7 M7; M7 M7 M11 M12; M7 M7 M15 M16]);
T = sparse([M7 M7 M3 M4; M7 M7 M7 M7; M9 M10 M7 M7; M7 M7 M7 M7]);

I = speye(4 * Nrz);

% Construct the augmented matrices for linear GEP
A = [M, T; sparse(4 * Nrz, 4 * Nrz), speye(4 * Nrz)];
B = [sparse(4 * Nrz, 4 * Nrz), speye(4 * Nrz); speye(4 * Nrz), sparse(4 * Nrz, 4 * Nrz)];

% Expand the boundary mask for each field component
bc_mask = repmat(boundary_idx, 4, 1);

% Apply Dirichlet BCs by zeroing rows/cols and inserting identity on diagonal
A(bc_mask, :) = 0;
A(:, bc_mask) = 0;
A(sub2ind(size(A), find(bc_mask), find(bc_mask))) = 1;

B(bc_mask, :) = 0;
B(:, bc_mask) = 0;

% Solve the linear GEP using eigs (for sparse matrices)
opts.tol = 1e-12;
opts.maxit = 1000;
[V, D] = eigs(A, B, 3*nmodes,beta_guess,opts);

% Extract eigenvalues (β) and eigenvectors (E)
betas = real(diag(D));
alphas = imag(diag(D));
psis = V(1:4 * Nrz, :);  % First 4*Nrz rows correspond to E

% Filter solutions within the desired neff range
neffs = real(betas) / k0;
[neffs, sort_index]  = sort(neffs, 'descend');
betas = betas(sort_index);
alphas = alphas(sort_index);
psis = psis(:, sort_index);

valid_modes = (neffs >= neff_min) & (neffs <= neff_max);
if length(valid_modes)<1
    disp('no solution')
end
betas = betas(valid_modes);
alphas = alphas(valid_modes);
psis = psis(:, valid_modes);
neffs = neffs(valid_modes)
neff = neffs;

nmodes_original = nmodes;
nmodes=length(neff);

Erho = zeros(Nz, Nr, nmodes);
Ephi = zeros(Nz, Nr, nmodes);
Ez = zeros(Nz, Nr, nmodes);
Hrho = zeros(Nz, Nr, nmodes);
Hphi = zeros(Nz, Nr, nmodes);
Hz = zeros(Nz, Nr, nmodes);

first_termE = inv_eps_diag./ (1i * omega*eps0);
first_termH = 1./ (1i * omega*mu0);

spurious = [];
mvals= [];
for imo = 1:min([nmodes, nmodes_original])
        m = round(real(betas(imo))*R_central);
        mvals = [mvals m];
   
        % Extract electric field components
        Erho1 = psis(1:Nrz,imo);        
        Ez1 = psis(Nrz+1:2*Nrz,imo);
        Hrho1 = psis(2*Nrz+1:3*Nrz,imo);
        Hz1 = psis(3*Nrz+1:4*Nrz,imo);
        
        Ephi1 = first_termE * (DZ1* Hrho1 - DR1*Hz1);        
        Hphi1 = first_termH * (DR1* Ez1 - DZ1*Erho1);

        % Reshape results for output
        Hphi(:,:,imo) = reshape(Hphi1, [Nz, Nr]);
        Hrho(:,:,imo) = reshape(Hrho1, [Nz, Nr]);
        Hz(:,:,imo) = reshape(Hz1, [Nz, Nr]);
        Ephi(:,:,imo) = reshape(Ephi1, [Nz, Nr]);
        Erho(:,:,imo) = reshape(Erho1, [Nz, Nr]);
        Ez(:,:,imo) = reshape(Ez1, [Nz, Nr]);

        % power along the direction of propagation
        Power_phi = real(sum(sum(Ez(:,:,imo).*conj(Hrho(:,:,imo))-Erho(:,:,imo).*conj(Hz(:,:,imo)))))*drho*dz/lambda^2;
        % these should be very small numbers
        Power_r = real(sum(sum(Ephi(:,:,imo).*conj(Hz(:,:,imo))-Ez(:,:,imo).*conj(Hphi(:,:,imo)))))*drho*dz/lambda^2;
        Power_z = real(sum(sum(Erho(:,:,imo).*conj(Hphi(:,:,imo))-Ephi(:,:,imo).*conj(Hrho(:,:,imo)))))*drho*dz/lambda^2;

        normalizer = sqrt(abs(Power_phi));
        Hrho(:,:,imo) = Hrho(:,:,imo)/normalizer;
        Hphi(:,:,imo) = Hphi(:,:,imo)/normalizer;
        Hz(:,:,imo) = Hz(:,:,imo)/normalizer;
        Erho(:,:,imo) = Erho(:,:,imo)/normalizer;
        Ephi(:,:,imo) = Ephi(:,:,imo)/normalizer;
        Ez(:,:,imo) = Ez(:,:,imo)/normalizer;        

        % check spurious or not
        % idea: check where the total power is strongest
        % if it is left or right upper or lower corner, then mark the mode
        % spurious!
        Etotal = sqrt(abs(Erho(:,:,imo)).^2+abs(Ez(:,:,imo)).^2+abs(Ephi(:,:,imo).^2));
        z9 = round(Nz/3);
        r9 = round(Nr/3);
        E1 = sum(sum( Etotal(1:z9,1:r9)));
        E2 = sum(sum( Etotal(1:z9,r9+1:2*r9)));
        E3 = sum(sum( Etotal(1:z9,2*r9+1:end)));
        E4 = sum(sum( Etotal(z9+1:2*z9,1:r9)));
        E5 = sum(sum( Etotal(z9+1:2*z9,r9+1:2*r9)));
        E6 = sum(sum( Etotal(z9+1:2*z9,2*r9+1:end)));
        E7 = sum(sum( Etotal(2*z9+1:end,1:r9)));
        E8 = sum(sum( Etotal(2*z9+1:end,r9+1:2*r9)));
        E9 = sum(sum( Etotal(2*z9+1:end,2*r9+1:end)));
        [~, whichone]= max([E1 E2 E3 E4 E5 E6 E7 E8 E9]);
        fig_message = ' Looks Normal ';
        spurious = [spurious 0];
        if whichone == 1 || whichone == 3  || whichone == 7  || whichone == 9
            disp('**** WARNING: probably a spurious mode ****')
            fig_message = 'Spurious Mode Risk';
            spurious(end) = 1;
        end
        
        if plot_fields == 1
            figure('position',[200+20*imo, 1000-20*imo, 1200, 600],'Name',fig_message);
            subplot(231); pcolor(rho*1e6,z*1e6,(abs(Erho(:,:,imo)))); colorbar; shading interp; title('|E_\rho|');            
            subplot(232); pcolor(rho*1e6,z*1e6,(abs(Ephi(:,:,imo)))); colorbar; shading interp; title('|E_\phi|')
            subplot(233); pcolor(rho*1e6,z*1e6,(abs(Ez(:,:,imo)))); colorbar; shading interp; title('|E_z|')
            subplot(234); pcolor(rho*1e6,z*1e6,(abs(Hrho(:,:,imo)))); colorbar; shading interp; title('|H_\rho|')
            subplot(235); pcolor(rho*1e6,z*1e6,(abs(Hphi(:,:,imo)))); colorbar; shading interp; title('|H_\phi|')
            subplot(236); pcolor(rho*1e6,z*1e6,(abs(Hz(:,:,imo)))); colorbar; shading interp; title('|H_z|')
            print('-dpng','-r100',[filename int2str(imo)])
        end
    end
end