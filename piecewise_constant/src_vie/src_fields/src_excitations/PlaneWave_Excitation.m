function [E,H] = PlaneWave_Excitation(r,k,omega_mu,Eo)
%%    Function to Generate a Plane wave excitation
% _________________________________________________________________________
%
%       Generates the E and H fields of a plane wave
%       The electric field magnitude is 1V/m
% _________________________________________________________________________
%
%% INPUT
%   r               4D (LxMxNx3) array with domain voxelized grid coordinates
%   k               vector with wavenumbers in the medium [kx ky kz]
%   omega_mu        omega*mu
%   Eo              polarization vector for E-field
%
%% OUTPUT
%   E               Electric field (LxMxNx3)
%   H               Magnetic field (LxMxNx3)
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright ? 2014
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Obtain voxel resolution and domain size
% -------------------------------------------------------------------------

[L, M, N, ~] = size(r);

dx = r(2,1,1,1) - r(1,1,1,1);
% dy = r(1,2,1,2) - r(1,1,1,2);
% dz = r(1,1,2,3) - r(1,1,1,3);
% HACK: for one voxel case, the dx and dy expressions don't work
% dx = dz;
dy = dx;
dz = dx;

% -------------------------------------------------------------------------
% polarization H-field
% -------------------------------------------------------------------------

Ho = (1/omega_mu) * cross(k,Eo);

% -------------------------------------------------------------------------
% extract wavenumbers
% -------------------------------------------------------------------------

kx = k(1);
ky = k(2);
kz = k(3);

% -------------------------------------------------------------------------
% traveling direction of excitation
% -------------------------------------------------------------------------
if kx ~= 0
    Ex =  exp(-1j*kx* r(:,:,:,1) ) *  (exp(-1j*kx*dx/2) - exp(1j*kx*dx/2)) / (-1j*kx)  ;
else
    Ex = dx * ones(L,M,N);
end

if ky ~= 0
    Ey = exp(-1j*ky* r(:,:,:,2) ) *  (exp(-1j*ky*dy/2) - exp(1j*ky*dy/2)) / (-1j*ky)  ;
else
    Ey = dy * ones(L,M,N);
end

if kz ~= 0
    Ez = exp(-1j*kz* r(:,:,:,3) ) *  (exp(-1j*kz*dz/2) - exp(1j*kz*dz/2)) / (-1j*kz)  ;
else
    Ez = dz * ones(L,M,N);
end
%
Eexc = Ex .* Ey .* Ez;

% -------------------------------------------------------------------------
% apply polarization
% -------------------------------------------------------------------------

E = [Eo(1)*Eexc(:) ; Eo(2)*Eexc(:) ; Eo(3)*Eexc(:) ];
%
H = [Ho(1)*Eexc(:) ; Ho(2)*Eexc(:) ; Ho(3)*Eexc(:) ];

% -------------------------------------------------------------------------
% reshape and scale to 1V/m
% -------------------------------------------------------------------------

E = reshape(E, L, M, N, 3)./(dx*dy*dz);
H = reshape(H, L, M, N, 3)./(dx*dy*dz);
