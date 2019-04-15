%-------------------------------------------------------------------------%
%
%                   Frequency sweep for Bragg grating 
%
%   Calculate the transmission spectrum between 1500nm and 1600nm.
%
%   S.P.Groth 14th April 2019
%
%-------------------------------------------------------------------------%
close all;
clear;
format compact;

addpath(genpath('piecewise_constant')) % add the directory path
tStart = tic;

%% Bragg grating geometrical parameters 
% These parameters are chosen to be discretized exactly by 25nm voxels.
% Modify if you desire finer discretizations.
width = 500e-9;    % channel width: 510nm
d = 225e-9;        % channel depth: 210nm
LAM = 325e-9;      % period length
delW = 50e-9;      % corrugation depth
numPer = 100;       % number of periods
resolution = 25e-9;% voxel size

% Model order reduction parameters
nG = 5;            % number of wavelengths at which to solve for freq. sweep
order = 5;         % model order at each of nG points

lWG = 4.5e-6;      % straight portion of waveguide to left of grating
rWG = 2.25e-6;     % straight portion of waveguide to right of grating
perturb = 0;       % number of voxels by which to perturb width
len_abs = 2.25e-6; % absorber length

% Generate voxels
[r,idx,dom_x,dom_y,dom_z,xStart,res] = Bragg_geometry(LAM,numPer,width,delW,...
    d,lWG,rWG,len_abs,perturb,resolution);

% plot_geometry;  % uncomment to plot voxels in 3D

Si_in_SiO2_parameters; % Establish refractive indices and wavelengths at which to solve

% freq_mid_index = round(nG/2);

%% Generate permittivity profile for Bragg grating
% Temperature gradient profiles - hack within these functions to modify!!
T_diff = 0;   % temperature difference, 0 gives homogeneous grating
% Variation along length
[Mr_mid,Mc_mid,e_r_left,e_r_rate] = ...
    Bragg_params_heat_lengthways(s_mid,n_ext,'Si',len_abs,...
    r,dom_x,idx,lWG,rWG,T_diff);

% plot a 2D slice of the permittivity profile
figure
imagesc(rot90(abs(Mr_mid(:,:,2))))

Mop_mid = Mc_mid./Mr_mid;
% get the voxel dimension
dx = r(2,1,1,1) - r(1,1,1,1);
% Gram matrix (value)
Gram = dx^3; % volume of the voxel
%% Incident field - Gaussian beam (complex dipolem source)
Eo = [0,1,0]; % y-polarized incident E-field
lambda_temp = 440e-9;
DipCoord = [1/4*lambda_temp+len_abs - 1j*lambda_temp 0 0];
cutoff = 0.5*lambda_temp;

% Domain dimensions
[L,M,N,~] = size(r);
nD = L*M*N; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS = find(abs(Mop_mid(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components

%% Construct reduced model of size nG * order
% This generates the required nG * order operators and preconditioners
Bragg_MOR_script;

%% Now solve reduced model at finely sampled points in 1500nm-1600nm range
% Calculate the transmitted power through the Bragg grating
intPowerThru = zeros(size(free_wavelength));
x = zeros(length(idxS3),length(free_wavelength));
for i_lambda = 1:numLambda
    S = free_wavelength(i_lambda);
    freq=free_freq(i_lambda);
    EMconstants;
    % Incident field and right hand side
    [E] = E_field_DGF_mollified(Eo,DipCoord,Ocoord,ko,cutoff);
    Einc_temp = reshape(E,L,M,N,3);
    Einc=Einc_temp./max(max(max(abs(Einc_temp))));
 
    tau = 1j*omega*eo*Mop_mid;
    tau3 = [tau(idx); tau(idx); tau(idx)]; % 3 Cartesian components in vector form
    Vrhs = Gram.*tau3.*Einc(idxS3);
    RHS = P'*Vrhs;                     % project down right-hand side to reduced subspace
    
    x_temp = lagrange_int(Q,S,s)\RHS;  % solve reduced system  
    x(:,i_lambda) = P*x_temp;          % project solution back up
    
    Je = zeros(L,M,N,3);
    Je(idxS3) = x(:,i_lambda) ;
    fN_temp=lagrange_int(fN,S,s);
    [eOut] = E_field_Nop_const(Je,fN_temp,Gram,freq,Einc);
    
    %% Transmitted power calculation across voxel slice at end of Bragg
    xThru = L-round((len_abs+450e-9)/res); % approx one wavelength from absorber
    nW = round(width/res); % number of voxels across waveguide
    voxFromEdge = round(delW/2)/res;
    eThru = squeeze(eOut(xThru,1+voxFromEdge:M-voxFromEdge,:,:));
    powerThru = sum(eThru.*conj(eThru),3);
    intPowerThru(i_lambda) = sum(sum(powerThru)); % transmitted power   
end

tTotal = toc(tStart);
fprintf('Total time taken for Bragg simulation = %.2f [sec] \n',tTotal);

% In order to obtain the normalized transmission, we must perform the same
% frequency sweep calculation but now for a straight waveguide. This is 
% relatively inexpensive since we only require a small reduced model.
[intPowerThruStraight] = ...
    straight_WG_sweep_heat(nG,order,resolution,e_r_left,e_r_rate,...
    width,d,LAM,delW,numPer,perturb,lWG,len_abs,rWG);

%% Plot transmitted power
normalizedPower = intPowerThru./intPowerThruStraight;
normalizedPower = normalizedPower./max(normalizedPower);
figure
plot(free_wavelength*1e9,normalizedPower)
xlabel('Wavelength (nm)')
ylabel('Transmission')

%% Save output
save('Bragg_N50_test','normalizedPower','free_wavelength','intPowerThru',...
    'intPowerThruStraight')

