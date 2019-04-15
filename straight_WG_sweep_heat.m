function[intPowerThru] = ...
    straight_WG_sweep_heat(nG,order,resolution,e_r_left,e_r_rate,...
    width,d,LAM,delW,numPer,perturb,lWG,len_abs,rWG)
tStart=tic; 

%% Set up problem parameters
Si_in_SiO2_parameters; % Establish refractive indices and wavelengths at which to solve
%% Bragg grating parameters
% To get a straight waveguide, we set delW = 0
delW = 0e-9;     % corrugation depth

[r,idx,dom_x,dom_y,dom_z,xStart,res] = Bragg_geometry(LAM,numPer,width,delW,...
    d,lWG,rWG,len_abs,perturb,resolution);

% plot_geometry; % Uncomment to plot the voxels in 3D

freq_mid_index = round(nG/2);

% [Mr_mid,Mc_mid] = straight_params_heat(s_mid,n_ext,'Si',len_abs,r,dom_x,...
%     idx,e_r_left,e_r_rate,y_shift);

[Mr_mid,Mc_mid] = straight_params_heat_lengthways(s_mid,n_ext,'Si',len_abs,...
    r,dom_x,idx,e_r_left,e_r_rate,lWG,rWG);
% keyboard

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
%% Perform frequency sweep.
% Solves at each of nG wavelength points to construct the reduced model
% This requires assembling the integral operator and preconditioner at each
% of these points.
% "order" is the number of solves at each wavelength point, creating a
% reduced model of size nG * order.

for i=1:nG
    freq = freq_s(i);
    EMconstants;
    % Obtain coordinates of the domain
    Ocoord = reshape(r,L*M*N,3);
    
    % Generate E field
    [E] = E_field_DGF_mollified(Eo,DipCoord,Ocoord,ko,cutoff);
    Einc = reshape(E,L,M,N,3);
    Einc=Einc./max(max(max(abs(Einc))));
    
    tau = 1j*omega*eo*Mop_mid;
    tau3 = [tau(idx); tau(idx); tau(idx)]; % 3 Cartesian components in vector form
    
    VRHS{i} = Gram.*tau3.*Einc(idxS3);
    
    [fN{i},opToep{i}] = getOPERATORS_SAM(r,freq,'N',[],'DEMCEM');
    
end

% Calculate the derivatives of fN at the Legendre points

for i=1:nG
    fN_deriv{i} = lagrange_int_deriv(fN,s(i),s);
end

%% Set up G operator
% fN_bar = fN_mid - s_mid.*fN_deriv;
% opToep_bar = opToep_mid - s_mid.*opToep_deriv;

for i=1:nG
    MVP_E{i} = @(J) mv_AN_const(J,fN_deriv{i},ones(size(Mop_mid)),Mop_mid,Gram,...
        'notransp',idxS3,0) - Gram.*J;  
end

Mop_wide = zeros(size(Mop_mid));
x_prec = ceil(xStart/res)+1;
Mop_wide = Mop_mid(x_prec,:,:);

%------
for iG=1:nG
%     tic
%     [circ_inv{iG}] = level_1_parallel_func(opToep{iG},Mop_wide,ones(size(Mop_mid)),L,...
%         M,N,Gram,1);
%     disp('First level assembly');
%     toc
%     for i=1:L
%         circ_inv{iG}{i}=-circ_inv{iG}{i};
%     end
    
    %%
    tic
     [circ_N{iG},circ_L_opToep] = level_1_parallel_func_N(opToep{iG},L,...
        M,N,Gram,1,'on');
    disp('Circulant approximation of N operator');
    hey = real(squeeze(sum(Mr_mid,1)./L));
    face = (hey-1)./hey;
    face3 = [face(:);face(:);face(:)];
    face3mat = diag(face3);
    
    for i=1:L
        circ_temp = Gram*eye(3*M*N)-face3mat*circ_N{iG}{i};
        circ_inv{iG}{i} = inv(circ_temp);
    end
    toc
    disp('Averaged 1-level assembly');
    %%

    prec_A{iG} = @(J) chan_mvp_idx(circ_inv{iG},J,L,M,N,idxS3);
    MVP_A{iG} = @(J) mv_AN_const(J,fN{iG},ones(size(Mop_mid)),Mop_mid,Gram,...
        'notransp',idxS3,0);
end

%----------------------------------------------

%% System solves at interpolation points 
tol = 1e-4;

for iG=1:nG
    tini=tic;
    [xSol{iG},flag,relres,iter,resvec] = gmres(@(J)MVP_A{iG}(J), VRHS{iG},[], tol, 1000, @(J)prec_A{iG}(J));
    tSolve=toc(tini);
    nIts = length(resvec);
    fprintf('1-level preconditioner. Solve time = %.2f [sec] \n',tSolve)
    fprintf('Iteration count = %d \n',nIts);
end

%% Create orthogonal matrix V via Arnoldi
numIt = order-1;  % numIt+1 = number of columns of V

for iG=1:nG
    b = xSol{iG};
    Vtemp(:,1) = b./norm(b);
    for n = 1:numIt
        if mod(n,2)==1
            v_temp = -MVP_E{iG}(Vtemp(:,n));
        else
            v_temp = MVP_E{iG}(Vtemp(:,n));
        end
        tini=tic;
        [v,flag,relres,iter,resvec] = gmres(@(J)MVP_A{iG}(J), v_temp,[], tol, 1000, @(J)prec_A{iG}(J));
        tSolve=toc(tini);
        nIts = length(resvec);
        fprintf('1-level preconditioner. Solve time = %.2f [sec] \n',tSolve)
        fprintf('Iteration count = %d \n',nIts);
        
        for j = 1:n
            h(j,n) = Vtemp(:,j)'*v;
            v = v-h(j,n)*Vtemp(:,j);
        end
        h(n+1,n) = norm(v);
        Vtemp(:,n+1) = v./h(n+1,n);
    end
    V{iG} = Vtemp;
end

P=[];
for iG=1:nG
    P = [P V{iG}];
end

[U_svd,S_svd,V_svd]=svd(P,0);
P=U_svd;


% Calculate Q and R from my notes
for i=1:nG*(numIt+1)
    for iG=1:nG
        Qtemp{iG}(:,i) = MVP_A{iG}(P(:,i));
    end
end

for iG=1:nG
    Q{iG} = P'*Qtemp{iG};
end

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
fprintf('Total time taken = %.2f [sec] \n',tTotal);