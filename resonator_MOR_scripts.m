
%% Perform frequency sweep.
% Solves at each of nG wavelength points to construct the reduced model
% This requires assembling the integral operator and preconditioner at each
% of these points.
% "order" is the number of solves at each wavelength point, creating a
% reduced model of size nG * order.

for i=1:nG
    freq = freq_s(i);
    EMconstants;
    % [Einc_minus,~] = Dipole_Excitation(r,ko,Eo,DipCoord);
    % Obtain coordinates of the domain
    Ocoord = reshape(r,L*M*N,3);
    
    
    % cutoff = 5*lambda_temp;
    % Generate E field
    [E] = E_field_DGF_mollified(Eo,DipCoord,Ocoord,ko,cutoff);
    Einc = reshape(E,L,M,N,3);
    Einc=Einc./max(max(max(abs(Einc))));
    
    tau = 1j*omega*eo*Mop_mid;
    tau3 = [tau(:); tau(:); tau(:)]; % 3 Cartesian components in vector form
    
    VRHS{i} = Gram.*tau3(idxS3).*Einc(idxS3);
    
    [fN{i},opToep{i}] = getOPERATORS_SAM(r,freq,'N',[],'DEMCEM');
    
end

% Calculate the derivatives of fN at the Legendre points

for i=1:nG
    fN_deriv{i} = lagrange_int_deriv(fN,s(i),s);
end

for i=1:nG
    MVP_E{i} = @(J) mv_AN_const(J,fN_deriv{i},ones(size(Mop_mid)),Mop_mid,Gram,...
        'notransp',idxS3,0) - Gram.*J;  
end


%% Create conveniently ordered geometry etc. for top WG alone
% Domain of waveguide
idx_WG3 = [idx_WG;nD+idx_WG;2*nD+idx_WG];
rWG=r(idx_WG3);
num_vox_wide = length(idx_WG)/(L*N);

r_WG = reshape(rWG,L,num_vox_wide,N,3);

%[r_WG] = generatedomain_new(res,dom_WG_x/2,dom_WG_y/2,dom_WG_z/2);
xd_WG = r_WG(:,:,:,1);
yd_WG = r_WG(:,:,:,2);
zd_WG = r_WG(:,:,:,3);

% % compute the circulants
for i=1:nG
    [fN_WG{i},opToeplitz_WG{i}] = getOPERATORS_SAM(r_WG,freq_s(i),'N',[],'DEMCEM');
end
    
[L_WG,M_WG,N_WG,~] = size(r_WG);
                
freq = freq_mid;
EMconstants;

sigma_e_WG = zeros(L_WG,M_WG,N_WG); 
sigma_e_WG = sigma0*fun_absorb(xd_WG);

epsilon_r_WG = ones(L_WG,M_WG,N_WG);
epsilon_r_WG = e_r; 

mu_r_WG = ones(L_WG,M_WG,N_WG);
mu_r_WG = m_r;

epsilon_r_WG = epsilon_r_WG - 1j*sigma_e_WG/(eo*omega);
% % keyboard
% 
%% Assemble the entire matrix to look at e-values and condition number
% Compute the relative permittivity and suceptibility
Mr_WG = epsilon_r_WG;
Mc_WG = epsilon_r_WG - 1.0;

% Domain dimensions
nD_WG = L_WG*M_WG*N_WG; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS_WG = find(abs(Mc_WG(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3_WG = [idxS_WG; nD_WG+idxS_WG; 2*nD_WG+idxS_WG]; % the vector of non-air positions for 3 Cartesian components

Mc_WG(:) = Mc_WG(:)./Mr_WG(:);
% notice that Vrhs is only defined for non-air components
Mr_WG(:) = Mr_WG(:)./Mr_WG(:);


Mc_center = Mc_WG(round(L_WG/2),1,1)*ones(size(Mc_WG));

%keyboard

%% Assemble 1-level preconditioner
for iG=1:nG
    tic
    [circ_1_inv{iG}] = level_1_parallel_func(opToeplitz_WG{iG},Mc_center,Mr_WG,L_WG,...
    M_WG,N_WG,Gram,1);
    disp('First level assembly');
    toc
    for i=1:L
        circ_1_inv{iG}{i}=-circ_1_inv{iG}{i};
    end
end

%%
% Domain of disk
idx_square3 = [idx_square;nD+idx_square;2*nD+idx_square];
rSquare=r(idx_square3);

indy = find((abs(xd(:,1,1)-ring_centre(1))<=rad));
[dim_x,~] = size(indy);

indy = find(((yd(1,:,1)<=dom_y/2-w-gap)));
[~,dim_y] = size(indy);

%keyboard
num_vox_wide = length(idx_square)/(L*N);

r_Sq = reshape(rSquare,dim_x,dim_y,N,3);

%[r_WG] = generatedomain_new(res,dom_WG_x/2,dom_WG_y/2,dom_WG_z/2);
xd_Sq = r_Sq(:,:,:,1);
yd_Sq = r_Sq(:,:,:,2);
zd_Sq = r_Sq(:,:,:,3);

% % compute the circulants
for i=1:nG
    [fN_Sq{i},opToeplitz_Sq{i}] = getOPERATORS_SAM(r_Sq,freq_s(i),'N',[],'DEMCEM');
end

[L_Sq,M_Sq,N_Sq,~] = size(r_Sq);
                
epsilon_r_Sq = ones(L_Sq,M_Sq,N_Sq);
epsilon_r_Sq = e_r; 

% % keyboard
% 
%%
% Compute the relative permittivity and suceptibility
Mr_Sq = epsilon_r_Sq;
Mc_Sq = epsilon_r_Sq - 1.0;

% Domain dimensions
nD_Sq = L_Sq*M_Sq*N_Sq; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS_Sq = find(abs(Mc_Sq(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3_Sq = [idxS_Sq; nD_Sq+idxS_Sq; 2*nD_Sq+idxS_Sq]; % the vector of non-air positions for 3 Cartesian components


Mc_Sq(:) = Mc_Sq(:)./Mr_Sq(:);
% notice that Vrhs is only defined for non-air components
Mr_Sq(:) = Mr_Sq(:)./Mr_Sq(:);

for iG=1:nG
    tic
    [circ_2_inv_T{iG}] = level_2_parallel_func(opToeplitz_Sq{iG},Mc_Sq,Mr_Sq,L_Sq,M_Sq,N_Sq,Gram);
    disp('Second level assembly');
    toc
    for i=1:L_Sq
        for j=1:M_Sq
            circ_2_inv_T{iG}{i,j}=-circ_2_inv_T{iG}{i,j};
        end
    end

end


pointring_outer_in_sq = ring_outer(r_Sq);
pointring_inner_in_sq = ring_inner(r_Sq);

if strcmp(problem_name,'ring')==1
    idx_disk_in_sq = find((pointring_outer_in_sq(:)<=0).*(pointring_inner_in_sq(:)>=0));
elseif strcmp(problem_name,'disk')==1
    idx_disk_in_sq = find((pointring_outer_in_sq(:)<=0));
end
    
nD_sq=L_Sq*M_Sq*N_Sq;
idx_disk_in_sq3 = [idx_disk_in_sq;nD_sq+idx_disk_in_sq;2*nD_sq+idx_disk_in_sq];

idx_disk3 = [idx_disk; nD+idx_disk; 2*nD+idx_disk];
%%

idxSq3 = [idx_square; nD+idx_square; 2*nD+idx_square];
% % prec = @(J) prec_fun(J,mat_inv_WG,idxT3);

for iG=1:nG
%     prec_A{iG} = @(J) prec_fun_disk_idx(J,circ_1_inv{iG},idxT3,idxSq3,L_WG,M_WG,N_WG,circ_2_inv_T{iG},L_Sq,M_Sq,N_Sq);
    prec_A{iG} = @(J) prec_fun_disk_idx(J,circ_1_inv{iG},idx_WG3,idx_disk3,...
    L_WG,M_WG,N_WG,circ_2_inv_T{iG},L_Sq,M_Sq,N_Sq,L,M,N,idxS3,idx_disk_in_sq3);
    MVP_A{iG} = @(J) mv_AN_const(J,fN{iG},ones(size(Mc)),Mop_mid,Gram,...
    'notransp',idxS3,0);
end

%keyboard
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
numIt = 4;  % numIt+1 = number of columns of V

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