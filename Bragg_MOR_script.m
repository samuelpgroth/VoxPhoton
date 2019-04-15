
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