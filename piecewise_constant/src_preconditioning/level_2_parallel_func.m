function[circ_2_inv_T] = level_2_parallel_func(opToeplitz,Mc,Mr,L,M,N,Gram)

circ_2_inv_T = cell(L,M);

A(:,:,:) = -Mc.*opToeplitz(:,:,:,1);
% Fix up for 1st element
A(1,1,1) = Gram*Mr(1)+A(1,1,1);

clear c
% Now construct circulant approximation

c = zeros(L,M,N);

for i=2:L
    c(i,:,:) = (L+1-i)/L.*A(i,:,:)+(i-1)/L.*A(L-i+2,:,:);
end

% Fix up for 1st element
c(1,:,:) = A(1,:,:);

cc=fft(c);

%% 2nd block
B(:,:,:) = -Mc.*opToeplitz(:,:,:,2);
% Now construct circulant approximation

C = zeros(L,M,N);
for i=2:L
    C(i,:,:) = (L+1-i)/L.*B(i,:,:)-(i-1)/L.*B(L-i+2,:,:);
end
C(1,:,:) = B(1,:,:);

CC=fft(C);

%% 3rd block
D(:,:,:) = -Mc.*opToeplitz(:,:,:,3);
C3 = zeros(L,M,N);
% Now construct circulant approximation

for i=2:L
    C3(i,:,:) = (L+1-i)/L.*D(i,:,:)-(i-1)/L.*D(L-i+2,:,:);
end
C3(1,:,:) = D(1,:,:);
C3=fft(C3);

%% 4th block
E(:,:,:) = -Mc.*opToeplitz(:,:,:,4);
E(1,1,1) = Gram*Mr(1)+E(1,1,1);
C4 = zeros(L,M,N);
% Now construct circulant approximation

for i=2:L
    C4(i,:,:) = (L+1-i)/L.*E(i,:,:)+(i-1)/L.*E(L-i+2,:,:);
end
C4(1,:,:) = E(1,:,:);

C4 = fft(C4);

%% 5th block
clear F
F(:,:,:) = -Mc.*opToeplitz(:,:,:,5);
C5 = zeros(L,M,N);
% Now construct circulant approximation

for i=2:L
    C5(i,:,:) = (L+1-i)/L.*F(i,:,:)+(i-1)/L.*F(L-i+2,:,:);
end
C5(1,:,:) = F(1,:,:);

C5 = fft(C5);

%% 6th block
clear G
G(:,:,:) = -Mc.*opToeplitz(:,:,:,6);
G(1,1,1) = Gram*Mr(1)+G(1,1,1);
C6 = zeros(L,M,N);
% Now construct circulant approximation

for i=2:L
    C6(i,:,:) = (L+1-i)/L.*G(i,:,:)+(i-1)/L.*G(L-i+2,:,:);
end
C6(1,:,:) = G(1,:,:);

C6 = fft(C6);

for i_loop=1:L
    temp = zeros(3*M*N,3*M*N);
    level_2_diag_temp = zeros(3*M*N,3*M*N);
    
    chan = cell(N,1);
    % First block
    result=[];
    for j=1:1
        for i=1:N
            chan{i}=toeplitz(cc(i_loop,1:M,i,j),cc(i_loop,1:M,i,j));
        end
        result=cell2mat(chan(toeplitz(1:N)));
    end
    
    temp(1:M*N,1:M*N) = result;
    
    % Second block
    for j=1:1
        for i=1:N 
            chan{i}=toeplitz([-CC(i_loop,1,i,j) CC(i_loop,2:M,i,j)],-CC(i_loop,1:M,i,j));
        end
        result=cell2mat(chan(toeplitz(1:N)));
    end
    
    temp(1:M*N,M*N+1:2*M*N) = result;
    temp(M*N+1:2*M*N,1:M*N) = result;
    
    % Third block
    for j=1:1
        for i=1:N
            chan{i}=toeplitz(-C3(i_loop,1:M,i,j),-C3(i_loop,1:M,i,j));
        end
        result=cell2mat(chan(toeplitz(1:N)));
    end
    
    Upper = triu(result);
    result = Upper-Upper.'+diag(diag(Upper));
    
    temp(1:M*N,2*M*N+1:3*M*N) = result;
    temp(2*M*N+1:3*M*N,1:M*N) = result;
    
    % Fourth block
    for j=1:1
        for i=1:N
            chan{i}=toeplitz(C4(i_loop,1:M,i,j),C4(i_loop,1:M,i,j));
        end
        result=cell2mat(chan(toeplitz(1:N)));
    end
    
    temp(M*N+1:2*M*N,M*N+1:2*M*N) = result;
    
    % Fifth block
    for j=1:1
        for i=1:N
            chan{i}=toeplitz([C5(i_loop,1,i,j) -C5(i_loop,2:M,i,j)],C5(i_loop,1:M,i,j));
        end
        result=cell2mat(chan(toeplitz(1:N)));
    end
    
    
    Upper = triu(result);
    result = Upper+Upper.'-diag(diag(Upper));
    
    
    temp(M*N+1:2*M*N,2*M*N+1:3*M*N)=result;
    temp(2*M*N+1:3*M*N,M*N+1:2*M*N)=result;
    
    % Sixth block
    
    for j=1:1
        for i=1:N
            chan{i}=toeplitz(C6(i_loop,1:M,i,j),C6(i_loop,1:M,i,j));
        end
        result=cell2mat(chan(toeplitz(1:N)));
    end
    
    temp(2*M*N+1:3*M*N,2*M*N+1:3*M*N)=result;
    
    %% Now construct 2nd level. NB: this can be optimized!!
    c = zeros(1,M);
        for k=1:3*N
            for j=1:3*N
                row_vec = temp((j-1)*M+1,(k-1)*M+1:k*M);
                col_vec = temp((j-1)*M+1:j*M,(k-1)*M+1);      
                c(1) = col_vec(1);
                for i=2:M
                    %if i==1
                    %    c(i) = col_vec(i);
                    %else
                        
                        c(i) = (M+1-i)/M * col_vec(i) +...
                            (i-1)/M * row_vec(M-i+2);
                    %end
                end
                level_2_diag_temp((j-1)*M+1:j*M,(k-1)*M+1:k*M)=...
                    diag(fft(c));
            end
        end

   
        
        diag_block_collapse_temp = reshape(reshape(level_2_diag_temp.',3*M*M*N,3*N).',3*M*N,3*M*N);
%   keyboard
    
    for m=1:M
        circ_2_inv_T{i_loop,m} = inv(diag_block_collapse_temp((m-1)*3*N+1:m*3*N,(m-1)*3*N+1:m*3*N));
    end
end
