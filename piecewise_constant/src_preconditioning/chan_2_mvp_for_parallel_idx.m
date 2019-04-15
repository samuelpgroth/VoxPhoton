%% Function for multiplication by 2-level circulant preconditioner
function[mvp] = chan_2_mvp_for_parallel_idx(inv_blocks,JIn,L,M,N,idx)

Vrhs = zeros(3*L*M*N,1);
Vrhs(idx) = JIn;

temp =reshape(Vrhs,L,3*M*N);
temp= fft(temp).'; % transpose is application of permutation matrix

for i=1:L
    TEMP = reshape(temp(:,i),M,3*N);
    TEMP = fft(TEMP).';
    for j=1:M
        TEMP(:,j) = inv_blocks{i,j}*TEMP(:,j);
    end
    TEMP = ifft(TEMP.');
    temp(:,i) = reshape(TEMP,1,3*M*N);
end

temp = ifft(temp.'); % transpose is application of permutation matrix transpose
TEMP = reshape(temp,3*L*M*N,1);
mvp = TEMP(idx);
end