%% Function for multiplication by preconditioner
function[mvp] = chan_mvp(inv_blocks,Vrhs,L,M,N)
temp =reshape(Vrhs,L,3*M*N);
temp= fft(temp).'; % transpose is application of permutation matrix

for i=1:L
    temp(:,i) = inv_blocks{i}*temp(:,i);
end

temp = ifft(temp.'); % transpose is application of permutation matrix transpose
mvp = reshape(temp,3*L*M*N,1);
end