function[P] = lagrange_int(fN,S,s)

n = length(s);
P=0;
for i=1:n
    L=1;
    for j=1:n
        if j~=i
            L = L.*(S-s(j))./(s(i)-s(j));
        end
    end
    P = P+L.*fN{i};
end