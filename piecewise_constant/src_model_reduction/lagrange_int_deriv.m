function[P] = lagrange_int_deriv(fN,S,s)

n = length(s);
P=0;
for i=1:n
    L=0;
    for j=1:n
        if j~=i
            Ld = 1;
            for k=1:n
                if k~=i && k~=j
                    Ld = Ld .* (S-s(k))./(s(i)-s(k));
                end
            end
            L = L + 1./(s(i)-s(j)).* Ld;
        end
    end
    P = P+L.*fN{i};
end


        