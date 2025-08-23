function w = inverse_loewner(z, L)

    [M, ~] = size(L);
    zbar = conj(z);
    a = zeros(M,1);
    b = zeros(M,1);
    % w = a + b*i

    for k = 1:M
        b(k) = 1/(-2*1i)*L(k,k)*(zbar(k) - z(k));
    end
    b = real(b);
    
    a(1) = 0;
    for j = 2:M
        a(j) = a(1) - b(1)*1i - L(1,j)*(zbar(1) - z(j)) - b(j)*1i;
    end
    a = real(a);

    a = a - sum(a)/M;
    
    w = a + b*1i;

end