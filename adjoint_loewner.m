function w = adjoint_loewner(z, L)

    [M, ~] = size(L);
    a = zeros(M,1);
    b = zeros(M,1);
    % w = a + b*i

    E = cell(M,1);
    F = cell(M,1);

    for k = 1:M
        w = zeros(M, 1);
        w(k) = 1;
        E{k} = loewner(z, w);
        w = zeros(M, 1);
        w(k) = 1i;
        F{k} = loewner(z, w);
    end

    for k = 1:M
        a(k) = trace(E{k}'*L);
        b(k) = trace(F{k}'*L);
    end
    a = real(a);
    b = real(b);

    w = a + b*1i;

end