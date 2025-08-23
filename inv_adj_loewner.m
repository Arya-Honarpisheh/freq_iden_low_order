function X = inv_adj_loewner(z, w)

    M = length(w);
    a = real(w);
    b = imag(w);
    X = zeros(M,M);
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

    E_inn = zeros(M,M);
    F_inn = zeros(M,M);
    for i = 1:M
        for j = 1:M
            E_inn(i,j) = trace(E{i}'*E{j});
            F_inn(i,j) = trace(F{i}'*F{j});
        end
    end
    
    % s = svd(E_inn);
    % s(end-1)
    % s = svd(F_inn);
    % s(end)

    c = pinv(E_inn)*a;
    d = pinv(F_inn)*b;

    for j = 1:M
        X = X + c(j)*E{j} + d(j)*F{j};
    end

end