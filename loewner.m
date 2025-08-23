function L = loewner(z, w)

    M = length(z);
    zbar = conj(z);
    wbar = conj(w);
    
    L = zeros(M, M);
    
    % Construct the Loewner matrix
    for i = 1:M
        for j = 1:M
            L(i,j) = (wbar(i) - w(j)) / (zbar(i) - z(j));
        end
    end

end
