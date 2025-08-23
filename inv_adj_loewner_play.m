function [a_E, a_F] = inv_adj_loewner_play(z, w)

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
    G = zeros(M,M);
    for i = 1:M
        for j = 1:M
            E_inn(i,j) = trace(E{i}'*E{j});
            F_inn(i,j) = trace(F{i}'*F{j});
            G(i,j) = 2 / abs(conj(z(i)) - z(j))^2;
        end
    end

    c = pinv(E_inn)*a;
    d = pinv(F_inn)*b;

    for j = 1:M
        X = X + c(j)*E{j} + d(j)*F{j};
    end

    %%%%%%%%%% To check how the algebraic connectivity of E and F change
    %%%%%%%%%% with the dimension M

    sv_E = svd(E_inn);
    a_E = sv_E(end-1);
    sv_F = svd(F_inn);
    a_F = sv_F(end);
    
    G = F_inn;
    for k = 1:M
        G(k,k) = 2 / abs(conj(z(k)) - z(k))^2;
    end
    D = F_inn - G;
    % fprintf("G = \n");
    % disp(eig(G));
    % min(diag(D));
end

    % k = 1:M;
    % theta = angle(z);
    % delta = theta(2) - theta(1);
    % margin = (pi - delta*M) / 2;
    % row = 1 ./ (1 - cos((pi - 2*margin)*k/M + 2*margin));
    % hankel(row, flip(row));
    % eig(H)
    % 
    % Q = F_inn;
    % r = floor(M/2);
    % I = eye(r);
    % J = flip(I);
    
    % if rem(M,2)==0
    %     A = Q(1:r, 1:r);
    %     C = Q(r+1:end, 1:r);
    % 
    %     K = 1/sqrt(2)*[I -J; I J];
    %     Q2 = K*Q*K';
    %     Qm = A - J*C;
    %     Qp = A + J*C;
    % 
    %     % theta = angle(z);
    %     % delta = theta(2) - theta(1);
    %     % x = linspace(-pi, pi, 1000);
    %     % y = arrayfun(@(x) func(x, r, delta), x);
    %     % plot(x, y)
    %     % grid on
    %     % real(eig(J*C))
    % 
    %     % row = F_inn(end,1:r);
    %     % row = [row, 0, flip(row(2:end))];
    %     % min(real(fft(row)));
    %     % min(real(eig(J*C)));
    %     % min(eig(F_inn));
    %     % lower_bound = min(2*diag(Qm) - sum(Qm, 2)) + 2*min(real(fft(row)));
    % 
    %     % [Vq, Sq] = eig(Q);
    %     % [Vqm, Sqm] = eig(Qm);
    %     % [Vqp, Sqp] = eig(Qp);
    %     % 
    %     % T = J*C;
    %     % alpha = 0;
    %     % row = zeros(1,r);
    %     % row(1) = alpha;
    %     % row(2:end) = fliplr(T(1,2:end));
    %     % S = toeplitz(row);
    %     % Cir = [T, S; S, T];
    %     % F = 1/sqrt(2)*[I, I; I, -I];
    %     % Cir2 = F*Cir*F';
    %     % sort(real(fft(Cir(1,:))));
    %     % sort(real(eig(T)));
    % 
    % else
    %     A  = Q(1:r, 1:r);
    %     x = Q(1:r,r+1);
    %     C = Q(r+2:end,1:r);
    %     q = Q(r+1,r+1);
    % 
    %     K = 1/sqrt(2)*[I, zeros(r,1), -J
    %                    zeros(1,r), sqrt(2), zeros(1,r)
    %                    I, zeros(r,1), J];
    % 
    %     Q2 = K*Q*K';
    %     Qm = A - J*C;
    %     Qp = A + J*C;
    %     Qpxq = [q, sqrt(2)*x'; sqrt(2)*x, Qp];
    % 
    %     [Vq, Sq] = eig(Q);
    %     [Vqm, Sqm] = eig(Qm);
    %     [Vqpxq, Sqpxq] = eig(Qpxq);
    % 
    % end

% function out = func(x, r, delta)
%     k = 1:r-1;
%     terms = cos(k .* x) ./ (1 + cos(delta*k));
%     out = 0.5 + 2*sum(terms);
% end