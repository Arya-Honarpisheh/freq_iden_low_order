function [what] = LNNM(z, Wtil, tau)


[M, N] = size(Wtil);
zbar = conj(z);

w = sdpvar(M, 1, 'full', 'complex');
wbar = conj(w);
Lw = sdpvar(M, M, 'full');
for i = 1:M
    for j = 1:M
        Lw(i,j) = (wbar(i) - w(j)) / (zbar(i) - z(j));
    end
end

quad_term = 0;
for k = 1:N
    quad_term = quad_term + norm(Wtil(:,k) - w, 2)^2;
end

nuc_term = tau / M * norm(Lw, '*');

Options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug', 1, ...
                      'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3, ... % Default 1e-8
                      'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-3, ...  % Default 1e-8
                      'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-3);     % Default 1e-8

Objective = quad_term + nuc_term;

sol = optimize([], Objective, Options);

if sol.problem == 0
    % Extract and display the optimal w
    what = value(w);
    disp('optimization was succesful!');
else
    disp('Solver failed:');
    sol.info
end

end