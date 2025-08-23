clear
clc
close all

set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 2);

rng(30);

%%

max_M = 80;
G = drss(8);
M_values = 5:5:max_M;
num_M = length(M_values);
a_E = zeros(num_M, 1);
a_F = zeros(num_M, 1);
Linvadj_values = zeros(num_M, 1);

counter = 0;
for M = M_values
    counter = counter + 1;
    margin = 0.1; % the margin for frequency points
    delta = (pi - 2*margin) / M; % the arc distance between frequency points
    theta = (margin + delta/2):delta:(pi - margin - delta/2);
    z = exp(1i.*theta);
    zbar = conj(z);
    v = freqresp(G, z); v = squeeze(v(1,1,:));
    [a_E(counter), a_F(counter)] = inv_adj_loewner_play(z, v);
    Linvadj = inv_adj_loewner(z, v);
    Linvadj_values(counter) = norm(Linvadj);
    disp(M_values(counter));
end

figure;
hold on;
plot(M_values, Linvadj_values, 'LineWidth', 2)

figure;
hold on
plot(M_values, a_E, 'g-', 'LineWidth', 2)
plot(M_values, a_F, 'b-', 'LineWidth', 2)
plot(M_values, 0.5*(M_values), 'r--')
legend(["Algebraic Connectivity of E", "Algebraic Connectivity of F", "y=0.5x"])