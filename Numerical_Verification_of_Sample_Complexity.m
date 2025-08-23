
clear
clc
close all

set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 2);

red = [0.75,0,0];
green = [0,0.5,0];
blue = [0,0,0.65];

rng(0);

%%

Gbar = tf([0.12 0.18], [1 -1.4 1.443 -1.123 0.7729], 1); % true system

M = 32; % the number of frequency points
N = 30; % the number of samples at each frequency
delta = 0.1; % the margin for frequency points
eta = 0.5; % the bound on noise
arc_dist = (pi - 2*delta) / M; % the arc distance between frequency points
theta = (delta + arc_dist/2):arc_dist:(pi - delta - arc_dist/2);
z = exp(1i.*theta);
zbar = conj(z);

figure;
hold on
axis equal
axis([-1.1, 1.1, -1.1, 1.1]);
rectangle('Position', [-1,-1,2,2], 'Curvature', 1, 'LineStyle', '--', 'LineWidth', 1, 'EdgeColor', 'k');
plot([-1, 1], [0, 0], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
plot(real(z), imag(z), '+', 'MarkerSize', 8, 'LineWidth', 3, 'Color', red);
plot(real(zbar), imag(zbar), '+', 'MarkerSize', 8, 'LineWidth', 3, 'Color', red);
title("Frequency points");

wbar = freqresp(Gbar, z); wbar = squeeze(wbar(1,1,:));
Wbar = repmat(wbar, 1, N);
V = eta*(2*rand(M, N)-1) + 1i*eta*(2*rand(M, N)-1);
Wtil = Wbar + V;

frdatabar = idfrd(wbar, theta, 1);
[magbar, phasebar] = bode(frdatabar);
magbar = squeeze(magbar);
phasebar = squeeze(phasebar);

magtil = zeros(M, N);
phasetil = zeros(M, N);
for k = 1:N
    [magtil(:,k), phasetil(:,k)] = bode(idfrd(Wtil(:,k), theta, 1));
end
magtil_mean = mean(magtil, 2);
magtil_lower = prctile(magtil, 2.5, 2);  % 2.5th percentile (M x 1)
magtil_upper = prctile(magtil, 97.5, 2); % 97.5th percentile (M x 1)
err_lower_mag = magtil_mean - magtil_lower;    % Error bar (lower side)
err_upper_mag = magtil_upper - magtil_mean;    % Error bar (upper side)
phasetil_mean = mean(phasetil, 2);
phasetil_lower = prctile(phasetil, 2.5, 2);  % 2.5th percentile (M x 1)
phasetil_upper = prctile(phasetil, 97.5, 2); % 97.5th percentile (M x 1)
err_lower_phase = phasetil_mean - phasetil_lower;    % Error bar (lower side)
err_upper_phase = phasetil_upper - phasetil_mean;    % Error bar (upper side)

figure();
subplot(2,1,1);
hold on;
plot(theta, 20*log10(magbar));
errorbar(theta, 20*log10(magtil_mean), 20*log10(err_lower_mag), 20*log10(err_upper_mag), 'LineWidth', 2);
set(gca, 'XScale', 'log');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title("Bode Magnitude Plot");
grid on;
subplot(2,1,2);
hold on;
plot(theta, phasebar);
errorbar(theta, phasetil_mean, err_lower_phase, err_upper_phase, 'LineWidth', 2);
set(gca, 'XScale', 'log');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
title('Bode Phase Plot');
grid on;

Lwbar = loewner(z, wbar);
sv_Lwdot = svd(Lwbar);
Lwtil = loewner(z, Wtil(:,1));
sv_Lwtil = svd(Lwtil);

figure;
hold on
title('Singular Values of the Loewner Matrix');
bar([sv_Lwdot, sv_Lwtil]);
legend(["True", "Noisy"]);

%%
% This part is to verify the formal derivations of loewner-related
% matrices in the paper.
Gpert = tf([0.12 0.18]*(1 + 5e-2), [1 -1.4 1.443+0.12 -1.123-0.02 0.7729+0.001], 1); % tiny num perturb
v = squeeze(freqresp(Gpert, z)) - wbar;
v = v - mean(v);

Lv = loewner(z, v);
Lvadj = adjoint_loewner(z, Lv);
Lvinvadj = inv_adj_loewner(z, v);

v - real(mean(v))
adjoint_loewner(z, Lvinvadj)

Lwbar = loewner(z, wbar);
Lwbaradj = adjoint_loewner(z, Lwbar);
Lvadj_real = [real(Lvadj); imag(Lvadj)];
wbar_real = [real(wbar); imag(wbar)];
trace(Lwbar' * Lv) - trace(wbar_real' * Lvadj_real)

%%

% Hinf vs N
num_exp = 20;
num_N = 8;
tau = 7;
N_values = floor(logspace(1, log10(300), num_N));
M = 32; % the number of frequency points
delta = 0.1; % the margin for frequency points
eta = 0.5; % the bound on noise
arc_dist = (pi - 2*delta) / M; % the arc distance between frequency points
theta = (delta + arc_dist/2):arc_dist:(pi - delta - arc_dist/2);
z = exp(1i.*theta);
zbar = conj(z);
wbar = freqresp(Gbar, z); wbar = squeeze(wbar(1,1,:));

id_errors = zeros(num_N, num_exp);

for roll = 1:num_exp
    counter = 0;
    for N = N_values
        counter = counter + 1;
        fprintf("exp = %d, N = %.2d \n", roll, N);
        Wbar = repmat(wbar, 1, N);
        V = eta*(2*rand(M, N)-1) + 1i*eta*(2*rand(M, N)-1);
        Wtil = Wbar + V;
        what = LNNM(z, Wtil, tau);
        id_errors(counter, roll) = max(abs(what - wbar));
    end
end

filename = sprintf("results/N_dep/M%d_eta%d_delta%d_num_exp%d_num_N%d_tau_%d.mat", M, eta, delta, num_exp, num_N, tau);
save(filename, 'id_errors', 'N_values');

%%

load("results/N_dep/M32_eta5.000000e-01_delta1.000000e-01_num_exp20_num_N8_tau_7.mat")

error_mean = mean(id_errors, 2);
error_lower = prctile(id_errors, 2.5, 2);  % 2.5th percentile (M x 1)
error_upper = prctile(id_errors, 97.5, 2); % 97.5th percentile (M x 1)
errorbar_lower = error_mean - error_lower; % Error bar (lower side)
errorbar_upper = error_upper - error_mean; % Error bar (upper side)

figure;
loglog(N_values, 0.8*N_values.^(-0.5), 'Color', red, 'LineStyle', '--', 'LineWidth', 2);
hold on;
errorbar(N_values, error_mean, errorbar_lower, errorbar_upper, 'LineWidth', 2, 'CapSize', 10, 'Color', blue);
xlabel("Number of samples at each frequency $N$", 'Interpreter', 'latex');
ylabel("$\mathcal{H}_\infty$ identification error", 'Interpreter', 'latex');
exportgraphics(gcf, 'results/N_dep/N_dep.pdf', 'ContentType', 'vector');

%%
figure;
% Hinf vs M
num_exp = 20;
num_M = 1;
tau = 7;
M_values = floor(logspace(1, log10(55), num_M));
N = 30; % the number of frequency points
delta = 0.1; % the margin for frequency points
eta = 0.5; % the bound on noise
id_errors = zeros(num_M, num_exp);

for roll = 1:num_exp
    counter = 0;
    for M = M_values
        counter = counter + 1;
        fprintf("exp = %d, M = %.2d \n", roll, M);
        arc_dist = (pi - 2*delta) / M; % the arc distance between frequency points
        theta = (delta + arc_dist/2):arc_dist:(pi - delta - arc_dist/2);
        z = exp(1i.*theta);
        zbar = conj(z);
        wbar = freqresp(Gbar, z); wbar = squeeze(wbar(1,1,:));
        Wbar = repmat(wbar, 1, N);
        V = eta*(2*rand(M, N)-1) + 1i*eta*(2*rand(M, N)-1);
        Wtil = Wbar + V;
        tic;
        what = LNNM(z, Wtil, tau);
        elapsed_time = toc;
        fprintf('LNNM took %.2f seconds.\n', elapsed_time);
        print("_______________________________")
        id_errors(counter, roll) = max(abs(what - wbar));
    end
end

filename = sprintf("results/M_dep/N%d_eta%d_delta%d_num_exp%d_num_M%d_tau_%d.mat", N, eta, delta, num_exp, num_M, tau);
save(filename, 'id_errors', 'M_values');


%%

load("results/M_dep/N30_eta5.000000e-01_delta1.000000e-01_num_exp20_num_M10_tau_7.mat")

error_mean = mean(id_errors, 2);
error_lower = prctile(id_errors, 2.5, 2);  % 2.5th percentile (M x 1)
error_upper = prctile(id_errors, 97.5, 2); % 97.5th percentile (M x 1)
errorbar_lower = error_mean - error_lower;    % Error bar (lower side)
errorbar_upper = error_upper - error_mean;    % Error bar (upper side)

figure;
% loglog(M_values, 0.03*M_values.^(0.5), 'r--', 'LineWidth', 2);
hold on;
errorbar(M_values, error_mean,  errorbar_lower, errorbar_upper, 'LineWidth', 2, 'CapSize', 10, 'Color', blue);
xlabel("Number of frequency points $M$", 'Interpreter', 'latex');
ylabel("$\mathcal{H}_\infty$ identification error", 'Interpreter', 'latex');
exportgraphics(gcf, 'results/M_dep/M_dep.pdf', 'ContentType', 'vector');


%%

figure;
% Hinf vs eta
num_exp = 20;
N = 30;
M = 32;
eta_values = 0.25:0.25:2;
num_eta = length(eta_values);
tau = 20;
id_errors = zeros(num_eta, num_exp);
delta = 0.1; % the margin for frequency points

counter = 0;
for roll = 1:num_exp
    counter = 0;
    for eta = eta_values
        counter = counter + 1;
        fprintf("exp = %d, eta = %.2d \n", roll, eta);

        Wbar = repmat(wbar, 1, N);
        V = eta*(2*rand(M, N)-1) + 1i*eta*(2*rand(M, N)-1);
        Wtil = Wbar + V;
        what = LNNM(z, Wtil, tau);

        id_errors(counter, roll) = norm(what - wbar, 'inf');
    end
end

filename = sprintf("results/eta_dep/N%d_M%d_delta%d_num_exp%d_num_eta%d_tau_%d.mat", N, M, delta, num_exp, num_eta, tau);
save(filename, 'id_errors', 'eta_values');

%%

load("results/eta_dep/N30_M32_delta1.000000e-01_num_exp20_num_eta8_tau_20.mat")

error_mean = mean(id_errors, 2);
error_lower = prctile(id_errors, 2.5, 2);  % 2.5th percentile (M x 1)
error_upper = prctile(id_errors, 97.5, 2); % 97.5th percentile (M x 1)
errorbar_lower = error_mean - error_lower;    % Error bar (lower side)
errorbar_upper = error_upper - error_mean;    % Error bar (upper side)

figure;
% plot(eta_values, 0.3*eta_values, 'Color', red, 'LineStyle', '--', 'LineWidth', 2);
hold on;
errorbar(eta_values, error_mean,  errorbar_lower, errorbar_upper, 'LineWidth', 2, 'CapSize', 10, 'Color', blue);
xlabel("Noise level $\eta$", 'Interpreter', 'latex');
ylabel("$\mathcal{H}_\infty$ identification error", 'Interpreter', 'latex');
exportgraphics(gcf, 'results/eta_dep/eta_dep.pdf', 'ContentType', 'vector');



