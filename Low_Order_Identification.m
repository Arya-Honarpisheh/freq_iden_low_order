
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

Gbar= tf([0.12 0.18], [1 -1.4 1.443 -1.123 0.7729], 1); % true system

M = 32; % the number of frequency points
delta = 0.1; % the margin for frequency points
eta = 2; % the bound on noise
% Here we use logspace theta to have better resolution at higher
% frequencies
theta = logspace(-2, 0.49, M);
z = exp(1i.*theta);
zbar = conj(z);

wbar = freqresp(Gbar, z); wbar = squeeze(wbar(1,1,:));
frdatabar = idfrd(wbar, theta, 1);
[magbar, phasebar] = bode(frdatabar);
magbar = squeeze(magbar);
phasebar = squeeze(phasebar);

%%

% Hinf identification error comparison
num_exp = 20;
N = 30;
% For this experiemnt we just choose some tau by trial and error. It is not
% optimal and one can find the optimal one by cross validation.
tau = 20;
errors_reg = zeros(num_exp);
errors_noreg = zeros(num_exp);

for k  = 1:num_exp

    fprintf("exp = %d \n", k);

    Wbar = repmat(wbar, 1, N);
    V = eta*(2*rand(M, N)-1) + 1i*eta*(2*rand(M, N)-1);
    Wtil = Wbar + V;

    wtil_mean = mean(Wtil, 2);

    what = LNNM(z, Wtil, tau);

    errors_reg(k) = norm(what - wbar, 'inf');
    errors_noreg(k) = norm(wtil_mean - wbar, 'inf');

end

filename = sprintf("results/comparison/M%d_N%d_eta%d_num_exp%d_tau_%d.mat", M, N, eta, num_exp, tau);
save(filename, 'errors_reg', 'errors_noreg');

%%

% comparison of identification error
load("results/comparison/M32_N30_eta2_num_exp20_tau_20.mat")
num_exp = 20;
figure;
boxplot([errors_reg(:,1), errors_noreg(:,1)], 'Labels', ["LNNM", "Averaging"], 'Colors', [green; red]);
set(findobj(gca, 'Type', 'Line'), 'LineWidth', 2); % Set all lines thicker
% ylim([0.3 1.4]);
ylabel("$\mathcal{H}_\infty$ identification error", 'Interpreter', 'latex');
exportgraphics(gcf, 'results/comparison/comparsion_finite_sample.pdf', 'ContentType', 'vector');

%%

% comparison of bode plots
rng(1);
N = 30;

Wbar = repmat(wbar, 1, N);
V = eta*(2*rand(M, N)-1) + 1i*eta*(2*rand(M, N)-1);
Wtil = Wbar + V;

wtil_mean = mean(Wtil, 2);
frdatahat_noreg = idfrd(wtil_mean, theta, 1);
[maghat_noreg, phasehat_noreg] = bode(frdatahat_noreg);
maghat_noreg = squeeze(maghat_noreg);
phasehat_noreg = squeeze(phasehat_noreg);

tau = 20;
what = LNNM(z, Wtil, tau);
frdatahat= idfrd(what, theta, 1);
[maghat, phasehat] = bode(frdatahat);
maghat = squeeze(maghat);
phasehat = squeeze(phasehat);

figure;
subplot(2,1,1);
hold on;
plot(theta, 20*log10(magbar), 'Color', blue, 'LineStyle', '-');
plot(theta, 20*log10(maghat_noreg), 'Color', red, 'LineStyle', '--');
plot(theta, 20*log10(maghat), 'Color', green, 'LineStyle', '--');
% legend(["Baseline", "Averaging", "LNNM"], "Location", "eastoutside", 'FontSize', 8);
xlim([0 pi]);
set(gca, 'XScale', 'log');
ylabel('Magnitude (dB)');
% title("Bode Magnitude Plot");
grid on;
subplot(2,1,2);
hold on;
plot(theta, phasebar, 'Color', blue, 'LineStyle', '-');
plot(theta, phasehat_noreg, 'Color', red, 'LineStyle', '--');
plot(theta, phasehat, 'Color', green, 'LineStyle', '--');
% legend(["Baseline", "Averaging", "LNNM"], "Location", "eastoutside", 'FontSize', 8);
xlim([0 pi]);
set(gca, 'XScale', 'log');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
% title('Bode Phase Plot');
grid on;

exportgraphics(gcf, 'results/comparison/comparsion_bode.pdf', 'ContentType', 'vector');

%%

% singular values of Loewner matrix
figure;
hold on
bar_handle = bar([svd(loewner(z, wbar)), svd(loewner(z, wtil_mean)), svd(loewner(z, what))]);
bar_handle(1).FaceColor = blue; bar_handle(2).FaceColor = red; bar_handle(3).FaceColor = green;
for k = 1:3
    bar_handle(k).BarWidth = 1;
end
xlim([0, 10]);
legend(["True System", "Averaging", "LNNM"]);

exportgraphics(gcf, 'results/comparison/comparsion_svd.pdf', 'ContentType', 'vector');
