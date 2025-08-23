
clear
clc
close all

set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, 'DefaultLineLineWidth', 2);

rng(6);

%%
delta = 0.1; % the margin for frequency points
M_values = 5:5:100;
num_M = length(M_values);
int_app = zeros(num_M, 1);
counter = 0;

for M = M_values
    counter = counter + 1;
    arc_dist = (pi - 2*delta) / M; % the arc distance between frequency points
    theta = (delta + arc_dist/2):arc_dist:(pi - delta - arc_dist/2);
    z = exp(1i.*theta);
    zbar = conj(z);
    sum = 0;
    for i = 1:M
        for j = 1:M
            sum = sum + 1 / abs(zbar(i) - z(j))^2;
        end
    end
    int_app(counter) = sum / M^2;
end

int = -2*log(sin(delta)) / (pi - 2*delta)^2;

figure
hold on;
plot(M_values, int_app);
yline(int, '--r');