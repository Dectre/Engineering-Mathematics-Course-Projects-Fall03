clc,clearvars, close all;
%%
f = @(t) 3*((abs(t/2)<=0.5)+(abs(t)<=0.5));
omega=linspace(-50,50,1000);
F_f = general_fourier_transform(f, omega);

figure;
plot(omega, abs(F_f), Color = "blue", LineWidth = 2);
xlabel('W');
ylabel('F');
%%
g = @(t) exp((-1i)*2*pi*3.*t) .* f(t);
F_g = general_fourier_transform(g,omega);
figure;
plot(omega, abs(F_g), Color = "magenta", LineWidth = 2);
xlabel('W');
ylabel('F');
%%
figure;
subplot(2, 1, 1);
plot(omega, abs(F_f), Color = "blue", LineWidth = 2);
title('F{f(t)}');

subplot(2, 1, 2);
plot(omega, abs(F_g), Color = "magenta", LineWidth = 2);
title('F{g(t)}');




