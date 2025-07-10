clc, clearvars, close all;
%%
f1 = @(t) sinc(t).^2;
f2 = @(t) exp(-3 * abs(t));
w = @(t) abs(t) <= 0.5;

omega = linspace(-50, 50, 1000);
F1 = general_fourier_transform((@(t) f1(t) .* w(t)), omega);
F2 = general_fourier_transform((@(t) f2(t) .* w(t)), omega);

figure;
subplot(2, 1, 1);

plot(omega, abs(F1), Color = 'Red', LineWidth= 2);
%% 
title('F\{f1(t)\}');

subplot(2, 1, 2);
plot(omega, abs(F2), Color = 'Magenta', LineWidth= 2);
title('F\{f2(t)\}');
%%
F1 = general_fourier_transform(f1, omega);
F2 = general_fourier_transform(f2, omega);
W = general_fourier_transform(w, omega);

F1Convolution = conv(F1, W,'same') * (omega(2) - omega(1));
F2Convolution = conv(F2, W,'same') * (omega(2) - omega(1));
figure;
subplot(2, 1, 1);
plot(omega, abs(F1Convolution), Color = 'Red', LineWidth= 2);
title('Convolution F1');

subplot(2, 1, 2);
plot(omega, abs(F2Convolution), Color = 'Magenta', LineWidth= 2);
title('Convolution F2');

%%
figure;
subplot(2,1,1)
plot(omega, abs(F1Convolution / (2 * pi)), Color = 'Red', LineWidth= 2);
title('F\{f1(t)\} using convolution');

subplot(2, 1, 2);
plot(omega, abs(F2Convolution / (2 * pi)), Color = 'Magenta', LineWidth= 2);
title('F\{f1(t)\} using convolution');
