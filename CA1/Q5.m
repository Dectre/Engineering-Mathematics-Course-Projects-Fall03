clc,clearvars,close all;
omega=linspace(-50,50,1000);
t = linspace(-10,10,200);
%% f
w = @(t) abs(t) <= 0.5;
f = @(t) w(t/5).*(1-cos(2*pi*t/5));
Ff = @(omega) 5 * (sinc(5/(2*pi) * omega) - 1/2 * sinc(5/(2*pi)*(omega - 2*pi/5)) - 1/2 * sinc(5/(2*pi)*(omega + 2*pi/5)));

FfTheory = Ff(omega);
FfGeneral=general_fourier_transform(f , omega);
FfApproximate=approximate_fourier_transform(f,t,omega);

figure;
subplot(3,1,1)
plot(omega,abs(FfTheory), Color = 'Red', LineWidth = 1)

subplot(3,1,2)
plot(omega,abs(FfGeneral), Color = 'Green', LineWidth = 1)

subplot(3,1,3)
plot(omega,abs(FfApproximate), Color = 'Blue', LineWidth = 1)
%% g 
Fg = @(omega) sinc(omega / 2 / (pi)) .^2 .* (1 + exp(-1j * omega * 2) + exp(-1j * omega * 4) + exp(-1j * omega * 6));
g = @(t) triangular_function(t,1,0)+triangular_function(t,1,2)+triangular_function(t,1,4)+triangular_function(t,1,6);
t = linspace(-10,10,200);
FgTheory = abs(Fg(omega));
FgGeneral=abs(general_fourier_transform(g , omega));
FgApproximate=abs(approximate_fourier_transform(g,t,omega));

figure;
subplot(3,1,1)
plot(omega,abs(FgTheory), Color = 'Red', LineWidth = 1);

subplot(3,1,2)
plot(omega,abs(FgGeneral), Color = 'Green', LineWidth = 1);

subplot(3,1,3)
plot(omega,abs(FgApproximate), Color = 'Blue', LineWidth = 1);