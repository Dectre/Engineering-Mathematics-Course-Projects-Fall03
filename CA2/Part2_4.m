clc,clearvars,close all;
alpha = 0.1;
W = 1;
H = 1;
Nx = 100;
Ny = 100;
dx = W/(Nx-1);
dy = H/(Ny-1);
dt = 0.0001;
totalTime = 1;
T_temp = 25 * ones(Ny, Nx);
T_temp(1, 1:end-1) = 75;
T_temp(end, 1:end-1) = 130;
T_temp(1:end-1, 1) = 110;
T_temp(1:end-1, end) = 45;
T = zeros(Ny, Nx, 4);
T(:, :, 1) = T_temp;
n = totalTime/dt;
for k = 1:n
    T_old = T_temp;
    for i = 2:Ny-1
        for j = 2:Nx-1
            d2T_dx2 = (T_old(i, j+1) - 2*T_old(i, j) + T_old(i, j-1)) / dx^2;
            d2T_dy2 = (T_old(i+1, j) - 2*T_old(i, j) + T_old(i-1, j)) / dy^2;
            T_temp(i, j) = T_old(i, j) + alpha * dt * (d2T_dx2 + d2T_dy2);
        end
    end
    T_temp(1, 1:end-1) = 75;
    T_temp(end, 1:end-1) = 130;
    T_temp(1:end-1, 1) = 110;
    T_temp(1:end-1, end) = 45;
    if k == 0.25/dt
        T(:, :, 2) = T_temp;
    elseif k == 0.5/dt
        T(:, :, 3) = T_temp;
    elseif k == 1/dt
        T(:, :, 4) = T_temp;
    end
end
timeLabels = {'t = 0', 't = 0.25', 't = 0.5', 't = 1'};
X = linspace(0, W, Nx);
Y = linspace(0, H, Ny);
figure;
for i = 1:4
    subplot(2, 2, i);
    pcolor(X, Y, T(:, :, i));
    shading interp;
    colorbar;
    title(timeLabels{i});
end
sgtitle('Transient-State Heat Equation using FDM');