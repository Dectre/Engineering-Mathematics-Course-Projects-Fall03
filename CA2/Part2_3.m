clc,clearvars, close all;
W = 1;
H = 1;
Nx = 100;
Ny = 100;
X = linspace(0, W, Nx);
Y = linspace(0, H, Ny);
T = 25 * ones(Ny, Nx);
T(1, 1:end-1) = 75;
T(end, 1:end-1) = 130;
T(1:end-1, 1) = 110;
T(1:end-1, end) = 45;
dx = 1;
n = 10000;
tol = 5e-4;
for k = 1:n
    T_old = T;
    for i = 2:Nx-1
        for j = 2:Ny-1
            T(i,j) = 0.25 * (T_old(i+1,j) + T_old(i-1,j) + dx ^ 2 * (T_old(i,j+1) + T_old(i,j-1)));
        end
    end
    if abs(T(2:end-1) - T_old(2:end-1)) < tol
        break;
    end
end
figure;
pcolor(X,Y,T);
shading interp;
colorbar;
title("Steady-State Heat Equation using FDM");
