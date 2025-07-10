clc, clearvars, close all;

dataFile = load('data.mat'); 
data = dataFile.time_series;
frequency = 1000;
lenData = length(data);
normalizedData = (data - min(data)) / (max(data) - min(data)); 
time = (1:lenData) / frequency;

figure; 
plot(time, normalizedData, Color = 'Red'); 
xlabel("Time(Second)")
ylabel("Intensity of the earthquake")
grid on;
%% Part 3
coeffs = [];
centerData = normalizedData(lenData / 3 + 1 : lenData * 2 / 3);
N = length(centerData);
coeffs = [coeffs, abs(calculateFourierCoeff(centerData, 1:10))];

figure; 
plot(coeffs, Color = 'Green', LineWidth = 1, Marker = 'X', MarkerEdgeColor='Black', MarkerSize=10); 
xlabel("Frequency")
ylabel("Fourier Series Coefficient")
grid on;
%% Functions
function fourierCoeff = calculateFourierCoeff(data, k)
    N = length(data);
    n = (1:N)';
    sigma = sum((exp(-1 * 1i * 2 * pi * n * k / N) .* data));
    fourierCoeff = sigma / N;
end
