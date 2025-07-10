%% Load
clc, clearvars, close all;
song = audioread('song.wav'); 

%% Process
fs = 44100;
segmentSize = 0.4 * fs;
numOfSegments = floor(length(song) / segmentSize);
segments = reshape(song(1:numOfSegments * segmentSize), segmentSize, numOfSegments);
noteFrequencies = [523.25, 587.33, 659.26, 698.46, 784.00];
notes = {'C', 'D', 'E', 'F', 'G'};
recognizedNotes = cell(numOfSegments, 1);
coeffs = zeros(length(noteFrequencies), numOfSegments);

for i = 1:numOfSegments
    segment = segments(:,i);
    meanAmplitude = mean(abs(segment));
    for j = 1:length(noteFrequencies)
        k = round(noteFrequencies(j) * segmentSize / fs);
        coeff = calculateFourierCoeff(segment, k);
        coeffs(j, i) = abs(coeff);
        if abs(coeff) > 0.5 * meanAmplitude
            recognizedNotes{i} = notes{j};
        end
    end
end
%% Plot
colors = [1 0 0;
          1 1 0;
          0 1 0;
          0 0 1;
          1 0 1];
for page = 1:2
    figure(page);
    for plotIdx = 1:16
        subplot(4, 4, plotIdx);
        segmentIdx = (page - 1) * 16 + plotIdx;
        if segmentIdx <= numOfSegments
            b = bar(coeffs(:, segmentIdx));
            b.FaceColor = 'flat';
            for barIdx = 1:length(notes)
                b.CData(barIdx, :) = colors(mod(barIdx - 1, size(colors, 1)) + 1, :);
            end
            set(gca, 'XTickLabel', notes);
        end
    end
end

%% Display recognized notes
disp(recognizedNotes');

%% Functions 
function fourierCoeff = calculateFourierCoeff(data, k)
    N = length(data);
    n = (1:N)';
    sigma = sum((exp(-1 * 1i * 2 * pi * (n-1) * k / N) .* data));
    fourierCoeff = sigma / N;
end