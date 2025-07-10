function F = general_fourier_transform(f, omega)
    F = arrayfun(@(w) integral(@(t) f(t) .* exp(-1j * w * t), -10, 10), omega);
end
