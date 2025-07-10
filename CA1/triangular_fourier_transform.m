function F = triangular_fourier_transform(T, omega)
    A = @(t) triangular_function(t, T, 0);
    F = arrayfun(@(w) integral(@(t) A(t) .* exp(-1j * w * t), -T, T), omega);
end
