function A = triangular_function(t, T, Shift)
    A = (abs(t - Shift) < T) .* (1 - abs(t - Shift) / T);
end