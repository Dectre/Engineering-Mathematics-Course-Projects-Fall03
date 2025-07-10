function F = approximate_fourier_transform(f, t, w_vals)
    idx = f(t) ~= 0;
    T = t(2) - t(1);
    tNew = t(idx):T:t(end);
    FTriangular = triangular_fourier_transform(T, w_vals);
    F = arrayfun(@(w, i) sum(exp(-1j * w * tNew) .* f(tNew) .* FTriangular(i)), w_vals, 1:length(w_vals));
end

