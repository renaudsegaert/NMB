function [c, kappa] = interpolate(x, f)
%INTERPOLATE Calculates an interpolation using the Chebyshev polynomes.

n = length(x);
M = zeros(n);
b = zeros(n, 1);

for col=1:n
    M(:, col) = chebyshevT(col-1, x)';
    b(col) = f(x(col));
end

c = linsolve(M, b);
kappa = cond(M);

% M, b, c

end

