function [coefs, residual] = fit_quadratic(x, y)
%FIT_QUADRATIC Fits convex quadratic function for given 3 points x, y
%   x - points, y - values, coefs = [a, b, c]: y = f(x) = a*x*x + b*x + c
% x = [1; 2; 3];
% y = [-1; -4; -9];
[n, ~] = size(x);
A = [x.*x x ones(n, 1)];
b = y;
[coefs, ~] = linsolve(A, b);
residual = norm(A * coefs - b) / norm(b);
%{
if r < 1e-8 || coefs(1) < 0
    A = [x(1) 1; x(end) 1];
    b = [y(1); y(end)];
    [lin_coefs, ~] = linsolve(A, b);
    coefs = [0; lin_coefs];
end
%}

end

