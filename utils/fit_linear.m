function [coefs] = fit_linear(x, y)
%FIT_LINEAR 
[n, ~] = size(x);
A = [x ones(n, 1)];
b = y;
[coefs, ~] = linsolve(A, b);
end
