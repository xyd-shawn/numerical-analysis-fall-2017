function y = my_func_1(x)
% Runge function
y = 1 ./ (1 + 25 * x .^2);
n = length(y);
if size(y, 2) == n
    y = y';
end
end