function [p, yy] = my_newton_interpolation(f, x, xx)
% 实现Newton多项式插值法
% 输入: f为被插函数，x是插值基点，xx是要求取的插值函数值点
% 输出: p是均差系数(由低次到高次)，yy是要求取的插值函数值
n = length(x);
mat = zeros(n);
if size(x, 2) == n
    x = x';
end
mat(:, 1) = f(x);
for i = 2:n
    mat(i:n, i) = (mat(i:n, i-1) - mat((i-1):(n-1), i-1)) ./ (x(i:n) - x(1:(n-i+1)));    % 构造均差表
end
p = diag(mat);
nn = length(xx);
if size(xx, 2) == nn
    xx = xx';
end
yy = zeros(nn, 1);
for i = n:(-1):2
    yy = (yy + p(i)) .* (xx - x(i - 1));
end
yy = yy + p(1);
end
