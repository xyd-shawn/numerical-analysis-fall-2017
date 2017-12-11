function [m, y, yy] = my_cubic_spline_interpolation_1(f, x, xx)
% 实现三次自然样条插值函数
% 输入:f为被插函数,x为插值基点,xx为要求的插值函数值点,要求x升序排列
% 输出:m为插值基点对应的一阶导数,y为插值基点的函数值,yy为插值函数对应xx的函数值
n = length(x);
if size(x, 2) == n
    x = x';
end
y = f(x);
h = x(2:n) - x(1:(n-1));
df = (y(2:n) - y(1:(n-1))) ./ h;
mu = ones(n, 1);
la = ones(n, 1);
d = zeros(n, 1);
mu(2:(n-1)) = h(2:(n-1)) ./ (h(1:(n-2)) + h(2:(n-1)));
la(2:(n-1)) = 1 - mu(2:(n-1));
d(2:(n-1)) = 3 * (df(1:(n-2)) .* mu(2:(n-1)) + df(2:(n-1)) .* la(2:(n-1)));
d(1) = 3 * df(1);
d(n) = 3 * df(n-1);
A = zeros(n);    % 构造三对角系数矩阵
for i = 1:(n-1)
    A(i, i+1) = la(i);
    A(i+1, i) = mu(i+1);
    A(i, i) = 2;
end
A(n, n) = 2;
m = my_chase(A, d);
nn = length(xx);
if size(xx, 2) == nn
    xx = xx';
end
yy = zeros(nn, 1);
for i = 1:nn
    ind = (xx(i) >= x);
    pos = length(x(ind));
    if xx(i) == x(n)
        pos = pos -1;
    end
    temp1 = (xx(i) - x(pos)) / (x(pos+1) - x(pos));
    temp2 = (xx(i) - x(pos+1)) / (x(pos) - x(pos+1));
    yy(i) = yy(i) + y(pos) * (1 + 2 * temp1) * temp2 * temp2;
    yy(i) = yy(i) + y(pos+1) * (1 + 2 * temp2) * temp1 * temp1;
    yy(i) = yy(i) + m(pos) * (xx(i) - x(pos)) * temp2 * temp2;
    yy(i) = yy(i) + m(pos+1) * (xx(i) - x(pos+1)) * temp1 * temp1;
end
end