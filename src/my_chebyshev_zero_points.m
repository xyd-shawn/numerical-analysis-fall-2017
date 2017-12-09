function x = my_chebyshev_zero_points(a, b, n)
% 生成chebyshev多项式的零点,并变换到区间[a, b]上
% 输入:[a, b]是区间,n是零点的数目
% 输出:x是变换到区间[a, b]上的基点
t = 1:2:(2*n-1);
res = cos(t * pi / (2 * n));
if size(res, 2) == n
    res = res';
end
x = (b + a) / 2 + ((b - a) / 2) * res;
end