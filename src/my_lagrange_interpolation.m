function yy = my_lagrange_interpolation(f, x, xx)
% 实现Lagrange多项式插值法
% 输入:f为被插函数,x是插值基点,xx是所要求得函数值点
% 输出:yy为所要求的函数值,与xx一一对应
n = length(x);
nn = length(xx);
if size(xx, 2) == nn
    xx = xx';
end
yy = zeros(nn, 1);
for i = 1:n
    li = ones(nn, 1);
    for j = 1:n
        if j ~= i
           li =  (li .* (xx - x(j))) / (x(i) - x(j));
        end
    end
    yy = yy + li * f(x(i));
end
end