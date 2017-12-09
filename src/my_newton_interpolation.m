function [p, yy] = my_newton_interpolation(f, x, xx)
% ʵ��Newton����ʽ��ֵ��
% ����:f�Ǳ��庯��,x�ǲ�ֵ����,xx��Ҫ��ȡ�Ĳ�ֵ����ֵ��
% ���:p�Ǿ���ϵ��(�ɵʹε��ߴ�),yy��Ҫ��ȡ�Ĳ�ֵ����ֵ
n = length(x);
mat = zeros(n);
if size(x, 2) == n
    x = x';
end
mat(:, 1) = f(x);
for i = 2:n
    mat(i:n, i) = (mat(i:n, i-1) - mat((i-1):(n-1), i-1)) ./ (x(i:n) - x(1:(n-i+1)));    % ��������
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