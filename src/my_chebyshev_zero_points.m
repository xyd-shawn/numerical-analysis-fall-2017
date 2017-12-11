function x = my_chebyshev_zero_points(a, b, n)
% ����chebyshev����ʽ�����,���任������[a, b]��
% ����:[a, b]������,n��������Ŀ
% ���:x�Ǳ任������[a, b]�ϵĻ���
t = 1:2:(2*n-1);
res = cos(t * pi / (2 * n));
if size(res, 2) == n
    res = res';
end
x = (b + a) / 2 + ((b - a) / 2) * res;
end