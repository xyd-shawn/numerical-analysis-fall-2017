function yy = my_lagrange_interpolation(f, x, xx)
% ʵ��Lagrange����ʽ��ֵ��
% ����:fΪ���庯��,x�ǲ�ֵ����,xx����Ҫ��ú���ֵ��
% ���:yyΪ��Ҫ��ĺ���ֵ,��xxһһ��Ӧ
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