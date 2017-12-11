function x = my_chase(A, b)
% ʵ��׷�Ϸ�������Խ����Է�����
% ����:AΪ���ԽǾ���,bΪ�Ҷ�ֵ
% ���:xΪ���Է�����Ľ�
n = length(b);
u = zeros(n, 1);
l = zeros(n, 1);
% LU�ֽ�
u(1) = A(1, 1);
for i = 2:n
    l(i) = A(i, i-1) / u(i-1);
    u(i) = A(i, i) - A(i-1, i) * l(i);
end
y = zeros(n, 1);
y(1) = b(1);
for i = 2:n
    y(i) = b(i) - l(i) * y(i-1);
end
x = zeros(n, 1);
x(n) = y(n) / u(n);
for i = (n-1):-1:1
    x(i) = (y(i) - A(i, i+1) * x(i+1)) / u(i);
end
end