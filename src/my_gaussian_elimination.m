function x = my_gaussian_elimination(A, b, tol)
    % 实现高斯消去法求解线性方程组
    % 输入: A为系数矩阵，b为右端项，tol判断系数矩阵能否继续消去的阈值
    % 输出: x为线性方程组的解
    n = length(b);
    B = [A b];
    x = zeros(n, 1);
    for k = 1:(n-1)
        for i = (k+1):n
            if abs(B(k, k)) < tol
                disp('Method failed!');
                return;
            end
            pik = B(i, k) / B(k, k);
            B(i, (k+1):end) = B(i, (k+1):end) - pik * B(k, (k+1):end);    % 化为上三角矩阵
        end
    end
    if abs(B(n, n)) < tol
        disp('Method failed!');
        return;
    end
    x(n) = B(n, n+1) / B(n, n);
    for i = (n-1):(-1):1
        x(i) = (B(i, n+1) - B(i, (i+1):n) * x((i+1):n)) / B(i, i);
    end
end
