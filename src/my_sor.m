function [x, niters] = my_sor(A, b, x0, w, max_iter, tol)
    % 实现SOR方法求解线性方程组
    % 输入: A为系数矩阵，b为右端项，x0为初始向量, w为权重，max_iter为最大迭代次数，tol为判断迭代终止的阈值
    % 输出: x为线性方程组的解
    n = length(x0);
    x1 = zeros(n, 1);
    iter = 1;
    while iter <= max_iter
        niters = iter;
        % x_k=w*D^(-1)*(L*x_k+U*x_(k-1)+b)+(1-w)*x_(k-1)
        for k = 1:n
            temp = (b(k) - A(k, 1:(k-1)) * x1(1:(k-1)) - A(k, (k+1):n) * x0((k+1):n)) / A(k, k);
            x1(k) = w * temp + (1 - w) * x0(k);
        end
        if norm(x1 - x0, 2) < tol
            break;
        end
        x0 = x1;
        iter = iter + 1;
    end
    x = x1;
    if iter > max_iter
        fprintf('Method does not converge in %d iterations!\n', max_iter);
    end
end
