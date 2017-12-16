function [lam, x] = my_power_method(A, v0, max_iter, tol)
    % 实现乘幂法求解矩阵主特征值和相应的特征向量
    % 输入: A为矩阵，v0为初始的向量，max_iter为最大迭代次数，tol为判断终止条件的阈值
    % 输出: lam为主特征值，x为对应的主特征向量
    n = length(v0);
    if size(v0, 2) == n
        v0 = v0';
    end
    m0 = 0;
    for i = 1:n
        if abs(v0(i)) > abs(m0)
            m0 = v0(i);
        end
    end
    v0 = v0 / m0;
    iter = 1;
    while iter <= max_iter
        v1 = A * v0;
        m1 = 0;
        for i = 1:n
            if abs(v1(i)) > abs(m1)
                m1 = v1(i);
            end
        end
        v1 = v1 / m1;
        if abs(m1 - m0) <= tol
            break
        end
        m0 = m1;
        v0 = v1;
        iter = iter + 1;
    end
    if iter > max_iter
        fprintf('Method does not converge in %d iterations!\n', max_iter);
    end
    lam = m1;
    x = v1;
end
