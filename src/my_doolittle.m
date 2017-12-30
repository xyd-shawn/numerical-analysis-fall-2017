function [B] = my_doolittle(A)
  % 实现矩阵的Doolittle分解
  B = A;
  n = size(B, 1);
  B(2:end, 1) = B(2:end, 1) / B(1, 1);
  for k = 2:n
    for j = k:n
      B(k, j) = B(k, j) - B(k, 1:(k - 1)) * B(1:(k - 1), j);
    end
    for j = (k + 1):n
      B(j, k) = (B(j, k) - B(j, 1:(k - 1)) * B(1:(k - 1), k)) / B(k, k);
    end
  end
end
