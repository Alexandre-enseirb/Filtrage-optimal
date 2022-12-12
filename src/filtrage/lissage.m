function X_smooth = lissage(X_hat, P, phi)
    N = length(Y);
    X_smooth = zeros(size(X_hat));
    X_smooth(:, N) = X_hat(:, N);
    P_smooth = P.aposter(:, :, N);
    for k = N-1:-1:1
        A = P.aposter(:, :, k) * phi.' / P.apriori(:, :, k+1);
        P_smooth = P.aposter(:, :, k) + A*(P_smooth - P.apriori(:, :, k+1))*A.';
        X_smooth(:, k) = X_hat(:, k) + A*(X_smooth(:, k+1) - phi*X_hat(:, k));
    end
end
